// =============================================================================
//  CADET
//  
//  Copyright © 2008-present: The CADET-Core Authors
//            Please see the AUTHORS.md file.
//  
//  All rights reserved. This program and the accompanying materials
//  are made available under the terms of the GNU Public License v3.0 (or, at
//  your option, any later version) which accompanies this distribution, and
//  is available at http://www.gnu.org/licenses/gpl.html
// =============================================================================

#include "model/binding/BindingModelBase.hpp"
#include "model/ExternalFunctionSupport.hpp"
#include "model/ModelUtils.hpp"
#include "cadet/Exceptions.hpp"
#include "model/Parameters.hpp"
#include "LocalVector.hpp"
#include "SimulationTypes.hpp"

#include <functional>
#include <unordered_map>
#include <string>
#include <vector>

/*<codegen>
{
	"name": "MultiComponentLDFFreundlichParamHandler",
	"externalName": "ExtMultiComponentLDFFreundlichParamHandler",
	"parameters":
		[
			{ "type": "ScalarComponentDependentParameter", "varName": "kLDF", "confName": "MCLDFFRL_KLDF"},
			{ "type": "ScalarComponentDependentParameter", "varName": "kF", "confName": "MCLDFFRL_KF"},
			{ "type": "ScalarComponentDependentParameter", "varName": "n", "confName": "MCLDFFRL_EXP"},
			{ "type": "ComponentDependentComponentVectorParameter", "varName": "a", "confName": "MCLDFFRL_A"}
		],
	"constantParameters":
		[
			{ "type": "ScalarParameter", "varName": "tau", "confName": "MCLDFFRL_TAU"}
		]
}
</codegen>*/

namespace cadet
{
namespace model
{

class FieldMultiComponentLDFFreundlichParamHandler : public FieldParamHandlerBase
{
public:
	struct ConstParams
	{
		typename ScalarParameter::storage_t tau;
		typename ComponentDependentComponentVectorParameter::storage_t a;
	};

	struct VariableParams
	{
		util::LocalVector<active> kLDF;
		util::LocalVector<active> kF;
		util::LocalVector<active> n;
		ComponentDependentComponentVectorParameter::storage_t a;
		ScalarParameter::storage_t tau;
	};

	typedef VariableParams params_t;
	typedef ConstBufferedScalar<params_t> ParamsHandle;

	static inline const char* identifier() CADET_NOEXCEPT;

	FieldMultiComponentLDFFreundlichParamHandler() CADET_NOEXCEPT : _tau(&_constParams.tau), _a(&_constParams.a)
	{ }

	inline bool configure(IParameterProvider& paramProvider, unsigned int nComp, unsigned int const* nBoundStates)
	{
		_kLDF.configure("MCLDFFRL_KLDF", paramProvider, nComp, nBoundStates);
		_kF.configure("MCLDFFRL_KF", paramProvider, nComp, nBoundStates);
		_n.configure("MCLDFFRL_EXP", paramProvider, nComp, nBoundStates);
		_a.configure("MCLDFFRL_A", paramProvider, nComp, nBoundStates);
		_tau.configure("MCLDFFRL_TAU", paramProvider, nComp, nBoundStates);

		FieldParamHandlerBase::configure(paramProvider, {"MCLDFFRL_KLDF", "MCLDFFRL_KF", "MCLDFFRL_EXP"}, {"TIME", "AXIAL", "RADIAL", "PARTICLE"});

		return validateConfig(nComp, nBoundStates);
	}

	inline void registerParameters(std::unordered_map<ParameterId, active*>& parameters, UnitOpIdx unitOpIdx, ParticleTypeIdx parTypeIdx, unsigned int nComp, unsigned int const* nBoundStates)
	{
		_kLDF.registerParam("MCLDFFRL_KLDF", parameters, unitOpIdx, parTypeIdx, nComp, nBoundStates);
		_kF.registerParam("MCLDFFRL_KF", parameters, unitOpIdx, parTypeIdx, nComp, nBoundStates);
		_n.registerParam("MCLDFFRL_EXP", parameters, unitOpIdx, parTypeIdx, nComp, nBoundStates);
		_a.registerParam("MCLDFFRL_A", parameters, unitOpIdx, parTypeIdx, nComp, nBoundStates);
		_tau.registerParam("MCLDFFRL_TAU", parameters, unitOpIdx, parTypeIdx, nComp, nBoundStates);
	}

	inline void reserve(unsigned int numElem, unsigned int numSlices, unsigned int nComp, unsigned int const* nBoundStates) \
	{
		_kLDF.reserve(numElem, numSlices, nComp, nBoundStates);
		_kF.reserve(numElem, numSlices, nComp, nBoundStates);
		_n.reserve(numElem, numSlices, nComp, nBoundStates);
	}

	inline ParamsHandle update(double t, unsigned int secIdx, const ColumnPosition& colPos, unsigned int nComp, unsigned int const* nBoundStates, LinearBufferAllocator& workSpace) const
	{
		LOG(Debug) << "MultiComponentLDFFreundlichParamHandler update(" << t << ", " << secIdx << ", "
			<< colPos.axial << "/" << colPos.radial << "/" << colPos.particle << ", "
			<< nComp << ", ...)";
		BufferedScalar<params_t> localParams = workSpace.scalar<params_t>();
		BufferedArray<double> fieldBuffer = workSpace.array<double>(_fields.size());

		evaluateField({t, colPos.axial, colPos.radial, colPos.particle}, static_cast<double*>(fieldBuffer));

		localParams->a = _constParams.a;
		localParams->tau = _constParams.tau;

		_kLDF.prepareCache(localParams->kLDF, workSpace);
		_kLDF.update(cadet::util::dataOfLocalVersion(localParams->kLDF), &fieldBuffer[0], nComp, nBoundStates);
		_kF.prepareCache(localParams->kF, workSpace);
		_kF.update(cadet::util::dataOfLocalVersion(localParams->kF), &fieldBuffer[0], nComp, nBoundStates);
		_n.prepareCache(localParams->n, workSpace);
		_n.update(cadet::util::dataOfLocalVersion(localParams->n), &fieldBuffer[0], nComp, nBoundStates);

		return localParams;
	}

	inline std::tuple<ParamsHandle, ParamsHandle> updateTimeDerivative(double t, unsigned int secIdx, const ColumnPosition& colPos, unsigned int nComp, unsigned int const* nBoundStates, LinearBufferAllocator& workSpace) const
	{
		BufferedScalar<params_t> localParams = workSpace.scalar<params_t>();
		BufferedScalar<params_t> p = workSpace.scalar<params_t>();

		BufferedArray<double> fieldBuffer = workSpace.array<double>(_fields.size());
		BufferedArray<double> fieldDerivBuffer = workSpace.array<double>(_fields.size());

		evaluateField({t, colPos.axial, colPos.radial, colPos.particle}, static_cast<double*>(fieldBuffer));
		evaluateTimeDerivativeField({t, colPos.axial, colPos.radial, colPos.particle}, static_cast<double*>(fieldDerivBuffer));

		_kLDF.prepareCache(localParams->kLDF, workSpace);
		_kLDF.update(cadet::util::dataOfLocalVersion(localParams->kLDF), &fieldBuffer[0], nComp, nBoundStates);
		_kLDF.prepareCache(p->kLDF, workSpace);
		_kLDF.updateTimeDerivative(cadet::util::dataOfLocalVersion(p->kLDF), &fieldBuffer[0], &fieldDerivBuffer[0], nComp, nBoundStates);

		_kF.prepareCache(localParams->kF, workSpace);
		_kF.update(cadet::util::dataOfLocalVersion(localParams->kF), &fieldBuffer[0], nComp, nBoundStates);
		_kF.prepareCache(p->kF, workSpace);
		_kF.updateTimeDerivative(cadet::util::dataOfLocalVersion(p->kF), &fieldBuffer[0], &fieldDerivBuffer[0], nComp, nBoundStates);

		_n.prepareCache(localParams->n, workSpace);
		_n.update(cadet::util::dataOfLocalVersion(localParams->n), &fieldBuffer[0], nComp, nBoundStates);
		_n.prepareCache(p->n, workSpace);
		_n.updateTimeDerivative(cadet::util::dataOfLocalVersion(p->n), &fieldBuffer[0], &fieldDerivBuffer[0], nComp, nBoundStates);

		return std::make_tuple<ParamsHandle, ParamsHandle>(std::move(localParams), std::move(p));
	}

	inline std::size_t cacheSize(unsigned int nComp, unsigned int totalNumBoundStates, unsigned int const* nBoundStates) const CADET_NOEXCEPT
	{
		return 2 * sizeof(params_t) + alignof(params_t) + 2 * 3 * sizeof(double) + alignof(double) + 2 * (
			_kLDF.additionalDynamicMemory(nComp, totalNumBoundStates, nBoundStates)
			+ _kF.additionalDynamicMemory(nComp, totalNumBoundStates, nBoundStates)
			+ _n.additionalDynamicMemory(nComp, totalNumBoundStates, nBoundStates)
		);
	}

protected:
	inline bool validateConfig(unsigned int nComp, unsigned int const* nBoundStates);
	ConstParams _constParams;
	FieldScalarComponentDependentParameter _kLDF;
	FieldScalarComponentDependentParameter _kF;
	FieldScalarComponentDependentParameter _n;

	ComponentDependentComponentVectorParameter _a;
	ScalarParameter _tau;
};

} // namespace model
} // namespace cadet

namespace cadet
{

namespace model
{

inline const char* MultiComponentLDFFreundlichParamHandler::identifier() CADET_NOEXCEPT { return "MULTI_COMPONENT_LDF_FREUNDLICH"; }

inline bool MultiComponentLDFFreundlichParamHandler::validateConfig(unsigned int nComp, unsigned int const* nBoundStates)
{
	if ((_kLDF.size() != _kF.size()) || (_kLDF.size() != _n.size()) || (_kLDF.size() < nComp))
		throw InvalidParameterException("MCLDFFRL_KLDF, MCLDFFRL_KF, and MCLDFFRL_N have to have the same size");

	if ((_a.slices() != nComp) || (_a.size() != nComp * nComp))
		throw InvalidParameterException("MCLDFFRL_A has to have NCOMP^2 elements");

	return true;
}

inline const char* ExtMultiComponentLDFFreundlichParamHandler::identifier() CADET_NOEXCEPT { return "EXT_MULTI_COMPONENT_LDF_FREUNDLICH"; }

inline bool ExtMultiComponentLDFFreundlichParamHandler::validateConfig(unsigned int nComp, unsigned int const* nBoundStates)
{
	if ((_kLDF.size() != _kF.size()) || (_kLDF.size() != _n.size()) || (_kLDF.size() < nComp))
		throw InvalidParameterException("EXT_MCLDFFRL_KLDF, EXT_MCLDFFRL_KF, and EXT_MCLDFFRL_N have to have the same size");

	if ((_a.slices() != nComp) || (_a.size() != nComp * nComp))
		throw InvalidParameterException("EXT_MCLDFFRL_A has to have NCOMP^2 elements");

	return true;
}

inline const char* FieldMultiComponentLDFFreundlichParamHandler::identifier() CADET_NOEXCEPT { return "FIELD_MULTI_COMPONENT_LDF_FREUNDLICH"; }

inline bool FieldMultiComponentLDFFreundlichParamHandler::validateConfig(unsigned int nComp, unsigned int const* nBoundStates)
{
	if ((_kLDF.size() != _kF.size()) || (_kLDF.size() != _n.size()) || (_kLDF.size() < nComp))
		throw InvalidParameterException("FIELD_MCLDFFRL_KLDF, FIELD_MCLDFFRL_KF, and FIELD_MCLDFFRL_N have to have the same size");

	if ((_a.slices() != nComp) || (_a.size() != nComp * nComp))
		throw InvalidParameterException("FIELD_MCLDFFRL_A has to have NCOMP^2 elements");

	return true;
}


/**
 * @brief Defines the multi component Freundlich binding model with linear driving force
 * @details Implements the multi component linear driving force Freundlich adsorption model: \f[ \begin{align} 
 *              \frac{\mathrm{d}q_i}{\mathrm{d}t} &= k_{ldf,i} \left( q_i^* - q_i \right) \\
 *              q_i^* &= k_{f,i} c_{p,i} \left( \sum_j a_{ij} c_{p,j} \right)^{n_i-1}
 *          \end{align} \f]
 *          where @f$ a_{ii} = 1 @f$ for all i.
 *          Multiple bound states are not supported. 
 *          Components without bound state (i.e., non-binding components) are supported.
 * @tparam ParamHandler_t Type that can add support for external function dependence
 */
template <class ParamHandler_t>
class MultiComponentLDFFreundlichBindingBase : public ParamHandlerBindingModelBase<ParamHandler_t>
{
public:

	MultiComponentLDFFreundlichBindingBase() { }
	virtual ~MultiComponentLDFFreundlichBindingBase() CADET_NOEXCEPT { }

	static const char* identifier() { return ParamHandler_t::identifier(); }

	virtual bool implementsAnalyticJacobian() const CADET_NOEXCEPT { return true; }

	CADET_BINDINGMODELBASE_BOILERPLATE

protected:
	using ParamHandlerBindingModelBase<ParamHandler_t>::_paramHandler;
	using ParamHandlerBindingModelBase<ParamHandler_t>::_reactionQuasistationarity;
	using ParamHandlerBindingModelBase<ParamHandler_t>::_nComp;
	using ParamHandlerBindingModelBase<ParamHandler_t>::_nBoundStates;

	template <typename StateType, typename CpStateType, typename ResidualType, typename ParamType>
	int fluxImpl(double t, unsigned int secIdx, const ColumnPosition& colPos, StateType const* y,
		CpStateType const* yCp, ResidualType* res, LinearBufferAllocator workSpace) const
	{
		using CpStateParamType = typename DoubleActivePromoter<CpStateType, ParamType>::type;
		using StateParamType = typename DoubleActivePromoter<StateType, ParamType>::type;
		typename ParamHandler_t::ParamsHandle const p = _paramHandler.update(t, secIdx, colPos, _nComp, _nBoundStates, workSpace);

		// Fluxes: k_{ldf,i} ( q_i - k_{f,i} c_{p,i} \left( tau + \sum_j a_{ij} c_{p,j} \right)^{n_i-1} )
		unsigned int bndIdx = 0;
		for (int i = 0; i < _nComp; ++i)
		{
			// Skip components without bound states (bound state index bndIdx is not advanced)
			if (_nBoundStates[i] == 0)
				continue;

			active const* const aSlice = p->a[i];
			ResidualType sum = static_cast<double>(p->tau);
			for (int j = 0; j < _nComp; ++j)
			{
				sum += static_cast<ParamType>(aSlice[j]) * yCp[j];
			}

			// Residual
			res[bndIdx] = static_cast<ParamType>(p->kLDF[i]) * ( y[bndIdx] - static_cast<ParamType>(p->kF[i]) * yCp[i] * pow(sum, 1 / static_cast<ParamType>(p->n[i]) - 1.0) );

			// Next bound component
			++bndIdx;
		}

		return 0;
	}

	template <typename RowIterator>
	void jacobianImpl(double t, unsigned int secIdx, const ColumnPosition& colPos, double const* y, double const* yCp, int offsetCp, RowIterator jac, LinearBufferAllocator workSpace) const
	{
		typename ParamHandler_t::ParamsHandle const p = _paramHandler.update(t, secIdx, colPos, _nComp, _nBoundStates, workSpace);

		// Fluxes: k_{ldf,i} ( q_i - k_{f,i} c_{p,i} \left( tau + \sum_j a_{ij} c_{p,j} \right)^{n_i-1} )

		int bndIdx = 0;
		for (int i = 0; i < _nComp; ++i)
		{
			// Skip components without bound states (bound state index bndIdx is not advanced)
			if (_nBoundStates[i] == 0)
				continue;

			active const* const aSlice = p->a[i];

			double cpSum = static_cast<double>(p->tau);
			for (int j = 0; j < _nComp; ++j)
			{
				cpSum += static_cast<double>(aSlice[j]) * yCp[j];
			}

			const double kldf = static_cast<double>(p->kLDF[i]);
			const double kf_ldf = static_cast<double>(p->kF[i]) * kldf;
			const double n = static_cast<double>(p->n[i]);

			// dres_i / dc_{p,i}
			jac[i - bndIdx - offsetCp] = -kf_ldf * pow(cpSum, 1 / n - 1.0);
			// Getting to c_{p,i}: -bndIdx takes us to q_0, another -offsetCp to c_{p,0} and a +i to c_{p,i}.
			//                     This means jac[i - bndIdx - offsetCp] corresponds to c_{p,i}.

			const double factor = kf_ldf * yCp[i] * ( 1 / n - 1.0) * pow(cpSum, 1 / n - 2.0);

			// Fill dres_i / dc_{p,j}
			for (int j = 0; j < _nComp; ++j)
			{
				// dres_i / dc_{p,j}
				jac[j - bndIdx - offsetCp] -= factor * static_cast<double>(aSlice[j]);
				// Getting to c_{p,j}: -bndIdx takes us to q_0, another -offsetCp to c_{p,0} and a +j to c_{p,j}.
				//                     This means jac[j - bndIdx - offsetCp] corresponds to c_{p,j}.
			}

			// Add to dres_i / dq_i
			jac[0] += kldf;

			// Advance to next flux and Jacobian row
			++bndIdx;
			++jac;
		}
	}	
};

typedef MultiComponentLDFFreundlichBindingBase<MultiComponentLDFFreundlichParamHandler> MultiComponentLDFFreundlichBinding;
typedef MultiComponentLDFFreundlichBindingBase<ExtMultiComponentLDFFreundlichParamHandler> ExternalMultiComponentLDFFreundlichBinding;
typedef MultiComponentLDFFreundlichBindingBase<FieldMultiComponentLDFFreundlichParamHandler> FieldMultiComponentLDFFreundlichBinding;

namespace binding
{
	void registerMultiComponentLDFFreundlichModel(std::unordered_map<std::string, std::function<model::IBindingModel*()>>& bindings)
	{
		bindings[MultiComponentLDFFreundlichBinding::identifier()] = []() { return new MultiComponentLDFFreundlichBinding(); };
		bindings[ExternalMultiComponentLDFFreundlichBinding::identifier()] = []() { return new ExternalMultiComponentLDFFreundlichBinding(); };
		bindings[FieldMultiComponentLDFFreundlichBinding::identifier()] = []() { return new FieldMultiComponentLDFFreundlichBinding(); };
	}
}  // namespace binding

}  // namespace model

}  // namespace cadet
