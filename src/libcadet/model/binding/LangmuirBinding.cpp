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
	"name": "LangmuirParamHandler",
	"externalName": "ExtLangmuirParamHandler",
	"parameters":
		[
			{ "type": "ScalarComponentDependentParameter", "varName": "kA", "confName": "MCL_KA"},
			{ "type": "ScalarComponentDependentParameter", "varName": "kD", "confName": "MCL_KD"},
			{ "type": "ScalarComponentDependentParameter", "varName": "qMax", "confName": "MCL_QMAX"}
		]
}
</codegen>*/

/* Parameter description
 ------------------------
 kA = Adsorption rate
 kD = Desorption rate
 qMax = Capacity
*/

namespace cadet
{

namespace model
{

class FieldLangmuirParamHandler : public FieldParamHandlerBase
{
public:
	struct VariableParams
	{
		util::LocalVector<active> kA;
		util::LocalVector<active> kD;
		util::LocalVector<active> qMax;
	};

	typedef VariableParams params_t;
	typedef ConstBufferedScalar<params_t> ParamsHandle;

	static inline const char* identifier() CADET_NOEXCEPT;

	FieldLangmuirParamHandler() CADET_NOEXCEPT { }

	inline bool configure(IParameterProvider& paramProvider, unsigned int nComp, unsigned int const* nBoundStates)
	{
		_kA.configure("MCL_KA", paramProvider, nComp, nBoundStates);
		_kD.configure("MCL_KD", paramProvider, nComp, nBoundStates);
		_qMax.configure("MCL_QMAX", paramProvider, nComp, nBoundStates);

		FieldParamHandlerBase::configure(paramProvider, { "TIME", "AXIAL", "RADIAL", "PARTICLE" });

		return validateConfig(nComp, nBoundStates);
	}

	inline void registerParameters(std::unordered_map<ParameterId, active*>& parameters, UnitOpIdx unitOpIdx, ParticleTypeIdx parTypeIdx, unsigned int nComp, unsigned int const* nBoundStates)
	{
		_kA.registerParam("MCL_KA", parameters, unitOpIdx, parTypeIdx, nComp, nBoundStates);
		_kD.registerParam("MCL_KD", parameters, unitOpIdx, parTypeIdx, nComp, nBoundStates);
		_qMax.registerParam("MCL_QMAX", parameters, unitOpIdx, parTypeIdx, nComp, nBoundStates);
	}

	inline void reserve(unsigned int numElem, unsigned int numSlices, unsigned int nComp, unsigned int const* nBoundStates) \
	{
		_kA.reserve(numElem, numSlices, nComp, nBoundStates);
		_kD.reserve(numElem, numSlices, nComp, nBoundStates);
		_qMax.reserve(numElem, numSlices, nComp, nBoundStates);
	}

	inline ParamsHandle update(double t, unsigned int secIdx, const ColumnPosition& colPos, unsigned int nComp, unsigned int const* nBoundStates, LinearBufferAllocator& workSpace) const
	{
		LOG(Debug) << "LangmuirParamHandler update(" << t << ", " << secIdx << ", "
			<< colPos.axial << "/" << colPos.radial << "/" << colPos.particle << ", "
			<< nComp << ", ...)";
		BufferedScalar<params_t> localParams = workSpace.scalar<params_t>();
		BufferedArray<double> fieldBuffer = workSpace.array<double>(_fields.size());

		evaluateField({t, colPos.axial, colPos.radial, colPos.particle}, static_cast<double*>(fieldBuffer));

		_kA.prepareCache(localParams->kA, workSpace);
		_kA.update(cadet::util::dataOfLocalVersion(localParams->kA), &fieldBuffer[0], nComp, nBoundStates);
		_kD.prepareCache(localParams->kD, workSpace);
		_kD.update(cadet::util::dataOfLocalVersion(localParams->kD), &fieldBuffer[0], nComp, nBoundStates);
		_qMax.prepareCache(localParams->qMax, workSpace);
		_qMax.update(cadet::util::dataOfLocalVersion(localParams->qMax), &fieldBuffer[0], nComp, nBoundStates);
		LOG(Debug) << "  kA,kD,qMax = " << static_cast<double>(localParams->kA[0]) << ", " << static_cast<double>(localParams->kD[0]) << ", " << static_cast<double>(localParams->qMax[0]);

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

		_kA.prepareCache(localParams->kA, workSpace);
		_kA.update(cadet::util::dataOfLocalVersion(localParams->kA), &fieldBuffer[0], nComp, nBoundStates);

		_kA.prepareCache(p->kA, workSpace);
		_kA.updateTimeDerivative(cadet::util::dataOfLocalVersion(p->kA), &fieldBuffer[0], &fieldDerivBuffer[0], nComp, nBoundStates);

		_kD.prepareCache(localParams->kD, workSpace);
		_kD.update(cadet::util::dataOfLocalVersion(localParams->kD), &fieldBuffer[0], nComp, nBoundStates);

		_kD.prepareCache(p->kD, workSpace);
		_kD.updateTimeDerivative(cadet::util::dataOfLocalVersion(p->kD), &fieldBuffer[0], &fieldDerivBuffer[0], nComp, nBoundStates);

		_qMax.prepareCache(localParams->qMax, workSpace);
		_qMax.update(cadet::util::dataOfLocalVersion(localParams->qMax), &fieldBuffer[0], nComp, nBoundStates);

		_qMax.prepareCache(p->qMax, workSpace);
		_qMax.updateTimeDerivative(cadet::util::dataOfLocalVersion(p->qMax), &fieldBuffer[0], &fieldDerivBuffer[0], nComp, nBoundStates);

		return std::make_tuple<ParamsHandle, ParamsHandle>(std::move(localParams), std::move(p));
	}

	inline std::size_t cacheSize(unsigned int nComp, unsigned int totalNumBoundStates, unsigned int const* nBoundStates) const CADET_NOEXCEPT
	{
		return 2 * sizeof(params_t) + alignof(params_t) + 2 * 3 * sizeof(double) + alignof(double) + 2 * (
			_kA.additionalDynamicMemory(nComp, totalNumBoundStates, nBoundStates)
			+ _kD.additionalDynamicMemory(nComp, totalNumBoundStates, nBoundStates)
			+ _qMax.additionalDynamicMemory(nComp, totalNumBoundStates, nBoundStates)
		);
	}

protected:
	inline bool validateConfig(unsigned int nComp, unsigned int const* nBoundStates);
	FieldScalarComponentDependentParameter _kA;
	FieldScalarComponentDependentParameter _kD;
	FieldScalarComponentDependentParameter _qMax;
};

} // namespace model
} // namespace cadet

namespace cadet
{

namespace model
{

inline const char* LangmuirParamHandler::identifier() CADET_NOEXCEPT { return "MULTI_COMPONENT_LANGMUIR"; }

inline bool LangmuirParamHandler::validateConfig(unsigned int nComp, unsigned int const* nBoundStates)
{
	if ((_kA.size() != _kD.size()) || (_kA.size() != _qMax.size()) || (_kA.size() < nComp))
		throw InvalidParameterException("MCL_KA, MCL_KD, and MCL_QMAX have to have the same size");

	return true;
}

inline const char* ExtLangmuirParamHandler::identifier() CADET_NOEXCEPT { return "EXT_MULTI_COMPONENT_LANGMUIR"; }

inline bool ExtLangmuirParamHandler::validateConfig(unsigned int nComp, unsigned int const* nBoundStates)
{
	if ((_kA.size() != _kD.size()) || (_kA.size() != _qMax.size()) || (_kA.size() < nComp))
		throw InvalidParameterException("EXT_MCL_KA, EXT_MCL_KD, and EXT_MCL_QMAX have to have the same size");

	return true;
}

inline const char* FieldLangmuirParamHandler::identifier() CADET_NOEXCEPT { return "FIELD_MULTI_COMPONENT_LANGMUIR"; }

inline bool FieldLangmuirParamHandler::validateConfig(unsigned int nComp, unsigned int const* nBoundStates)
{
	if ((_kA.size() != _kD.size()) || (_kA.size() != _qMax.size()) || (_kA.size() < nComp))
		throw InvalidParameterException("FIELD_MCL_KA, FIELD_MCL_KD, and FIELD_MCL_QMAX have to have the same size");

	return true;
}


/**
 * @brief Defines the multi component Langmuir binding model
 * @details Implements the Langmuir adsorption model: \f[ \begin{align} 
 *              \frac{\mathrm{d}q_i}{\mathrm{d}t} &= k_{a,i} c_{p,i} q_{\text{max},i} \left( 1 - \sum_j \frac{q_j}{q_{\text{max},j}} \right) - k_{d,i} q_i
 *          \end{align} \f]
 *          Multiple bound states are not supported. 
 *          Components without bound state (i.e., non-binding components) are supported.
 *          
 *          See @cite Langmuir1916.
 * @tparam ParamHandler_t Type that can add support for external function dependence
 */
template <class ParamHandler_t>
class LangmuirBindingBase : public ParamHandlerBindingModelBase<ParamHandler_t>
{
public:

	LangmuirBindingBase() { }
	virtual ~LangmuirBindingBase() CADET_NOEXCEPT { }

	static const char* identifier() { return ParamHandler_t::identifier(); }

	virtual void timeDerivativeQuasiStationaryFluxes(double t, unsigned int secIdx, const ColumnPosition& colPos, double const* yCp, double const* y, double* dResDt, LinearBufferAllocator workSpace) const
	{
		if (!this->hasQuasiStationaryReactions())
			return;

		if (!ParamHandler_t::dependsOnTime())
			return;

		// Update external function and compute time derivative of parameters
		typename ParamHandler_t::ParamsHandle p;
		typename ParamHandler_t::ParamsHandle dpDt;
		std::tie(p, dpDt) = _paramHandler.updateTimeDerivative(t, secIdx, colPos, _nComp, _nBoundStates, workSpace);

		// Protein fluxes: -k_{a,i} * c_{p,i} * q_{max,i} * (1 - \sum_j q_j / q_{max,j}) + k_{d,i} * q_i
		double qSum = 1.0;
		double qSumT = 0.0;
		unsigned int bndIdx = 0;
		for (int i = 0; i < _nComp; ++i)
		{
			// Skip components without bound states (bound state index bndIdx is not advanced)
			if (_nBoundStates[i] == 0)
				continue;

			const double summand = y[bndIdx] / static_cast<double>(p->qMax[i]);
			qSum -= summand;
			qSumT += summand / static_cast<double>(p->qMax[i]) * static_cast<double>(dpDt->qMax[i]);

			// Next bound component
			++bndIdx;
		}

		bndIdx = 0;
		for (int i = 0; i < _nComp; ++i)
		{
			// Skip components without bound states (bound state index bndIdx is not advanced)
			if (_nBoundStates[i] == 0)
				continue;

			// Residual
			dResDt[bndIdx] = static_cast<double>(dpDt->kD[i]) * y[bndIdx] 
				- yCp[i] * (static_cast<double>(dpDt->kA[i]) * static_cast<double>(p->qMax[i]) * qSum
				           + static_cast<double>(p->kA[i]) * static_cast<double>(dpDt->qMax[i]) * qSum
				           + static_cast<double>(p->kA[i]) * static_cast<double>(p->qMax[i]) * qSumT);

			// Next bound component
			++bndIdx;
		}
	}

	CADET_BINDINGMODELBASE_BOILERPLATE

protected:
	using ParamHandlerBindingModelBase<ParamHandler_t>::_paramHandler;
	using ParamHandlerBindingModelBase<ParamHandler_t>::_reactionQuasistationarity;
	using ParamHandlerBindingModelBase<ParamHandler_t>::_nComp;
	using ParamHandlerBindingModelBase<ParamHandler_t>::_nBoundStates;

	virtual bool implementsAnalyticJacobian() const CADET_NOEXCEPT { return true; }

	template <typename StateType, typename CpStateType, typename ResidualType, typename ParamType>
	int fluxImpl(double t, unsigned int secIdx, const ColumnPosition& colPos, StateType const* y,
		CpStateType const* yCp, ResidualType* res, LinearBufferAllocator workSpace) const
	{
		typename ParamHandler_t::ParamsHandle const p = _paramHandler.update(t, secIdx, colPos, _nComp, _nBoundStates, workSpace);

		// Protein fluxes: -k_{a,i} * c_{p,i} * q_{max,i} * (1 - \sum_j q_j / q_{max,j}) + k_{d,i} * q_i
		ResidualType qSum = 1.0;
		unsigned int bndIdx = 0;
		for (int i = 0; i < _nComp; ++i)
		{
			// Skip components without bound states (bound state index bndIdx is not advanced)
			if (_nBoundStates[i] == 0)
				continue;

			qSum -= y[bndIdx] / static_cast<ParamType>(p->qMax[i]);

			// Next bound component
			++bndIdx;
		}

		bndIdx = 0;
		for (int i = 0; i < _nComp; ++i)
		{
			// Skip components without bound states (bound state index bndIdx is not advanced)
			if (_nBoundStates[i] == 0)
				continue;

			// Residual
			res[bndIdx] = static_cast<ParamType>(p->kD[i]) * y[bndIdx] - static_cast<ParamType>(p->kA[i]) * yCp[i] * static_cast<ParamType>(p->qMax[i]) * qSum;

			// Next bound component
			++bndIdx;
		}

		return 0;
	}

	template <typename RowIterator>
	void jacobianImpl(double t, unsigned int secIdx, const ColumnPosition& colPos, double const* y, double const* yCp, int offsetCp, RowIterator jac, LinearBufferAllocator workSpace) const
	{
		typename ParamHandler_t::ParamsHandle const p = _paramHandler.update(t, secIdx, colPos, _nComp, _nBoundStates, workSpace);

		// Protein fluxes: -k_{a,i} * c_{p,i} * q_{max,i} * (1 - \sum_j q_j / q_{max,j}) + k_{d,i} * q_i
		double qSum = 1.0;
		int bndIdx = 0;
		for (int i = 0; i < _nComp; ++i)
		{
			// Skip components without bound states (bound state index bndIdx is not advanced)
			if (_nBoundStates[i] == 0)
				continue;

			qSum -= y[bndIdx] / static_cast<double>(p->qMax[i]);

			// Next bound component
			++bndIdx;
		}

		bndIdx = 0;
		for (int i = 0; i < _nComp; ++i)
		{
			// Skip components without bound states (bound state index bndIdx is not advanced)
			if (_nBoundStates[i] == 0)
				continue;

			const double ka = static_cast<double>(p->kA[i]);
			const double kd = static_cast<double>(p->kD[i]);

			// dres_i / dc_{p,i}
			jac[i - bndIdx - offsetCp] = -ka * static_cast<double>(p->qMax[i]) * qSum;
			// Getting to c_{p,i}: -bndIdx takes us to q_0, another -offsetCp to c_{p,0} and a +i to c_{p,i}.
			//                     This means jac[i - bndIdx - offsetCp] corresponds to c_{p,i}.

			// Fill dres_i / dq_j
			int bndIdx2 = 0;
			for (int j = 0; j < _nComp; ++j)
			{
				// Skip components without bound states (bound state index bndIdx is not advanced)
				if (_nBoundStates[j] == 0)
					continue;

				// dres_i / dq_j
				jac[bndIdx2 - bndIdx] = ka * yCp[i] * static_cast<double>(p->qMax[i]) / static_cast<double>(p->qMax[j]);
				// Getting to q_j: -bndIdx takes us to q_0, another +bndIdx2 to q_j. This means jac[bndIdx2 - bndIdx] corresponds to q_j.

				++bndIdx2;
			}

			// Add to dres_i / dq_i
			jac[0] += kd;

			// Advance to next flux and Jacobian row
			++bndIdx;
			++jac;
		}
	}	
};

typedef LangmuirBindingBase<LangmuirParamHandler> LangmuirBinding;
typedef LangmuirBindingBase<ExtLangmuirParamHandler> ExternalLangmuirBinding;
typedef LangmuirBindingBase<FieldLangmuirParamHandler> FieldLangmuirBinding;

namespace binding
{
	void registerLangmuirModel(std::unordered_map<std::string, std::function<model::IBindingModel*()>>& bindings)
	{
		bindings[LangmuirBinding::identifier()] = []() { return new LangmuirBinding(); };
		bindings[ExternalLangmuirBinding::identifier()] = []() { return new ExternalLangmuirBinding(); };
		bindings[FieldLangmuirBinding::identifier()] = []() { return new FieldLangmuirBinding(); };
	}
}  // namespace binding

}  // namespace model

}  // namespace cadet
