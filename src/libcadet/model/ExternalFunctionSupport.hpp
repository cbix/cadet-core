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

/**
 * @file 
 * Provides support functions for defining externally dependent models.
 */

#ifndef LIBCADET_EXTFUNSUPPORT_HPP_
#define LIBCADET_EXTFUNSUPPORT_HPP_

#include "cadet/ExternalFunction.hpp"
#include "cadet/Exceptions.hpp"
#include "cadet/Field.hpp"
#include "Logging.hpp"
#include "LoggingUtils.hpp"
#include "SimulationTypes.hpp"

#include <vector>
#include <algorithm>
#include <string>

namespace cadet
{

template <typename ValType>
inline bool readScalarParameterOrArray(std::vector<ValType>& dest, IParameterProvider& paramProvider, const std::string& dataSet, unsigned int nExpand);

namespace model
{

	/**
	 * @brief Base class of model parameter storage classes that do not depend on external functions
	 */
	struct ConstParamHandlerBase
	{
		inline void setFields(Field** fields, unsigned int size) { }

		/**
		 * @brief Sets external functions for this model
		 * @param [in] extFuns Pointer to array of IExternalFunction objects of size @p size
		 * @param [in] size Number of elements in the IExternalFunction array @p extFuns
		 */
		inline void setExternalFunctions(IExternalFunction** extFuns, unsigned int size) { }

		/**
		 * @brief Returns whether the model parameters depend on time
		 * @details Model parameters that do not use external functions do not depend on time.
		 * @return @c true if the model parameters depends on time, otherwise @c false
		 */
		static bool dependsOnTime() CADET_NOEXCEPT { return false; }

		/**
		 * @brief Returns whether workspace is required for the parameters
		 * @details Model parameters that do not use external functions do not require a workspace for parameters.
		 * @return @c true if the model parameters require a workspace, otherwise @c false
		 */
		static bool requiresWorkspace() CADET_NOEXCEPT { return false; }

		/**
		 * @brief Returns how much memory is required for caching in bytes
		 * @details Memory size in bytes.
		 * @param [in] nComp Number of components
		 * @param [in] totalNumBoundStates Total number of bound states
		 * @param [in] nBoundStates Array with bound states for each component
		 * @return Memory size in bytes
		 */
		inline std::size_t cacheSize(unsigned int nComp, unsigned int totalNumBoundStates, unsigned int const* nBoundStates) const CADET_NOEXCEPT { return 0; }

		/**
		 * @brief Returns how much memory is required for caching in bytes
		 * @details Memory size in bytes.
		 * @param [in] nReactions Number of reactions
		 * @param [in] nComp Number of components
		 * @param [in] totalNumBoundStates Total number of bound states
		 * @return Memory size in bytes
		 */
		inline std::size_t cacheSize(unsigned int nReactions, unsigned int nComp, unsigned int totalNumBoundStates) const CADET_NOEXCEPT { return 0; }
	};

	/**
	 * @brief Base class for externally dependent model parameter classes
	 * @details Configures and stores the external function used for the model parameters.
	 */
	struct ExternalParamHandlerBase
	{	
	public:

		inline void setFields(Field** field, int size)
		{
			// only implemented for FieldParamHandlerBase
		}

		/**
		 * @brief Sets external functions for this model
		 * @param [in] extFuns Pointer to array of IExternalFunction objects of size @p size
		 * @param [in] size Number of elements in the IExternalFunction array @p extFuns
		 */
		inline void setExternalFunctions(IExternalFunction** extFuns, int size)
		{
			_extFun.clear();
			_extFun.resize(_extFunIndex.size(), nullptr);
			for (std::size_t i = 0; i < _extFunIndex.size(); ++i)
			{
				if ((_extFunIndex[i] >= 0) && (_extFunIndex[i] < size))
					_extFun[i] = extFuns[_extFunIndex[i]];
				else
				{
					_extFun[i] = nullptr;
					LOG(Warning) << "Index " << _extFunIndex[i] << " exceeds number of passed external functions (" << size << "), external dependence is ignored";
				}
			}
		}

		/**
		 * @brief Returns whether the model parameters depend on time
		 * @details Model parameters that do not use external functions do not depend on time.
		 * @return @c true if the model parameters depends on time, otherwise @c false
		 */
		static bool dependsOnTime() CADET_NOEXCEPT { return true; }

		/**
		 * @brief Returns whether workspace is required for the parameters
		 * @details Model parameters that do not use external functions do not require a workspace for parameters.
		 * @return @c true if the model parameters require a workspace, otherwise @c false
		 */
		static bool requiresWorkspace() CADET_NOEXCEPT { return true; }

	protected:

		std::vector<IExternalFunction*> _extFun; //!< Pointer to the external function
		std::vector<int> _extFunIndex; //!< Index to the external function

		ExternalParamHandlerBase() : _extFun(), _extFunIndex() { }
		
		/**
		 * @brief Configures the external data source of this externally dependent parameter set
		 * @param [in] paramProvider Parameter provider
		 * @param [in] nParams Number of externally dependent parameters (also size of buffer)
		 */
		inline void configure(IParameterProvider& paramProvider, unsigned int nParams)
		{			
			std::vector<int> idx;
			if (paramProvider.exists("EXTFUN"))
				idx = paramProvider.getIntArray("EXTFUN");

			if (idx.size() >= nParams)
				_extFunIndex = idx;
			else
			{
				_extFunIndex.resize(nParams);
				if (!idx.empty())
				{
					// Use one external function for all parameters
					std::fill(_extFunIndex.begin(), _extFunIndex.end(), idx[0]);
				}
				else
				{
					// There is no external dependence configured
					std::fill(_extFunIndex.begin(), _extFunIndex.end(), -1);
				}
			}
		}

		/**
		 * @brief Evaluates the external functions for the different parameters
		 * @param [in] t Current time
		 * @param [in] z Axial coordinate in the column
		 * @param [in] r Radial coordinate in the bead
		 * @param [in] secIdx Index of the current section
		 * @param [in] nParams Number of externally dependent parameters (also size of buffer)
		 * @param [out] buffer Buffer that holds function evaluations
		 */
		inline void evaluateExternalFunctions(double t, unsigned int secIdx, const ColumnPosition& colPos, unsigned int nParams, double* buffer) const
		{
			for (unsigned int i = 0; i < nParams; ++i)
			{
				IExternalFunction* const fun = _extFun[i];
				if (fun)
					buffer[i] = fun->externalProfile(t, colPos.axial, colPos.radial, colPos.particle, secIdx);
				else
					buffer[i] = 0.0;
			}
		}

		/**
		 * @brief Evaluates the time derivative of the external functions for the different parameters
		 * @param [in] t Current time
		 * @param [in] z Axial coordinate in the column
		 * @param [in] r Radial coordinate in the bead
		 * @param [in] secIdx Index of the current section
		 * @param [in] nParams Number of externally dependent parameters (also size of buffer)
		 * @param [out] buffer Buffer that holds time derivatives of each external function
		 */
		inline void evaluateTimeDerivativeExternalFunctions(double t, unsigned int secIdx, const ColumnPosition& colPos, unsigned int nParams, double* buffer) const
		{
			for (unsigned int i = 0; i < nParams; ++i)
			{
				IExternalFunction* const fun = _extFun[i];
				if (fun)
					buffer[i] = fun->timeDerivative(t, colPos.axial, colPos.radial, colPos.particle, secIdx);
				else
					buffer[i] = 0.0;
			}
		}
	};

	/**
	 * @brief Base class of model parameter storage classes that depend on a field
	 */
	struct FieldParamHandlerBase
	{
	public:

		inline void setFields(Field** fields, unsigned int size) { }

		/*
		inline void setExternalFunctions(IExternalFunction** extFuns, int size)
		{
			_extFun.clear();
			_extFun.resize(_extFunIndex.size(), nullptr);
			for (std::size_t i = 0; i < _extFunIndex.size(); ++i)
			{
				if ((_extFunIndex[i] >= 0) && (_extFunIndex[i] < size))
					_extFun[i] = extFuns[_extFunIndex[i]];
				else
				{
					_extFun[i] = nullptr;
					LOG(Warning) << "Index " << _extFunIndex[i] << " exceeds number of passed external functions (" << size << "), external dependence is ignored";
				}
			}
		}
		*/

		inline void setFields(Field** fields, int size)
		{
			_fields.clear();
			_fields.resize(_fieldIndexes.size(), nullptr);
			for (size_t i = 0; i < _fieldIndexes.size(); ++i)
			{
				if (_fieldIndexes[i] >= 0 && _fieldIndexes[i] < size)
					_fields[i] = fields[_fieldIndexes[i]];
				else
				{
					_fields[i] = nullptr;
					LOG(Warning) << "Index " << _fieldIndexes[i] << " exceeds number of passed fields (" << size << "), field dependence is ignored";
				}
			}
		}

		static bool dependsOnTime() CADET_NOEXCEPT { return true; }
		static bool requiresWorkspace() CADET_NOEXCEPT { return true; }

	protected:

		std::vector<std::vector<Field*>> _fields; //!< Pointer to the field, by parameter/component
		std::vector<std::vector<int>> _fieldIndexes; //!< Index to the field, by parameter/component
		std::vector<std::vector<int>> _dimensionMaps; //!< Mapping model to field dimensions, by parameter

		FieldParamHandlerBase() : _fields(), _fieldIndexes() { }

		/**
		 * @brief Configures the external data source of this externally dependent parameter set
		 * @param [in] paramProvider Parameter provider
		 * @param [in] nParams Number of externally dependent parameters (also size of buffer)
		 * @param [in] dimensions 
		 */
		inline void configure(IParameterProvider& paramProvider, std::vector<std::string> params, std::vector<std::string> dimensions)
		{			
			for (size_t i = 0; i < params.size(); ++i)
			{
				std::vector<int> idx;
				std::string param = params[i] + "_EXTFUN";
				if (!paramProvider.exists(param))
					idx = paramProvider.getIntArray(param);

				if (idx.size() >= params.size())
					_fieldIndexes[i] = idx;
				else
				{
					_fieldIndexes.resize(params.size());
					if (!idx.empty())
					{
						// Use one external function for all components
						std::fill(_fieldIndexes[i].begin(), _fieldIndexes[i].end(), idx[0]);
					}
					else
					{
						// There is no external dependence configured
						std::fill(_fieldIndexes[i].begin(), _fieldIndexes[i].end(), -1);
					}
				}
			}


			for (unsigned int i = 0; i < _fields.size()) {
				Field* field = _fields[i];
				if (field) {
					_dimensionMaps[i] = field->dimensionMap(dimensions);
				}
			}
		}

		inline void evaluateField(std::vector<double> coords, unsigned int nParams, double* buffer) const
		{

			for (unsigned int i = 0; i < nParams; ++i)
			{
				Field* const field = _fields[i];
				std::vector<int> dimMap = _dimensionMaps[i];
				if (field)
				{
					auto mappedCoords = std::vector<double>();
					for (unsigned int dimIdx = 0; dimIdx < coords.size(); ++dimIdx)
					{
						mappedCoords[dimIdx] = dimMap[dimIdx];
					}
					buffer[i] = field->interpolateValue(mappedCoords);
				}
				else
				{
					buffer[i] = 0.0;
				}
			}
		}

		/**
		 * @brief Evaluates the time derivative of the external functions for the different parameters
		 * @param [in] t Current time
		 * @param [in] z Axial coordinate in the column
		 * @param [in] r Radial coordinate in the bead
		 * @param [in] secIdx Index of the current section
		 * @param [in] nParams Number of externally dependent parameters (also size of buffer)
		 * @param [out] buffer Buffer that holds time derivatives of each external function
		 */
		inline void evaluateTimeDerivativeExternalFunctions(double t, unsigned int secIdx, const ColumnPosition& colPos, unsigned int nParams, double* buffer) const
		{
			for (unsigned int i = 0; i < nParams; ++i)
			{
				IExternalFunction* const fun = _extFun[i];
				if (fun)
					buffer[i] = fun->timeDerivative(t, colPos.axial, colPos.radial, colPos.particle, secIdx);
				else
					buffer[i] = 0.0;
			}
		}

	};

}  // namespace model

}  // namespace cadet

#endif  // LIBCADET_EXTFUNSUPPORT_HPP_
