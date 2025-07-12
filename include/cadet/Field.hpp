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
 * Fields are external functions of arbitrary dimension and can be used for any
 * parameter
 */

#ifndef LIBCADET_FIELD_HPP_
#define LIBCADET_FIELD_HPP_

#include "cadet/LibExportImport.hpp"
#include "cadet/cadetCompilerInfo.hpp"

#include <map>
#include <string>
#include <vector>

namespace cadet {

class IParameterProvider;

// typedef std::pair<std::string, std::vector<field_val_t>> field_dim_t;

/**
 * @brief Field class
 */
class CADET_API Field {
public:
  Field() {}
  ~Field() CADET_NOEXCEPT {}

  /**
   * @brief Configures the field by extracting all parameters from the given @p
   * paramProvider
   * @details The scope of the cadet::IParameterProvider is left unchanged on
   * return.
   *
   * @param [in] paramProvider Pointer to parameter provider (may be @c nullptr)
   * @return @c true if the configuration was successful, otherwise @c false
   */
  bool configure(IParameterProvider *paramProvider);

  std::vector<int> dimensionMap(std::vector<std::string> dims);

  double valueAtIndex(std::vector<size_t> idx);

  double interpolateValue(std::vector<double> coords);

  static const char *identifier() { return "LINEAR_INTERP_FIELD"; }

private:
  std::vector<double> _time;
  std::vector<size_t> _shape;
  //std::vector<std::pair<std::string, std::vector<double>>> _dimensions;
  std::vector<std::vector<double>> _dimensions;
  std::vector<std::string> _dimNames;
  std::map<std::string, size_t> _dimNamesMap;
  std::vector<double> _data;
};

} // namespace cadet

#endif // LIBCADET_FIELD_HPP_
