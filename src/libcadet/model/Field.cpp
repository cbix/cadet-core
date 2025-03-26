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
 * Provides a linearly interpolated external function from moving data points
 * without radial dependence.
 */

#include "cadet/Field.hpp"
#include "cadet/ParameterProvider.hpp"
#include "common/CompilerSpecific.hpp"

#include <algorithm>
#include <functional>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

namespace cadet {

bool Field::configure(IParameterProvider *paramProvider) {
  if (!paramProvider)
    return false;

  std::vector<std::string> dimNames =
      paramProvider->getStringArray("DIMENSIONS");
  for (std::string dim : dimNames) {
    std::vector dimCoords = paramProvider->getDoubleArray("COORD_" + dim);
    _shape.push_back(dimCoords.size());
    _dimensions.push_back(dimCoords);
  }
  std::vector<double> _data = paramProvider->getDoubleArray("DATA");
  // TODO get other parameters

  return true;
}

double Field::interpolateValue(std::vector<double> coords) {
  size_t n = _shape.size();
  if (coords.size() != n) {
    throw std::domain_error("Number of dimensions doesn't match");
  }
  for (size_t i = 0; i < n; i++) {
    double coord = coords.at(i);
    std::vector<double> dim = _dimensions.at(i);
    if (coord < dim.front() || coord > dim.back()) {
      throw std::range_error("Coordinate out of bounds");
    }
  }

  typedef struct {
    size_t lower, upper;
    double weight;
  } field_interp_bounds_t;
  std::vector<field_interp_bounds_t> bounds;
  for (size_t i = 0; i < n; i++) {
    double coord = coords.at(i);
    std::vector<double> dim = _dimensions.at(i);
    // find bounds for coordinate
    auto lower = std::upper_bound(dim.begin(), dim.end(), coord) - 1;
    size_t l = lower - dim.begin();
    size_t r = l + 1 >= dim.size() ? l : l + 1;

    double coord_l = dim.at(l);
    double coord_r = dim.at(r);
    double weight = 0.0;
    if (coord_r - coord_l > 0) {
      weight = (coord - coord_l) / (coord_r - coord_l);
    }
    bounds.push_back({l, r, weight});
  }

  // TODO sample bounds

  double result = 0.0;

  // TODO interpolate

  return result;
}

} // namespace cadet
