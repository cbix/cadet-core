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
#include "LoggingUtils.hpp"
#include "Logging.hpp"

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
  LOG(Debug) << "Field::configure";
  LOG(Debug) << "field dims: " << dimNames;
  size_t size = 1;
  for (std::string dim : dimNames) {
    std::vector dimCoords = paramProvider->getDoubleArray("COORD_" + dim);
    size *= dimCoords.size();
    _shape.push_back(dimCoords.size());
    _dimensions.push_back(dimCoords);
    LOG(Debug) << dim << " size=" << dimCoords.size() << ", coords=" << dimCoords;
  }
  _data = paramProvider->getDoubleArray("DATA");
  if (_data.size() != size) {
    throw std::domain_error("DATA is " + std::to_string(_data.size()) + ", expected " + std::to_string(size));
  }
  // TODO get other parameters

  return true;
}

double Field::valueAtIndex(std::vector<size_t> idx) {
  LOG(Debug) << "valueAtIndex" << idx;
  size_t n = _dimensions.size();
  if (idx.size() != n || n == 0) {
    throw std::domain_error("Number of dimensions doesn't match or is zero");
  }
  size_t index = 0;
  size_t size = 1;
  for (size_t i = 0; i < n; i++) {
    size_t currentIndex = idx[n - i - 1];
    size_t currentDimSize = _dimensions[n - i - 1].size();
    if (currentIndex >= currentDimSize) {
      throw std::out_of_range("Index out of bounds");
    }
    index += currentIndex * size;
    size *= currentDimSize;
  }
  return _data[index];
}

double Field::interpolateValue(std::vector<double> coords) {
  LOG(Debug) << "Field::interpolateValue " << coords;
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

  typedef struct {
    std::vector<size_t> idx;
    double weight;
  } field_bound_sample_t;
  std::vector<field_bound_sample_t> boundSamples;
  for (const field_interp_bounds_t &bound : bounds) {
    if (boundSamples.size() == 0) {
      boundSamples.push_back({{bound.lower}, bound.weight});
      boundSamples.push_back({{bound.upper}, 1 - bound.weight});
      continue;
    }
    size_t n = boundSamples.size();
    for (size_t i = 0; i < n; i++) {
      field_bound_sample_t sampleUpper = boundSamples[i];
      boundSamples[i].idx.push_back(bound.lower);
      boundSamples[i].weight *= bound.weight;
      sampleUpper.idx.push_back(bound.upper);
      sampleUpper.weight *= 1 - bound.weight;
      boundSamples.push_back(sampleUpper);
    }
  }

  double result = 0.0;
  for (field_bound_sample_t s : boundSamples) {
    result += valueAtIndex(s.idx) * s.weight;
  }
  LOG(Debug) << "interpolation result = " << result;
  return result;
}

} // namespace cadet
