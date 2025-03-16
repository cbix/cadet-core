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
#include <string>
#include <unordered_map>
#include <vector>

namespace cadet {

bool Field::configure(IParameterProvider *paramProvider) {
  if (!paramProvider)
    return false;

  std::vector<std::string> dimNames =
      paramProvider->getStringArray("DIMENSIONS");
  std::vector<size_t> shape;
  for (std::string dim : dimNames) {
    std::vector dimCoords = paramProvider->getDoubleArray("COORD_" + dim);
    shape.push_back(dimCoords.size());
  }
  std::vector<double> _data = paramProvider->getDoubleArray("DATA");

  return true;
}

} // namespace cadet
