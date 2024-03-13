/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Utils

*/

#ifndef AMANZI_UTILS_LEAST_SQUARE_HH_
#define AMANZI_UTILS_LEAST_SQUARE_HH_

#include <vector>

namespace Amanzi {
namespace Utils {

double
bestLSfit(const std::vector<double>& h, const std::vector<double>& error);

} // namespace Utils
} // namespace Amanzi

#endif
