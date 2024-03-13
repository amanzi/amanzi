/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#include <algorithm>

#include "FunctionGridColor.hh"

namespace Amanzi {

int
FunctionGridColor::operator()(const double* x) const
{
  int offset = 0;
  for (int k = dim_ - 1; k >= 0; --k) {
    int i = (int)((x[k] - x0_[k]) / dx_[k]);
    i = std::min(i, count_[k] - 1);
    i = std::max(i, 0);
    offset = i + count_[k] * offset; // offset = 0 on first pass
  }
  return array_[offset];
}

} // namespace Amanzi
