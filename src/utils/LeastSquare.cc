/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Konstantin Lipnikov (lipnikov@lanl.gov)
*/

//! <MISSING_ONELINE_DOCSTRING>

#include <cmath>
#include <vector>

#include "LeastSquare.hh"

namespace Amanzi {
namespace Utils {

/* *******************************************************************
 * Calculates best least square fit for data (h[i], error[i]).
 ******************************************************************* */
double
bestLSfit(const std::vector<double>& h, const std::vector<double>& error)
{
  double a = 0.0, b = 0.0, c = 0.0, d = 0.0, tmp1, tmp2;

  int n = h.size();
  for (int i = 0; i < n; i++) {
    tmp1 = log(h[i]);
    tmp2 = log(error[i]);
    a += tmp1;
    b += tmp2;
    c += tmp1 * tmp1;
    d += tmp1 * tmp2;
  }

  return (a * b - n * d) / (a * a - n * c);
}

} // namespace Utils
} // namespace Amanzi
