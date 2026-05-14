/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  WhetStone, Version 2.2
  Release name: naka-to.

  Not-a-knot 2D tensor-product cubic spline.
*/

#ifndef AMANZI_UTILS_SPLINE_CUBIC_NOT_A_KNOT_2D_HH_
#define AMANZI_UTILS_SPLINE_CUBIC_NOT_A_KNOT_2D_HH_

#include <array>
#include <cmath>
#include <vector>

#include "SplineCubicNotAKnot1D.hh"

namespace Amanzi {
namespace WhetStone {

class SplineCubicNotAKnot2D {
 public:
  SplineCubicNotAKnot2D() {};
  SplineCubicNotAKnot2D(const std::vector<double>& x,
                        const std::vector<double>& y,
                        const std::vector<std::vector<double>>& values) {
    build(x, y, values);
  }
  ~SplineCubicNotAKnot2D() {};

  int build(const std::vector<double>& x,
            const std::vector<double>& y,
            const std::vector<std::vector<double>>& values);

  // returns spline value and five derivatives
  std::array<double, 6> evaluate(double x, double y) const;

 private:
  int findInterval_(const std::vector<double>& x, double v) const;

 private:
  std::vector<double> x_;
  std::vector<double> y_;

  // for each x-interval and polynomial power p (p <= 4), we need
  // to store a spline in x.
  std::vector<SplineCubicNotAKnot1D> y_splines_;
};

} // namespace WhetStone
} // namespace Amanzi

#endif
