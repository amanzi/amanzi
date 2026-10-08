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

  Not-a-knot 1D cubic spline.
*/

#ifndef AMANZI_UTILS_SPLINE_CUBIC_NOT_A_KNOT_1D_HH_
#define AMANZI_UTILS_SPLINE_CUBIC_NOT_A_KNOT_1D_HH_

#include <array>
#include <cmath>
#include <vector>

namespace Amanzi {
namespace WhetStone {

struct SplineCoef {
  double a; // constant
  double b; // linear
  double c; // quadratic
  double d; // cubic
};

class SplineCubicNotAKnot1D {
 public:
  SplineCubicNotAKnot1D() = default;
  SplineCubicNotAKnot1D(const std::vector<double>& x, const std::vector<double>& values) {
    build(x, values);
  }
  ~SplineCubicNotAKnot1D() {};

  int build(const std::vector<double>& x, const std::vector<double>& values);

  // returns spline value, 1st and 2nd derivatives
  std::array<double, 3> evaluate(double x) const;

  // access
  const SplineCoef& coeff(int interval) const { return coeff_.at(interval); }

 private:
  int findInterval_(double x) const;

  bool checkMonotonicity_(const std::vector<double>& x) const {
    for (int i = 1; i < x.size(); ++i) {
      if (x[i] <= x[i - 1]) return false;
    }
    return true;
  }

 private:
  std::vector<double> x_;
  std::vector<SplineCoef> coeff_;
};

} // namespace WhetStone
} // namespace Amanzi

#endif
