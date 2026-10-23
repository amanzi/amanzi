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

  A faster implementation of not-a-knot 2D tensor-product cubic spline.
*/

#include <array>
#include <vector>
#include <cmath>

#include "lapack.hh"

namespace Amanzi {
namespace WhetStone {

class TensorCubicBSpline2D {
 public:
  TensorCubicBSpline2D(std::vector<double> x, std::vector<double> y) 
    : x_(std::move(x)), y_(std::move(y)) {};

  int build(const std::vector<std::vector<double>>& values);

  std::array<double, 6> evaluate(double x, double y) const;

  // checks
  int checkIncreasing(const std::vector<double>& x) {
    for (int i = 1; i < x.size(); ++i) {
      if (!(x[i] > x[i - 1])) return -1;
    }
    return 0;
  }

 private:
  // Cubic not-a-knot knot vector for interpolation points x[0],...,x[n-1].
  // For cubic splines, not-a-knot removes the first and last interior knots.
  // Knot vector:
  //
  //   x0,x0,x0,x0, x2,x3,...,x[n-3], x[n-1],x[n-1],x[n-1],x[n-1]
  //
  // Number of basis functions is n.
  // Knot vector size is n + p + 1 = n + 4.
  std::vector<double> makeNotAKnotKnots_(const std::vector<double>& x);

  // Build dense collocation matrix A(i,j) = B_j(x_i).
  // Returned in column-major layout for LAPACK.
  std::vector<double> buildCollocationMatrix_(const std::vector<double>& x,
                                              const std::vector<double>& t);

  // Find span for n basis functions of degree p.
  int findSpan_(int nBasis, int p, double u, const std::vector<double>& U) const;

  // DersBasisFuns from The NURBS Book. Returns ders[k][j], where:
  //   k = derivative order, 0 <= k <= nDeriv
  //   j = local basis index, 0 <= j <= p
  //
  // The global basis index is span - p + j.
  std::array<std::array<double, 4>, 3>
  basisDers_(int span, double u, int p, int nDeriv, const std::vector<double>& U) const;

 private:
  int nx_ = 0;
  int ny_ = 0;

  std::vector<double> x_, y_;
  std::vector<double> tx_, ty_;
  std::vector<double> coeff_;  // coeff_[j * nx_ + i]
};

}  // namespace WhetStone
}  // namespace Amanzi
