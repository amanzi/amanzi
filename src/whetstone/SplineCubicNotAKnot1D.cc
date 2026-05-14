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

#include <iostream>

#include "lapack.hh"
#include "SplineCubicNotAKnot1D.hh"

namespace Amanzi {
namespace WhetStone {

/* ******************************************************************
* Build not-a-knot 1D cubic spline
****************************************************************** */
int
SplineCubicNotAKnot1D::build(const std::vector<double>& x, const std::vector<double>& values)
{
  if (x.size() != values.size()) return -1;

  int n(x.size());
  if (n < 4) return -2; // not-a-knot cubic spline requires at least 4 nodes

  x_ = x;
  bool ok = checkMonotonicity_(x_);
  if (!ok) return -3; // invalid 1D mesh

  std::vector<double> h(n - 1);
  for (int i = 0; i < n - 1; ++i) {
    h[i] = x_[i + 1] - x_[i];
  }

  // boundary unknowns are eliminated using not-a-knot conditions.
  int r(n - 2); // number of reduced unknowns
  std::vector<double> dl(r - 1), d(r), du(r - 1), rhs(r);

  // the normal interior equations for i = 1, ..., n-2
  for (int i = 1; i <= n - 2; ++i) {
    int k = i - 1;
    d[k] = 2.0 * (h[i - 1] + h[i]);

    if (k > 0) dl[k - 1] = h[i - 1];
    if (k < r - 1) du[k] = h[i];

    double left = (values[i] - values[i - 1]) / h[i - 1];
    double right = (values[i + 1] - values[i]) / h[i];

    rhs[k] = 6.0 * (right - left);
  }

  // left not-a-knot equation: (M1 - M0) / h0 = (M2 - M1) / h1
  // express M0 and substitute to the first interior equation
  d[0] += h[0] * (h[0] + h[1]) / h[1];
  if (r > 1) du[0] += -h[0] * h[0] / h[1];

  // right not-a-knot equation: (M[n-1] - M[n-2]) / h[n-2]
  // express xM[n-1] and substitute into the last interior equation
  d[r - 1] += h[n - 2] * (h[n - 3] + h[n - 2]) / h[n - 3];
  if (r > 1) dl[r - 2] += -h[n - 2] * h[n - 2] / h[n - 3];

  int nrhs(1), ldb(r), info(0);
  DGTSV_F77(&r, &nrhs, dl.data(), d.data(), du.data(), rhs.data(), &ldb, &info);

  if (info < 0) return -4;  // LAPACK illegal argument
  if (info > 0) return -5;  // LAPACK singular tridiagonal system

  // reconstruct the vector of second derivatives M[0], ..., M[n-1].
  std::vector<double> M(n);

  for (int k = 0; k < r; ++k) M[k + 1] = rhs[k];
  M[0] = ((h[0] + h[1]) * M[1] - h[0] * M[2]) / h[1];
  M[n - 1] = ((h[n - 3] + h[n - 2]) * M[n - 2] - h[n - 2] * M[n - 3]) / h[n - 3];

  // convert second derivatives into local cubic coefficients on interval [x_i, x_{i+1}]
  // S_i(u) = a_i + b_i u + c_i u^2 + d_i u^3, where u = x - x_i.
  coeff_.resize(n - 1);

  for (int i = 0; i < n - 1; ++i) {
    double hi = h[i];
    coeff_[i].a = values[i];
    coeff_[i].b = (values[i + 1] - values[i]) / hi - hi * (2.0 * M[i] + M[i + 1]) / 6.0;
    coeff_[i].c = M[i] / 2.0;
    coeff_[i].d = (M[i + 1] - M[i]) / (6.0 * hi);
  }

  return 0;
}


/* ******************************************************************
* Return value, 1st and 2nd derivatives
****************************************************************** */
std::array<double, 3>
SplineCubicNotAKnot1D::evaluate(double x) const
{
  int i = findInterval_(x);
  double u = x - x_[i];
  auto& q = coeff_[i];

  double value = q.a + u * (q.b + u * (q.c + u * q.d));
  double first = q.b + u * (2.0 * q.c + u * 3.0 * q.d);
  double second = 2.0 * q.c + 6.0 * q.d * u;

  return {value, first, second};
}


/* ******************************************************************
* 
****************************************************************** */
int 
SplineCubicNotAKnot1D::findInterval_(double x) const
{
  if (x <= x_.front()) return 0;
  if (x >= x_.back()) return static_cast<int>(x_.size()) - 2;

  auto it = std::upper_bound(x_.begin(), x_.end(), x);
  int i = static_cast<int>(it - x_.begin()) - 1;

  if (i < 0) i = 0;
  if (i > x_.size() - 2) i = x_.size() - 2;

  return i;
}

} // namespace WhetStone
} // namespace Amanzi

