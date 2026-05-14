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

#include <iostream>

#include "SplineCubicNotAKnot1D.hh"
#include "SplineCubicNotAKnot2D.hh"

namespace Amanzi {
namespace WhetStone {

/* ******************************************************************
*
****************************************************************** */
int
SplineCubicNotAKnot2D::build(const std::vector<double>& x,
                             const std::vector<double>& y,
                             const std::vector<std::vector<double>>& values)
{
  x_ = x;
  y_ = y;

  int nx = x_.size();
  int ny = y_.size();
  if (nx < 4 || ny < 4) return -1;  // both grid directions need at least 4 nodes

  if (values.size() != ny) return -2;  // values is mis-sized
  for (int j = 0; j < ny; ++j) {
    if (values[j].size() != nx) return -2;
  }

  // bool ok1 = checkMonotonicity_(x_);
  // bool ok2 = checkMonotonicity_(y_);
  // if (!ok1 || !ok2) return -3; // invalid 2D mesh

  int ncells_x = nx - 1;

  // coeff_values[ix * 4 + p][j], where
  //   ix = x interval
  //   p = local x-polynomial coefficient 0,1,2,3
  //   j = y node
  std::vector<std::vector<double>> coeff_values(ncells_x * 4, std::vector<double>(ny));

  // spline in x for every fixed y_j
  for (int j = 0; j < ny; ++j) {
    std::vector<double> row(nx);

    for (int i = 0; i < nx; ++i) {
      row[i] = values[j][i];
    }

    SplineCubicNotAKnot1D sx(x_, row);

    for (int ix = 0; ix < ncells_x; ++ix) {
      const auto& c = sx.coeff(ix);

      coeff_values[4 * ix][j] = c.a;
      coeff_values[4 * ix + 1][j] = c.b;
      coeff_values[4 * ix + 2][j] = c.c;
      coeff_values[4 * ix + 3][j] = c.d;
    }
  }

  // spline those x-polynomial coefficients in y
  y_splines_.resize(ncells_x * 4);

  for (int ix = 0; ix < ncells_x; ++ix) {
    for (int p = 0; p < 4; ++p) {
      int idx = 4 * ix + p;
      y_splines_[idx].build(y_, coeff_values[idx]);
    }
  }

  return 0;
}


/* ******************************************************************
* Compute spline value and its derivatice, s, ds/dx,
* ds/dy, d2s/dx2, d2s/dx/dy, d2s/dy2
****************************************************************** */
std::array<double, 6>
SplineCubicNotAKnot2D::evaluate(double x, double y) const
{
  int ix = findInterval_(x_, x);
  double u = x - x_[ix];

  double C[4], Ct[4], Ctt[4];

  for (int p = 0; p < 4; ++p) {
    auto e = y_splines_[ix * 4 + p].evaluate(y);
    C[p] = e[0];
    Ct[p] = e[1];
    Ctt[p] = e[2];
  }

  double u2, u3, f, fd, ft, fdd, fdt, ftt;
  u2 = u * u;
  u3 = u2 * u;

  f = C[0] + C[1] * u + C[2] * u2 + C[3] * u3;
  fd = C[1] + 2.0 * C[2] * u + 3.0 * C[3] * u2;
  ft = Ct[0] + Ct[1] * u + Ct[2] * u2 + Ct[3] * u3;
  fdd = 2.0 * C[2] + 6.0 * C[3] * u;
  fdt = Ct[1] + 2.0 * Ct[2] * u + 3.0 * Ct[3] * u2;
  ftt = Ctt[0] + Ctt[1] * u + Ctt[2] * u2 + Ctt[3] * u3;

  return {f, fd, ft, fdd, fdt, ftt};
}

/* ******************************************************************
* 
****************************************************************** */
int 
SplineCubicNotAKnot2D::findInterval_(const std::vector<double>& x, double v) const
{
  if (v <= x.front()) return 0;
  if (v >= x.back()) return x.size() - 2;

  auto it = std::upper_bound(x.begin(), x.end(), v);
  int i = static_cast<int>(it - x.begin()) - 1;

  if (i < 0) i = 0;
  if (i > x.size() - 2) i = x.size() - 2;

  return i;
}

} // namespace WhetStone
} // namespace Amanzi

