/*
  WhetStone, Version 2.2
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Partial specialization of operations with vector of polynomials.
*/ 

#include <vector>

#include "Polynomial.hh"
#include "SpaceTimePolynomial.hh"
#include "VectorObjects.hh"

namespace Amanzi {
namespace WhetStone {

/* ******************************************************************
* Gradient for polynomials
****************************************************************** */
VectorPolynomial Gradient(const Polynomial& p)
{
  int d = p.dimension();
  int order = std::max(0, p.order() - 1);

  VectorPolynomial poly(d, d, order);
  poly.set_origin(p.get_origin());

  int index[3];
  for (auto it = p.begin(); it < p.end(); ++it) {
    int k = it.MonomialSetOrder();
    if (k > 0) {
      const int* idx = it.multi_index();
      int n = it.PolynomialPosition();
      double val = p(n);

      for (int i = 0; i < d; ++i) {
        for (int j = 0; j < d; ++j) index[j] = idx[j];

        if (index[i] > 0) {
          index[i]--;
          int m = MonomialSetPosition(d, index);
          poly[i](k - 1, m) = val * idx[i];
        }
      }
    }
  }

  return poly;
}


/* ******************************************************************
* Gradient for space-time polynomials
****************************************************************** */
VectorSpaceTimePolynomial Gradient(const SpaceTimePolynomial& p)
{
  int d = p.dimension();
  int order = p.order();

  VectorSpaceTimePolynomial poly(d, d, order);

  for (int n = 0; n < order + 1; ++n) {
    auto tmp = Gradient(p[n]);
    for (int i = 0; i < d; ++i) poly[i][n] = tmp[i];
  }

  return poly;
}


/* ******************************************************************
* Divergence
****************************************************************** */
Polynomial Divergence(const VectorPolynomial& vp) 
{
  int d = vp[0].dimension();
  AMANZI_ASSERT(d == vp.size());

  int order = vp[0].order();
  order = std::max(0, order - 1);

  Polynomial div(d, order);
  div.set_origin(vp[0].get_origin());

  int index[3];
  for (int i = 0; i < d; ++i) {
    for (auto it = vp[i].begin(); it < vp[i].end(); ++it) {
      int k = it.MonomialSetOrder();
      if (k > 0) {
        const int* idx = it.multi_index();
        for (int j = 0; j < d; ++j) index[j] = idx[j];

        if (index[i] > 0) {
          int n = it.PolynomialPosition();
          double val = vp[i](n);

          index[i]--;
          int m = MonomialSetPosition(d, index);
          div(k - 1, m) += val * idx[i];
        }
      }
    }
  }

  return div;
}


/* ******************************************************************
* Projecton of gradient using Taylor expansion with k terms
****************************************************************** */
VectorPolynomial GradientOnUnitSphere(const Polynomial& poly, int k)
{
  int d = poly.dimension();
  AMANZI_ASSERT(d == 2);
  AMANZI_ASSERT(k < 3);

  VectorPolynomial out(d, d, k);
  out.set_origin(poly.get_origin());

  double a1, a2, a3, a4, a5, a6, a7, a8, a9;
  double len, len2, len3, len5, ux, uy, vx, vy;
  a1 = poly(1);
  a2 = poly(2);
  len = std::pow(a1 * a1 + a2 * a2, 0.5);

  out[0](0) = a1 / len;
  out[1](0) = a2 / len;

  if (k > 0) {
    a3 = poly(3);
    a4 = poly(4);
    a5 = poly(5);

    double tmp1 = a1 * a4 - 2 * a2 * a3;
    double tmp2 = a2 * a4 - 2 * a1 * a5;

    len3 = len * len * len;
    ux =-a2 * tmp1;
    uy = a2 * tmp2;

    vx = a1 * tmp1;
    vy =-a1 * tmp2;

    out[0](1) = ux / len3; 
    out[0](2) = uy / len3; 

    out[1](1) = vx / len3; 
    out[1](2) = vy / len3; 
  }

  if (k > 1) {
    a6 = poly(6);
    a7 = poly(7);
    a8 = poly(8);
    a9 = poly(9);

    len2 = len * len;
    len5 = len3 * len2;
    double uxx = 2 * a2 * (3 * a2 * a6 + a3 * a4) - a1 * (a4 * a4 + 2 * a2 * a7);
    double uxy = 2 * a2 * (a2 * a7 + 4 * a3 * a5 - a1 * a8) - a4 * (a2 * a4 + 2 * a1 * a5);
    double uyy = 2 * a2 * (a2 * a8 + a4 * a5) - a1 * (4 * a5 * a5 + 6 * a2 * a9);

    double dx = 2 * a1 * a3 + a2 * a4;
    double dy = 2 * a2 * a5 + a1 * a4;

    double vxx = 2 * a1 * (a1 * a7 + a3 * a4) - a2 * (4 * a3 * a3 + 6 * a1 * a6); 
    double vxy = 2 * a1 * (a1 * a8 + a4 * a4 - 2 * a3 * a5) - 2 * a2 * (a3 * a4 + a1 * a7); 
    double vyy = 2 * a1 * (3 * a1 * a9 + a4 * a5) - a2 * (a4 * a4 + 2 * a1 * a8); 

    out[0](3) = (uxx * len2 - 3 * ux * dx) / (2 * len5);
    out[0](4) = (uxy * len2 - 3 * ux * dy) / len5;
    out[0](5) = (uyy * len2 - 3 * uy * dy) / (2 * len5);

    out[1](3) = (vxx * len2 - 3 * vx * dx) / (2 * len5);
    out[1](4) = (vxy * len2 - 3 * vx * dy) / len5;
    out[1](5) = (vyy * len2 - 3 * vy * dy) / (2 * len5);
  }

  return out;
}

}  // namespace WhetStone
}  // namespace Amanzi


