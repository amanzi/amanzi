/*
  WhetStone, Version 2.2
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Miscalleneous tools for working with vector objects.
*/ 

#include "Monomial.hh"
#include "VectorObjects.hh"
#include "VectorObjectsUtils.hh"

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
      if (val == 0.0) continue;

      for (int i = 0; i < d; ++i) {
        for (int j = 0; j < d; ++j) index[j] = idx[j];

        if (index[i] > 0) {
          index[i]--;
          int m = PolynomialPosition(d, index);
          poly[i](m) = val * idx[i];
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
* Curl in 3D: vector function -> vector function
****************************************************************** */
VectorPolynomial Curl3D(const VectorPolynomial& vp)
{
  int d = vp[0].dimension();
  AMANZI_ASSERT(d == vp.NumRows() && d == 3);

  VectorPolynomial curl(d, d, 0);
  curl.set_origin(vp[0].get_origin());

  for (int i = 0; i < d; ++i) {
    int j = (i + 1) % d;
    int k = (i + 2) % d;

    int order = std::max(vp[j].order(), vp[k].order());
    order = std::max(0, order - 1);
    curl[i].Reshape(d, order, true);
  }

  int index[3];
  for (int i = 0; i < d; ++i) {
    int j = (i + 1) % d;
    int k = (i + 2) % d;

    for (auto it = vp[i].begin(); it < vp[i].end(); ++it) {
      int m = it.MonomialSetOrder();
      if (m > 0) {
        const int* idx = it.multi_index();
        for (int l = 0; l < d; ++l) index[l] = idx[l];

        // j-component is update by d u_i / d x_k
        if (index[k] > 0) {
          int n = it.PolynomialPosition();

          index[k]--;
          int pos = MonomialSetPosition(d, index);
          curl[j](m - 1, pos) += vp[i](n) * idx[k];
          
          index[k]++;
        }

        // k-component is update by -d u_i / d x_j
        if (index[j] > 0) {
          int n = it.PolynomialPosition();

          index[j]--;
          int pos = MonomialSetPosition(d, index);
          curl[k](m - 1, pos) -= vp[i](n) * idx[j];
        }
      }
    }
  }

  return curl;
}


/* ******************************************************************
* Curl in 2D: vector function -> scalar function
****************************************************************** */
Polynomial Curl2D(const VectorPolynomial& p)
{
  int d = p[0].dimension();
  int m = p.NumRows();
  AMANZI_ASSERT(d == m && d == 2);

  // gradient of vector polynomial: poor efficiency FIXME
  std::vector<VectorPolynomial> grad(d);
  for (int i = 0; i < d; ++i) { 
    grad[i] = Gradient(p[i]);
  }

  // assemble curl
  Polynomial tmp(d, 0);
  tmp = grad[1][0] - grad[0][1];

  return tmp;
}


/* ******************************************************************
* Rot in 2D: scalar function -> vector function
****************************************************************** */
VectorPolynomial Rot2D(const Polynomial& p)
{
  int d = p.dimension();
  AMANZI_ASSERT(d == 2);

  auto grad = Gradient(p);

  VectorPolynomial tmp(d, d, 0);
  tmp[0] = grad[1];
  tmp[1] = grad[0];
  tmp[1] *= -1.0;

  return tmp;
}


/* ******************************************************************
* Divergence
****************************************************************** */
Polynomial Divergence(const VectorPolynomial& vp) 
{
  int d = vp[0].dimension();
  AMANZI_ASSERT(d == vp.NumRows());

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
* Matrix of a 3D curl operator
****************************************************************** */
DenseMatrix Curl3DMatrix(int d, int order)
{
  // dimensions of the curl-matrix
  int nd0 = PolynomialSpaceDimension(d, order - 1);
  int nd1 = PolynomialSpaceDimension(d, order);

  int nrows(0), ncols(0);
  int offset_rows[3], offset_cols[3];

  for (int i = 0; i < d; ++i) {
    offset_rows[i] = nrows;
    offset_cols[i] = ncols;

    nrows += nd0;
    ncols += nd1;
  }

  DenseMatrix curl(nrows, ncols);
  curl.PutScalar(0.0);

  // populate with cofficients for partial derivatives
  PolynomialIterator it0(d), it1(d);
  it0.begin(0);
  it1.begin(order + 1);

  int index[3];
  for (int i = 0; i < d; ++i) {
    int j = (i + 1) % d;
    int k = (i + 2) % d;

    for (auto it = it0; it < it1; ++it) {
      int m = it.PolynomialPosition();
      if (m > 0) {
        int col = offset_cols[i] + m;

        const int* idx = it.multi_index();
        for (int l = 0; l < d; ++l) index[l] = idx[l];

        // j-component is update by d u_i / d x_k
        if (index[k] > 0) {
          index[k]--;
          int row = offset_rows[j] + PolynomialPosition(d, index);
          curl(row, col) += idx[k];
          index[k]++;
        }

        // k-component is update by -d u_i / d x_j
        if (index[j] > 0) {
          index[j]--;
          int row = offset_rows[k] + PolynomialPosition(d, index);
          curl(row, col) -= idx[j];
        }
      }
    }
  }

  return curl;
}


/* ******************************************************************
* Vector is created by setting one component to *this.
*   3D: q_k = curl(p_k ^ x) + x . p_{k-1}
****************************************************************** */
void VectorDecomposition3DCurl(const Monomial& q, int component,
                               VectorPolynomial& p1, Polynomial& p2)
{
  int d = q.dimension();
  const int* index = q.multi_index();

  double a(2.0);
  for (int i = 0; i < d; ++i) a += index[i];

  double coef = q.coefs()(0);
  p1.Reshape(d, d, 0, true);
  p1[component] = Polynomial(d, index, coef / a);
  p1.set_origin(q.get_origin());

  int idx[3];
  if (index[component] > 0) {
    for (int i = 0; i < d; ++i) idx[i] = index[i];
    idx[component]--;
    p2 = Polynomial(d, idx, coef * index[component] / a);
  } else {
    p2.Reshape(d, 0, true);
  }
  p2.set_origin(q.get_origin());
}


/* ******************************************************************
* 2D vector decomposition: q_k = tor(p_{k+1}) + x . p_{k-1}
****************************************************************** */
void VectorDecomposition2DRot(
    const VectorPolynomial& q, Polynomial& p1, Polynomial& p2)
{
  // reshape output
  int d = q[0].dimension();
  int order(0);
  for (int k = 0; k < d; ++k) order = std::max(order, q[k].order() - 1);

  p1.Reshape(d, order + 2, true);
  p1.set_origin(q[0].get_origin());

  p2.Reshape(d, order, true);
  p2.set_origin(q[0].get_origin());

  // calculate decomposition for each monomial of each component
  int idx[3];
  for (int k = 0; k < d; ++k) {
    for (auto it = q[k].begin(); it < q[k].end(); ++it) {
      int n = it.PolynomialPosition();
      const int* index = it.multi_index();

      double a(1.0);
      for (int i = 0; i < d; ++i) {
        a += index[i];
        idx[i] = index[i];
      }

      idx[1 - k]++;
      int l = PolynomialPosition(d, idx);

      double coef = q[k](n);
      if (k == 0)
        p1(l) += coef / a;
      else
        p1(l) -= coef / a;
      
      if (index[k] > 0) {
        idx[1 - k]--;
        idx[k]--;
        
        l = PolynomialPosition(d, idx);
        p2(l) += coef * index[k] / a;
      }
    }
  }
}


/* ******************************************************************
* 3D vector decomposition: q_k = grad(p_{k+1}) + x ^ p_{k-1}.
* Returns position of nonzero in monomial (!) p1.
****************************************************************** */
int VectorDecomposition3DGrad(const Monomial& q, int component,
                              Polynomial& p1, VectorPolynomial& p2)
{
  int d = q.dimension();
  const int* index = q.multi_index();
  double coef = q.coefs()(0);

  int idx[3], a(1);
  for (int i = 0; i < d; ++i) {
    idx[i] = index[i];
    a += idx[i];
  }

  idx[component]++;
  p1 = Polynomial(d, idx, coef / a);
  p1.set_origin(q.get_origin());
  int pos = PolynomialPosition(d, idx);

  p2.Reshape(d, d, 0, true);
  idx[component]--;
  
  int k1 = (component + 1) % d;
  int k2 = (component + 2) % d;

  if (idx[k2] > 0) {
    double tmp = (idx[k2] * coef) / a;
    idx[k2]--;
    p2[k1] = Polynomial(d, idx, -tmp);
    idx[k2]++;
  }

  if (idx[k1] > 0) {
    double tmp = (idx[k1] * coef) / a;
    idx[k1]--;
    p2[k2] = Polynomial(d, idx, tmp);
  }
  p2.set_origin(q.get_origin());

  return pos;
}


/* ******************************************************************
* Project a vector on manifold.
****************************************************************** */
VectorPolynomial ProjectVectorPolynomialOnManifold(
   const VectorPolynomial& vpoly, const AmanziGeometry::Point& x0,
   const std::vector<AmanziGeometry::Point>& tau)
{
  int d = tau.size();
  VectorPolynomial tmp(d, d, 0);    

  for (int i = 0; i < d; ++i) {
    tmp[i] = vpoly * tau[i];
    tmp[i].ChangeCoordinates(x0, tau);
  }

  return tmp;
}


/* ******************************************************************
* Return vector of coefficients. Blocks of fixed size are alloced
* for vector components which is due to discretization requirements.
****************************************************************** */
DenseVector ExpandCoefficients(VectorPolynomial& vp)
{
  int m(0);
  for (int k= 0; k < vp.NumRows(); ++k) {
    m = std::max(m, vp[k].size());
  }
  int mrows = m * vp.NumRows();

  DenseVector tmp(mrows);
  double* vdata = tmp.Values();

  for (int k= 0; k < vp.NumRows(); ++k) {
    const double* data = vp[k].coefs().Values();
    int size = vp[k].size();
    for (int i = 0; i < size; ++i) {
      *vdata = *data;
      vdata++;
      data++;
    }
    for (int i = size; i < m; ++i) {
      *vdata = 0.0;
      vdata++;
    }
  }

  return tmp;
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


