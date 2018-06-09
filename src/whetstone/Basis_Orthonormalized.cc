/*
  WhetStone, version 2.1
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  The normalized basis for dG methods that is also orthogonal to 1.
*/

#include "Basis_Orthonormalized.hh"
#include "NumericalIntegration.hh"
#include "PolynomialIterator.hh"

namespace Amanzi {
namespace WhetStone {

/* ******************************************************************
* Prepare scaling data for the orthonormalized basis.
****************************************************************** */
void Basis_Orthonormalized::Init(
    const Teuchos::RCP<const AmanziMesh::Mesh>& mesh, int c, int order)
{
  int d = mesh->space_dimension();
  monomial_scales_.Reshape(d, order);
  monomial_ortho_.Reshape(d, order);

  Polynomial integrals(d, 2 * order);
  NumericalIntegration numi(mesh);

  for (int k = 0; k <= 2 * order; ++k) {
    numi.IntegrateMonomialsCell(c, k, integrals);
  }
  
  double volume = integrals(0, 0); 

  monomial_scales_(0, 0) = 1.0;
  monomial_ortho_(0, 0) = 0.0;

  for (auto it = monomial_scales_.begin(); it.end() <= monomial_scales_.end(); ++it) {
    int k = it.MonomialSetPosition();
    const int* multi_index = it.multi_index();
    int index[d]; 

    int m(0);
    for (int i = 0; i < d; ++i) {
      m += multi_index[i];
      index[i] = 2 * multi_index[i];
    }

    if (m > 0) {
      double b = integrals(m, k) / volume;
      double norm = integrals(2 * m, integrals.MonomialSetPosition(index));

      norm -= b * b * volume;
      monomial_scales_(m, k) = std::pow(volume / norm, 0.5);
      monomial_ortho_(m, k) = b;
    }
  }
}


/* ******************************************************************
* Transformation of regularized basis to owned basis.
****************************************************************** */
void Basis_Orthonormalized::ChangeBasisNaturalToMy(DenseMatrix& A) const
{
  AMANZI_ASSERT(A.NumRows() == monomial_scales_.size());

  int order = monomial_scales_.order();
  int d = monomial_scales_.dimension();

  int nrows = A.NumRows();
  std::vector<double> a(nrows), b(nrows);

  PolynomialIterator it(d);
  for (it.begin(); it.end() <= order; ++it) {
    int n = it.PolynomialPosition();
    int m = it.MonomialSetOrder();
    int k = it.MonomialSetPosition();

    double ak = monomial_scales_(m, k);
    double bk = monomial_ortho_(m, k);

    a[n] = ak;
    b[n] = -ak * bk;
  }

  // calculate A * R
  for (int k = 1; k < nrows; ++k) {
    for (int i = 0; i < nrows; ++i) {
      A(i, k) = A(i, k) * a[k] + A(i, 0) * b[k];
    }
  }

  // calculate R^T * A * R
  for (int k = 1; k < nrows; ++k) {
    for (int i = 0; i < nrows; ++i) {
      A(k, i) = A(k, i) * a[k] + A(0, i) * b[k];
    }
  }
}


/* ******************************************************************
* Vector transformation between my and natural bases.
* Transformation is bi-diagonal: Bn_i = a_i (B_i - b_i B_0) and 
* v = Sum_i v_i Bn_i. Note that B_0 = 1.
****************************************************************** */
void Basis_Orthonormalized::ChangeBasisMyToNatural(DenseVector& v) const
{
  AMANZI_ASSERT(v.NumRows() == monomial_scales_.size());
 
  for (auto it = monomial_scales_.begin(); it.end() <= monomial_scales_.end(); ++it) {
    int n = it.PolynomialPosition();
    int m = it.MonomialSetOrder();
    int k = it.MonomialSetPosition();

    double a = monomial_scales_(m, k);
    double b = monomial_ortho_(m, k);
    v(n) *= a;
    v(0) -= b * v(n);
  }
}


void Basis_Orthonormalized::ChangeBasisNaturalToMy(DenseVector& v) const
{
  AMANZI_ASSERT(v.NumRows() == monomial_scales_.size());
 
  double v0 = v(0);
  for (auto it = monomial_scales_.begin(); it.end() <= monomial_scales_.end(); ++it) {
    int n = it.PolynomialPosition();
    int m = it.MonomialSetOrder();
    int k = it.MonomialSetPosition();

    double a = monomial_scales_(m, k);
    double b = monomial_ortho_(m, k);
    v(n) /= a;
    v(0) += b * v(n);
  }
}


/* ******************************************************************
* Transformation from natural to my basis for polynomial coefficents.
****************************************************************** */
void Basis_Orthonormalized::ChangeBasisNaturalToMy(
    std::shared_ptr<Basis> bl, std::shared_ptr<Basis> br, DenseMatrix& A) const
{
  int order = monomial_scales_.order();
  int d = monomial_scales_.dimension();

  int nrows = A.NumRows();
  int m(nrows / 2);
  std::vector<double> a1(m), a2(m), b1(m), b2(m);

  auto bll = std::dynamic_pointer_cast<Basis_Orthonormalized>(bl);
  auto brr = std::dynamic_pointer_cast<Basis_Orthonormalized>(br);

  PolynomialIterator it(d);
  for (it.begin(); it.end() <= order; ++it) {
    int n = it.PolynomialPosition();
    int m = it.MonomialSetOrder();
    int k = it.MonomialSetPosition();

    double ak = (bll->monomial_scales())(m, k);
    double bk = (bll->monomial_ortho())(m, k);

    a1[n] = ak;
    b1[n] = -ak * bk;

    ak = (brr->monomial_scales())(m, k);
    bk = (brr->monomial_ortho())(m, k);
    a2[n] = ak;
    b2[n] = -ak * bk;
  }

  // calculate A * R
  for (int k = 1; k < m; ++k) {
    for (int i = 0; i < m; ++i) {
      A(i, k) = A(i, k) * a1[k] + A(i, 0) * b1[k];
      A(i, k + m) = A(i, k + m) * a2[k] + A(i, m) * b2[k];

      A(i + m, k) = A(i + m, k) * a1[k] + A(i + m, 0) * b1[k];
      A(i + m, k + m) = A(i + m, k + m) * a2[k] + A(i + m, m) * b2[k];
    }
  }

  // calculate R^T * A * R
  for (int k = 1; k < m; ++k) {
    for (int i = 0; i < m; ++i) {
      A(k, i) = A(k, i) * a1[k] + A(0, i) * b1[k];
      A(k + m, i) = A(k + m, i) * a2[k] + A(m, i) * b2[k];

      A(k, i + m) = A(k, i + m) * a1[k] + A(0, i + m) * b1[k];
      A(k + m, i + m) = A(k + m, i + m) * a2[k] + A(m, i + m) * b2[k];
    }
  }
}


/* ******************************************************************
* Recover polynomial from data coeffieints. 
****************************************************************** */
Polynomial Basis_Orthonormalized::CalculatePolynomial(
    const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
    int c, int order, DenseVector& coefs) const
{
  int d = mesh->space_dimension();
  Polynomial poly(monomial_scales_);
  poly.set_origin(mesh->cell_centroid(c));

  int n(0);
  for (auto it = poly.begin(); it.end() <= poly.end(); ++it) {
    int m = it.MonomialSetOrder();
    int k = it.MonomialSetPosition();

    poly(m, k) *= coefs(n++); 
    poly(0, 0) -= poly(m, k) * monomial_ortho_(m, k);
  }

  return poly;
}

}  // namespace WhetStone
}  // namespace Amanzi

