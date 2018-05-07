/*
  WhetStone, version 2.1
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  The normalized basis for dG methods.
*/

#include "Basis_Normalized.hh"
#include "NumericalIntegration.hh"

namespace Amanzi {
namespace WhetStone {

/* ******************************************************************
* Prepare scaling data for the natural basis.
****************************************************************** */
void Basis_Normalized::Init(
    const Teuchos::RCP<const AmanziMesh::Mesh>& mesh, int c, int order)
{
  int d = mesh->space_dimension();
  monomial_scales_.Reshape(d, order);

  Polynomial integrals(d, 2 * order);
  NumericalIntegration numi(mesh, true);

  for (int k = 0; k <= 2 * order; ++k) {
    numi.IntegrateMonomialsCell(c, integrals.monomials(k));
  }
  
  double volume = integrals(0, 0); 

  monomial_scales_(0, 0) = 1.0;
  for (auto it = monomial_scales_.begin(); it.end() <= monomial_scales_.end(); ++it) {
    int k = it.MonomialPosition();
    const int* multi_index = it.multi_index();
    int index[d]; 

    int m(0);
    for (int i = 0; i < d; ++i) {
      m += multi_index[i];
      index[i] = 2 * multi_index[i];
    }

    if (m > 0) {
      const auto& aux2 = integrals.monomials(2 * m).coefs();
      double norm = aux2[integrals.MonomialPosition(index)];
      monomial_scales_(m, k) = std::pow(volume / norm, 0.5);
    }
  }
}


/* ******************************************************************
* Transfrmation of natural basis to owned basis.
****************************************************************** */
void Basis_Normalized::ChangeBasisMatrix(DenseMatrix& A) const
{
  ASSERT(A.NumRows() == monomial_scales_.size());

  int order = monomial_scales_.order();
  int d = monomial_scales_.dimension();

  int nrows = A.NumRows();
  std::vector<double> a(nrows);

  Iterator it(d);
  for (it.begin(); it.end() <= order; ++it) {
    int m = it.MonomialOrder();
    int k = it.MonomialPosition();

    a[k] = monomial_scales_(m, k);
  }

  // calculate R^T * A * R
  for (int k = 0; k < nrows; ++k) {
    for (int i = 0; i < nrows; ++i) {
      A(i, k) = A(i, k) * a[k] * a[i];
    }
  }
}


/* ******************************************************************
* Transfrmation of natural basis to owned basis.
****************************************************************** */
void Basis_Normalized::ChangeBasisVector(DenseVector& v) const
{
  ASSERT(v.NumRows() == monomial_scales_.size());

  for (auto it = monomial_scales_.begin(); it.end() <= monomial_scales_.end(); ++it) {
    int n = it.PolynomialPosition();
    int m = it.MonomialOrder();
    int k = it.MonomialPosition();

    v(n) /= monomial_scales_(m, k);
  }
}


/* ******************************************************************
* Transfrmation of natural basis to owned basis.
****************************************************************** */
void Basis_Normalized::ChangeBasisMatrix(
    std::shared_ptr<Basis> bl, std::shared_ptr<Basis> br, DenseMatrix& A) const
{
  int order = monomial_scales_.order();
  int d = monomial_scales_.dimension();

  int nrows = A.NumRows();
  int m(nrows / 2);
  std::vector<double> a1(m), a2(m);

  auto bll = std::dynamic_pointer_cast<Basis_Normalized>(bl);
  auto brr = std::dynamic_pointer_cast<Basis_Normalized>(br);

  Iterator it(d);
  for (it.begin(); it.end() <= order; ++it) {
    int n = it.PolynomialPosition();
    int m = it.MonomialOrder();
    int k = it.MonomialPosition();

    double ak = (bll->monomial_scales())(m, k);
    a1[n] = ak;

    ak = (brr->monomial_scales())(m, k);
    a2[n] = ak;
  }

  // calculate R^T * A * R
  for (int k = 0; k < m; ++k) {
    for (int i = 0; i < m; ++i) {
      A(i, k) = A(i, k) * a1[k] * a2[i];
      A(i, k + m) = A(i, k + m) * a2[k] * a1[i];

      A(i + m, k) = A(i + m, k) * a1[k] * a2[i];
      A(i + m, k + m) = A(i + m, k + m) * a2[k] * a1[i];
    }
  }
}


/* ******************************************************************
* Recover polynomial from data coeffieints. 
****************************************************************** */
Polynomial Basis_Normalized::CalculatePolynomial(
    const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
    int c, int order, DenseVector& coefs) const
{
  int d = mesh->space_dimension();
  Polynomial poly(monomial_scales_);
  poly.set_origin(mesh->cell_centroid(c));

  int n(0);
  for (auto it = poly.begin(); it.end() <= poly.end(); ++it) {
    int m = it.MonomialOrder();
    int k = it.MonomialPosition();

    poly(m, k) *= coefs(n++); 
  }

  return poly;
}

}  // namespace WhetStone
}  // namespace Amanzi

