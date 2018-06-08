/*
  WhetStone, version 2.1
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  The regularized basis for dG methods: x^k y^l / h^(k+l), where
  h is a measure of cell size.
*/

#include "Basis_Regularized.hh"

namespace Amanzi {
namespace WhetStone {

/* ******************************************************************
* Prepare scaling data for the regularized basis.
****************************************************************** */
void Basis_Regularized::Init(
    const Teuchos::RCP<const AmanziMesh::Mesh>& mesh, int c, int order)
{
  int k0 = monomial_scales_.size();

  if (k0 < order + 1) {
    order_ = order;
    d_ = mesh->space_dimension();
    double volume = mesh->cell_volume(c);

    monomial_scales_.resize(order + 1);
    for (int k = k0; k < order + 1; ++k) {
      monomial_scales_[k] = std::pow(volume, -(double)k / d_);
    }
  }
}


/* ******************************************************************
* Transformation from natural basis to owned basis.
****************************************************************** */
void Basis_Regularized::ChangeBasisMatrix(DenseMatrix& A) const
{
  AMANZI_ASSERT(A.NumRows() == monomial_scales_.size());

  int nrows = A.NumRows();
  std::vector<double> a(nrows);

  PolynomialIterator it(d_);
  for (it.begin(); it.end() <= order_; ++it) {
    int n = it.PolynomialPosition();
    int m = it.MonomialSetOrder();
    a[n] = monomial_scales_[m];
  }

  // calculate R^T * A * R
  for (int k = 0; k < nrows; ++k) {
    for (int i = 0; i < nrows; ++i) {
      A(i, k) = A(i, k) * a[k] * a[i];
    }
  }
}


/* ******************************************************************
* Transformation from natural basis to owned basis.
****************************************************************** */
void Basis_Regularized::ChangeBasisVector(DenseVector& v) const
{
  AMANZI_ASSERT(v.NumRows() == monomial_scales_.size());

  PolynomialIterator it(d_);
  for (it.begin(); it.end() <= order_; ++it) {
    int n = it.PolynomialPosition();
    int m = it.MonomialSetOrder();
    v(n) /= monomial_scales_[m];
  }
}


/* ******************************************************************
* Transformation from natural basis to owned basis.
****************************************************************** */
void Basis_Regularized::ChangeBasisMatrix(
    std::shared_ptr<Basis> bl, std::shared_ptr<Basis> br, DenseMatrix& A) const
{
  int nrows = A.NumRows();
  int m(nrows / 2);
  std::vector<double> a1(m), a2(m);

  auto bll = std::dynamic_pointer_cast<Basis_Regularized>(bl);
  auto brr = std::dynamic_pointer_cast<Basis_Regularized>(br);

  PolynomialIterator it(d_);
  for (it.begin(); it.end() <= order_; ++it) {
    int n = it.PolynomialPosition();
    int m = it.MonomialSetOrder();

    a1[n] = (bll->monomial_scales())[m];
    a2[n] = (brr->monomial_scales())[m];
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
Polynomial Basis_Regularized::CalculatePolynomial(
    const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
    int c, int order, DenseVector& coefs) const
{
  int d = mesh->space_dimension();
  Polynomial poly(d, order);

  poly.SetPolynomialCoefficients(coefs);
  poly.set_origin(mesh->cell_centroid(c));

  for (int k = 0; k < order + 1; ++k) {
    auto& mono = poly.MonomialSet(k);
    mono *= monomial_scales_[k];
  }

  return poly;
}

}  // namespace WhetStone
}  // namespace Amanzi

