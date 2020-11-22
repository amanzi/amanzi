/*
  WhetStone, Version 2.2
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  The partially orthonormalized basis for dG methods: |\psi| = 1 and
  (\psi, 1) = 0 for \psi != 1. Transformation matrix has the form
      | 1 -b1 -b2 -b3 -b4 |
      |    a1             |
  R = |        a2         |
      |            a3     |
      |               a4  |
*/

#ifndef AMANZI_DG_BASIS_ORTHONORMALIZED_HH_
#define AMANZI_DG_BASIS_ORTHONORMALIZED_HH_

#include "Basis.hh"
#include "NumericalIntegration.hh"
#include "Polynomial.hh"
#include "WhetStoneDefs.hh"

namespace Amanzi {
namespace WhetStone {

template<class MyMesh>
class Basis_Orthonormalized : public Basis<MyMesh> { 
 public:
  Basis_Orthonormalized() { id_ = TAYLOR_BASIS_NORMALIZED_ORTHO; }
  ~Basis_Orthonormalized() {};

  // initialization
  virtual void Init(const Teuchos::RCP<const MyMesh>& mesh,
                    int c, int order, Polynomial& integrals);

  // transformation of bilinear form
  virtual void BilinearFormNaturalToMy(DenseMatrix& A) const;
  virtual void BilinearFormNaturalToMy(std::shared_ptr<Basis<MyMesh> > bl,
                                       std::shared_ptr<Basis<MyMesh> > br, DenseMatrix& A) const;

  // transformation of linear form
  virtual void LinearFormNaturalToMy(DenseVector& v) const;

  // transformation of vector 
  virtual void ChangeBasisMyToNatural(DenseVector& v) const;
  virtual void ChangeBasisNaturalToMy(DenseVector& v) const;

  // Recover polynomial in the natural basis
  virtual Polynomial CalculatePolynomial(const Teuchos::RCP<const MyMesh>& mesh,
                                         int c, int order, DenseVector& coefs) const;

  // access
  const Polynomial& monomial_scales() const { return monomial_scales_; }
  const Polynomial& monomial_ortho() const { return monomial_ortho_; }

 private:
  using Basis<MyMesh>::id_;
  Polynomial monomial_scales_, monomial_ortho_;
};


/* ******************************************************************
* Prepare scaling data for the orthonormalized basis.
****************************************************************** */
template<class MyMesh>
void Basis_Orthonormalized<MyMesh>::Init(
    const Teuchos::RCP<const MyMesh>& mesh,
    int c, int order, Polynomial& integrals)
{
  int d = mesh->space_dimension();
  monomial_scales_.Reshape(d, order);
  monomial_ortho_.Reshape(d, order);

  NumericalIntegration<AmanziMesh::MeshLight> numi(mesh);
  numi.UpdateMonomialIntegralsCell(c, 2 * order, integrals);
  
  double volume = integrals(0); 

  monomial_scales_(0) = 1.0;
  monomial_ortho_(0) = 0.0;

  for (auto it = monomial_scales_.begin(); it < monomial_scales_.end(); ++it) {
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
      double norm = integrals(2 * m, MonomialSetPosition(d, index));

      norm -= b * b * volume;
      monomial_scales_(m, k) = std::pow(volume / norm, 0.5);
      monomial_ortho_(m, k) = b;
    }
  }
}


/* ******************************************************************
* Transformation from natural basis to my basis: A_new = R^T A_old R.
****************************************************************** */
template<class MyMesh>
void Basis_Orthonormalized<MyMesh>::BilinearFormNaturalToMy(DenseMatrix& A) const
{
  AMANZI_ASSERT(A.NumRows() == monomial_scales_.size());

  int nrows = A.NumRows();
  std::vector<double> a(nrows), b(nrows);

  for (int n = 0; n < nrows; ++n) {
    double ak = monomial_scales_(n);
    double bk = monomial_ortho_(n);

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
* Transformation of interface matrix from natural to my bases.
****************************************************************** */
template<class MyMesh>
void Basis_Orthonormalized<MyMesh>::BilinearFormNaturalToMy(
    std::shared_ptr<Basis<MyMesh> > bl,
    std::shared_ptr<Basis<MyMesh> > br, DenseMatrix& A) const
{
  int nrows = A.NumRows();
  int m(nrows / 2);
  std::vector<double> a1(m), a2(m), b1(m), b2(m);

  auto bll = std::dynamic_pointer_cast<Basis_Orthonormalized>(bl);
  auto brr = std::dynamic_pointer_cast<Basis_Orthonormalized>(br);

  for (auto it = bll->monomial_scales().begin(); it < bll->monomial_scales().end(); ++it) {
    int n = it.PolynomialPosition();
    int s = it.MonomialSetOrder();
    int k = it.MonomialSetPosition();

    double ak = (bll->monomial_scales())(s, k);
    double bk = (bll->monomial_ortho())(s, k);

    a1[n] = ak;
    b1[n] = -ak * bk;

    ak = (brr->monomial_scales())(s, k);
    bk = (brr->monomial_ortho())(s, k);
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
* Transformation from natural basis to my basis: f_new = R^T f_old.
****************************************************************** */
template<class MyMesh>
void Basis_Orthonormalized<MyMesh>::LinearFormNaturalToMy(DenseVector& f) const
{
  int nrows = f.NumRows();
  std::vector<double> a(nrows), b(nrows);

  for (auto it = monomial_scales_.begin(); it < monomial_scales_.end(); ++it) {
    int n = it.PolynomialPosition();
    int m = it.MonomialSetOrder();
    int k = it.MonomialSetPosition();

    double ak = monomial_scales_(m, k);
    double bk = monomial_ortho_(m, k);

    a[n] = ak;
    b[n] = -ak * bk;
  }

  // calculate R^T * f
  for (int k = 1; k < nrows; ++k) {
    f(k) = f(k) * a[k] + f(0) * b[k];
  }
}


/* ******************************************************************
* Transformation from my to natural bases: v_old = R * v_new.
****************************************************************** */
template<class MyMesh>
void Basis_Orthonormalized<MyMesh>::ChangeBasisMyToNatural(DenseVector& v) const
{
  AMANZI_ASSERT(v.NumRows() == monomial_scales_.size());
 
  for (auto it = monomial_scales_.begin(); it < monomial_scales_.end(); ++it) {
    int n = it.PolynomialPosition();
    int m = it.MonomialSetOrder();
    int k = it.MonomialSetPosition();

    double a = monomial_scales_(m, k);
    double b = monomial_ortho_(m, k);
    v(n) *= a;
    v(0) -= b * v(n);
  }
}


/* ******************************************************************
* Transformation from natural to my bases: v_new = inv(R) * v_old.
****************************************************************** */
template<class MyMesh>
void Basis_Orthonormalized<MyMesh>::ChangeBasisNaturalToMy(DenseVector& v) const
{
  AMANZI_ASSERT(v.NumRows() == monomial_scales_.size());
 
  for (auto it = monomial_scales_.begin(); it < monomial_scales_.end(); ++it) {
    int n = it.PolynomialPosition();
    int m = it.MonomialSetOrder();
    int k = it.MonomialSetPosition();

    double a = monomial_scales_(m, k);
    double b = monomial_ortho_(m, k);
    v(0) += b * v(n);
    v(n) /= a;
  }
}


/* ******************************************************************
* Recover polynomial in the natural basis from vector coefs of 
* coefficients in the orthonormalized basis. 
****************************************************************** */
template<class MyMesh>
Polynomial Basis_Orthonormalized<MyMesh>::CalculatePolynomial(
    const Teuchos::RCP<const MyMesh>& mesh,
    int c, int order, DenseVector& coefs) const
{
  Polynomial poly(monomial_scales_);
  poly.set_origin(mesh->cell_centroid(c));

  int n(0);
  for (auto it = poly.begin(); it < poly.end(); ++it) {
    int m = it.MonomialSetOrder();
    int k = it.MonomialSetPosition();

    poly(m, k) *= coefs(n++); 
    poly(0, 0) -= poly(m, k) * monomial_ortho_(m, k);
  }

  return poly;
}

}  // namespace WhetStone
}  // namespace Amanzi

#endif

