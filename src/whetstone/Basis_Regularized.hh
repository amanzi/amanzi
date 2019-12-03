/*
  WhetStone, Version 2.2
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  The regularized basis for dG methods with monomials of the form
  (x - x0)^m / h^m. The transformation matrix R is diagonal with
  entries h^m.
*/

#ifndef AMANZI_DG_BASIS_REGULARIZED_HH_
#define AMANZI_DG_BASIS_REGULARIZED_HH_

#include <vector>

#include "Basis.hh"
#include "WhetStoneDefs.hh"

namespace Amanzi {
namespace WhetStone {

template<class MyMesh>
class Basis_Regularized : public Basis<MyMesh> { 
 public:
  Basis_Regularized() { id_ = TAYLOR_BASIS_REGULARIZED; }
  ~Basis_Regularized() {};

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

  // special transformation of vectors for vector discretizations
  // it is assumed that v was assembled from same degree polynomials
  void ChangeBasisMyToNatural(DenseVector& v, int nrows) const;
  void ChangeBasisNaturalToMy(DenseVector& v, int nrows) const;

  // Recover polynomial in the natural basis
  virtual Polynomial CalculatePolynomial(const Teuchos::RCP<const MyMesh>& mesh,
                                         int c, int order, DenseVector& coefs) const;

  // access
  const std::vector<double>& monomial_scales() const { return monomial_scales_; }

 private:
  using Basis<MyMesh>::id_;
  std::vector<double> monomial_scales_;

  int d_, order_;
};


/* ******************************************************************
* Prepare scaling data for the regularized basis.
****************************************************************** */
template<class MyMesh>
void Basis_Regularized<MyMesh>::Init(
    const Teuchos::RCP<const MyMesh>& mesh,
    int c, int order, Polynomial& integrals)
{
  int dim, k0 = monomial_scales_.size();
  double volume;

  if (k0 < order + 1) {
    order_ = order;
    d_ = mesh->space_dimension();
    volume = mesh->cell_volume(c);

    monomial_scales_.resize(order + 1);
    for (int k = k0; k < order + 1; ++k) {
      monomial_scales_[k] = std::pow(volume, -(double)k / d_);
    }
  }
}


/* ******************************************************************
* Transformation from natural basis to my basis: A_new = R^T A_old R.
****************************************************************** */
template<class MyMesh>
void Basis_Regularized<MyMesh>::BilinearFormNaturalToMy(DenseMatrix& A) const
{
  int nrows = A.NumRows();
  std::vector<double> a(nrows);

  PolynomialIterator it(d_);
  for (it.begin(); it.MonomialSetOrder() <= order_; ++it) {
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
* Transformation from natural basis to my basis: f_new = R^T f_old.
****************************************************************** */
template<class MyMesh>
void Basis_Regularized<MyMesh>::LinearFormNaturalToMy(DenseVector& f) const
{
  PolynomialIterator it(d_);
  for (it.begin(); it.MonomialSetOrder() <= order_; ++it) {
    int n = it.PolynomialPosition();
    int m = it.MonomialSetOrder();
    f(n) *= monomial_scales_[m];
  }
}


/* ******************************************************************
* Transformation of interface matrix from natural to my bases.
****************************************************************** */
template<class MyMesh>
void Basis_Regularized<MyMesh>::BilinearFormNaturalToMy(
    std::shared_ptr<Basis<MyMesh> > bl,
    std::shared_ptr<Basis<MyMesh> > br, DenseMatrix& A) const
{
  int nrows = A.NumRows();
  int m(nrows / 2);
  std::vector<double> a1(m), a2(m);

  auto bll = std::dynamic_pointer_cast<Basis_Regularized>(bl);
  auto brr = std::dynamic_pointer_cast<Basis_Regularized>(br);

  PolynomialIterator it(d_);
  for (it.begin(); it.MonomialSetOrder() <= order_; ++it) {
    int n = it.PolynomialPosition();
    int m = it.MonomialSetOrder();

    a1[n] = (bll->monomial_scales())[m];
    a2[n] = (brr->monomial_scales())[m];
  }

  // calculate R^T * A * R
  for (int k = 0; k < m; ++k) {
    for (int i = 0; i < m; ++i) {
      A(i, k) = A(i, k) * a1[k] * a1[i];
      A(i, k + m) = A(i, k + m) * a1[i] * a2[k];

      A(i + m, k) = A(i + m, k) * a2[i] * a1[k];
      A(i + m, k + m) = A(i + m, k + m) * a2[i] * a2[k];
    }
  }
}


/* ******************************************************************
* Transformation from my to natural bases: v_old = R * v_new.
****************************************************************** */
template<class MyMesh>
void Basis_Regularized<MyMesh>::ChangeBasisMyToNatural(DenseVector& v) const
{
  PolynomialIterator it(d_);
  for (it.begin(); it.MonomialSetOrder() <= order_; ++it) {
    int n = it.PolynomialPosition();
    int m = it.MonomialSetOrder();
    v(n) *= monomial_scales_[m];
  }
}


/* ******************************************************************
* Transformation from natural to my bases: v_new = inv(R) * v_old.
****************************************************************** */
template<class MyMesh>
void Basis_Regularized<MyMesh>::ChangeBasisNaturalToMy(DenseVector& v) const
{
  PolynomialIterator it(d_);
  for (it.begin(); it.MonomialSetOrder() <= order_; ++it) {
    int n = it.PolynomialPosition();
    int m = it.MonomialSetOrder();
    v(n) /= monomial_scales_[m];
  }
}


/* ******************************************************************
* Special transformation for vector discretizations
****************************************************************** */
template<class MyMesh>
void Basis_Regularized<MyMesh>::ChangeBasisMyToNatural(DenseVector& v, int nrows) const
{
  int stride = v.NumRows() / nrows;

  PolynomialIterator it(d_);
  for (it.begin(); it.MonomialSetOrder() <= order_; ++it) {
    int n = it.PolynomialPosition();
    int m = it.MonomialSetOrder();
    for (int i = 0; i < nrows; ++i) {
      int k = i * stride + n;
      v(k) *= monomial_scales_[m];
    }
  }
}


/* ******************************************************************
* Special transformation for vector discretizations
****************************************************************** */
template<class MyMesh>
void Basis_Regularized<MyMesh>::ChangeBasisNaturalToMy(DenseVector& v, int nrows) const
{
  int stride = v.NumRows() / nrows;

  PolynomialIterator it(d_);
  for (it.begin(); it.MonomialSetOrder() <= order_; ++it) {
    int n = it.PolynomialPosition();
    int m = it.MonomialSetOrder();
    for (int i = 0; i < nrows; ++i) { 
      int k = i * stride + n;
      v(k) /= monomial_scales_[m];
    }
  }
}


/* ******************************************************************
* Recover polynomial in the natural basis from vector coefs of 
* coefficients in the regularized basis. 
****************************************************************** */
template<class MyMesh>
Polynomial Basis_Regularized<MyMesh>::CalculatePolynomial(
    const Teuchos::RCP<const MyMesh>& mesh,
    int c, int order, DenseVector& coefs) const
{
  int d = mesh->space_dimension();
  Polynomial poly(d, order, coefs);
  poly.set_origin(mesh->cell_centroid(c));

  int n(0);
  for (int k = 0; k < order + 1; ++k) {
    int mk = MonomialSpaceDimension(d, k);
    double scale = monomial_scales_[k];
    for (int m = 0; m < mk; ++m) poly(n++) *= scale;
  }

  return poly;
}

}  // namespace WhetStone
}  // namespace Amanzi

#endif
