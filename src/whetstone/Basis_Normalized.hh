/*
  WhetStone, Version 2.2
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  The normalized basis for dG methods with monomials of the form
  |a (x-x0)^m| = 1. The transformation matrix R is diagonal with
  entries a.
*/

#ifndef AMANZI_DG_BASIS_NORMALIZED_HH_
#define AMANZI_DG_BASIS_NORMALIZED_HH_

#include "Basis.hh"
#include "NumericalIntegration.hh"
#include "Polynomial.hh"
#include "WhetStoneDefs.hh"

namespace Amanzi {
namespace WhetStone {

template <class MyMesh>
class Basis_Normalized : public Basis<MyMesh> { 
 public:
  Basis_Normalized() { id_ = TAYLOR_BASIS_NORMALIZED; }
  ~Basis_Normalized() {};

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

  // recover polynomial in the natural basis
  virtual Polynomial CalculatePolynomial(const Teuchos::RCP<const MyMesh>& mesh,
                                         int c, int order, DenseVector& coefs) const;

  // access
  const Polynomial& monomial_scales() const { return monomial_scales_; }

 private:
  using Basis<MyMesh>::id_;
  Polynomial monomial_scales_;
};


/* ******************************************************************
* Prepare scaling data for the normalized basis.
****************************************************************** */
template<class MyMesh>
void Basis_Normalized<MyMesh>::Init(
    const Teuchos::RCP<const MyMesh>& mesh,
    int c, int order, Polynomial& integrals)
{
  int d = mesh->space_dimension();
  monomial_scales_.Reshape(d, order);

  NumericalIntegration numi(mesh);
  numi.UpdateMonomialIntegralsCell(c, 2 * order, integrals);
  
  double volume = integrals(0); 

  monomial_scales_(0) = 1.0;
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
      double norm = integrals(2 * m, MonomialSetPosition(d, index));
      monomial_scales_(m, k) = std::pow(volume / norm, 0.5);
    }
  }
}


/* ******************************************************************
* Transformation from natural basis to my basis: A_new = R^T A_old R.
****************************************************************** */
template<class MyMesh>
void Basis_Normalized<MyMesh>::BilinearFormNaturalToMy(DenseMatrix& A) const
{
  AMANZI_ASSERT(A.NumRows() == monomial_scales_.size());

  int nrows = A.NumRows();
  const DenseVector& a = monomial_scales_.coefs();

  // calculate R^T * A * R
  for (int k = 0; k < nrows; ++k) {
    for (int i = 0; i < nrows; ++i) {
      A(i, k) = A(i, k) * a(k) * a(i);
    }
  }
}


/* ******************************************************************
* Transformation of interface matrix from natural to my bases.
****************************************************************** */
template<class MyMesh>
void Basis_Normalized<MyMesh>::BilinearFormNaturalToMy(
    std::shared_ptr<Basis<MyMesh> > bl,
    std::shared_ptr<Basis<MyMesh> > br, DenseMatrix& A) const
{
  int nrows = A.NumRows();
  int m(nrows / 2);

  auto bll = std::dynamic_pointer_cast<Basis_Normalized>(bl);
  auto brr = std::dynamic_pointer_cast<Basis_Normalized>(br);

  const DenseVector& a1 = bll->monomial_scales().coefs();
  const DenseVector& a2 = brr->monomial_scales().coefs();

  // calculate R^T * A * R
  for (int k = 0; k < m; ++k) {
    for (int i = 0; i < m; ++i) {
      A(i, k) = A(i, k) * a1(i) * a1(k);
      A(i, k + m) = A(i, k + m) * a1(i) * a2(k);

      A(i + m, k) = A(i + m, k) * a2(i) * a1(k);
      A(i + m, k + m) = A(i + m, k + m) * a2(i) * a2(k);
    }
  }
}


/* ******************************************************************
* Transformation from natural basis to my basis: f_new = R^T f_old.
****************************************************************** */
template<class MyMesh>
void Basis_Normalized<MyMesh>::LinearFormNaturalToMy(DenseVector& f) const
{
  for (int n = 0; n < monomial_scales_.size(); ++n) {
    f(n) *= monomial_scales_(n);
  }
}


/* ******************************************************************
* Transformation from my to natural bases: v_old = R * v_new.
****************************************************************** */
template<class MyMesh>
void Basis_Normalized<MyMesh>::ChangeBasisMyToNatural(DenseVector& v) const
{
  AMANZI_ASSERT(v.NumRows() == monomial_scales_.size());

  for (int n = 0; n < monomial_scales_.size(); ++n) {
    v(n) *= monomial_scales_(n);
  }
}


/* ******************************************************************
* Transformation from natural to my bases: v_new = inv(R) * v_old.
****************************************************************** */
template<class MyMesh>
void Basis_Normalized<MyMesh>::ChangeBasisNaturalToMy(DenseVector& v) const
{
  AMANZI_ASSERT(v.NumRows() == monomial_scales_.size());

  for (int n = 0; n < monomial_scales_.size(); ++n) {
    v(n) /= monomial_scales_(n);
  }
}


/* ******************************************************************
* Recover polynomial in the natural basis from vector coefs of 
* coefficients in the normalized basis. 
****************************************************************** */
template<class MyMesh>
Polynomial Basis_Normalized<MyMesh>::CalculatePolynomial(
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
  }

  return poly;
}

}  // namespace WhetStone
}  // namespace Amanzi

#endif
