/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  WhetStone, Version 2.2
  Release name: naka-to.

  The regularized basis for dG methods with monomials of the form
  (x - x0)^m / h^m. The transformation matrix R is diagonal with
  entries h^m.
*/

#ifndef AMANZI_DG_BASIS_REGULARIZED_HH_
#define AMANZI_DG_BASIS_REGULARIZED_HH_

#include <vector>

#include "Mesh.hh"

#include "Basis.hh"
#include "WhetStoneDefs.hh"

namespace Amanzi {
namespace WhetStone {

class Basis_Regularized : public Basis {
 public:
  Basis_Regularized() { id_ = TAYLOR_BASIS_REGULARIZED; }
  ~Basis_Regularized(){};

  // initialization
  virtual void Init(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
                    int c,
                    int order,
                    Polynomial<>& integrals) override;

  // transformation of bilinear form
  virtual void BilinearFormNaturalToMy(DenseMatrix<>& A) const override;
  virtual void BilinearFormNaturalToMy(std::shared_ptr<Basis> bl,
                                       std::shared_ptr<Basis> br,
                                       DenseMatrix<>& A) const override;

  // transformation of linear form
  virtual void LinearFormNaturalToMy(DenseVector<>& v) const override;

  // transformation of vector
  virtual void ChangeBasisMyToNatural(DenseVector<>& v) const override;
  virtual void ChangeBasisNaturalToMy(DenseVector<>& v) const override;

  // Recover polynomial in the natural basis
  virtual Polynomial<> CalculatePolynomial(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
                                         int c,
                                         int order,
                                         const DenseVector<>& coefs) const override;

  // access
  const AmanziMesh::Double_List& monomial_scales() const { return monomial_scales_; }

 private:
  using Basis::id_;
  AmanziMesh::Double_List monomial_scales_;

  int d_, order_;
};

} // namespace WhetStone
} // namespace Amanzi

#endif
