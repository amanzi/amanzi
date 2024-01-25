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

  The normalized basis for dG methods with monomials of the form
  |a (x-x0)^m| = 1. The transformation matrix R is diagonal with
  entries a.
*/

#ifndef AMANZI_DG_BASIS_NORMALIZED_HH_
#define AMANZI_DG_BASIS_NORMALIZED_HH_

#include "Mesh.hh"

#include "Basis.hh"
#include "NumericalIntegration.hh"
#include "Polynomial.hh"
#include "WhetStoneDefs.hh"

namespace Amanzi {
namespace WhetStone {

class Basis_Normalized : public Basis {
 public:
  Basis_Normalized() { id_ = TAYLOR_BASIS_NORMALIZED; }
  ~Basis_Normalized(){};

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

  // recover polynomial in the natural basis
  virtual Polynomial<> CalculatePolynomial(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
                                         int c,
                                         int order,
                                         const DenseVector<>& coefs) const override;

  // access
  const Polynomial<>& monomial_scales() const { return monomial_scales_; }

 private:
  using Basis::id_;
  Polynomial<> monomial_scales_;
};

} // namespace WhetStone
} // namespace Amanzi

#endif
