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

#include "MeshLight.hh"

#include "Basis.hh"
#include "NumericalIntegration.hh"
#include "Polynomial.hh"
#include "WhetStoneDefs.hh"

namespace Amanzi {
namespace WhetStone {

class Basis_Normalized : public Basis { 
 public:
  Basis_Normalized() { id_ = TAYLOR_BASIS_NORMALIZED; }
  ~Basis_Normalized() {};

  // initialization
  virtual void Init(const Teuchos::RCP<const AmanziMesh::MeshLight>& mesh,
                    int c, int order, Polynomial& integrals);

  // transformation of bilinear form
  virtual void BilinearFormNaturalToMy(DenseMatrix& A) const;
  virtual void BilinearFormNaturalToMy(std::shared_ptr<Basis> bl,
                                       std::shared_ptr<Basis> br, DenseMatrix& A) const;

  // transformation of linear form
  virtual void LinearFormNaturalToMy(DenseVector& v) const;

  // transformation of vector 
  virtual void ChangeBasisMyToNatural(DenseVector& v) const;
  virtual void ChangeBasisNaturalToMy(DenseVector& v) const;

  // recover polynomial in the natural basis
  virtual Polynomial CalculatePolynomial(const Teuchos::RCP<const AmanziMesh::MeshLight>& mesh,
                                         int c, int order, DenseVector& coefs) const;

  // access
  const Polynomial& monomial_scales() const { return monomial_scales_; }

 private:
  using Basis::id_;
  Polynomial monomial_scales_;
};

}  // namespace WhetStone
}  // namespace Amanzi

#endif
