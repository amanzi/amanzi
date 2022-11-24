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

#include "MeshLight.hh"

#include "Basis.hh"
#include "NumericalIntegration.hh"
#include "Polynomial.hh"
#include "WhetStoneDefs.hh"

namespace Amanzi {
namespace WhetStone {

class Basis_Orthonormalized : public Basis {
 public:
  Basis_Orthonormalized() { id_ = TAYLOR_BASIS_NORMALIZED_ORTHO; }
  ~Basis_Orthonormalized(){};

  // initialization
  virtual void Init(const Teuchos::RCP<const AmanziMesh::MeshLight>& mesh,
                    int c,
                    int order,
                    Polynomial& integrals);

  // transformation of bilinear form
  virtual void BilinearFormNaturalToMy(DenseMatrix& A) const;
  virtual void BilinearFormNaturalToMy(std::shared_ptr<Basis> bl,
                                       std::shared_ptr<Basis> br,
                                       DenseMatrix& A) const;

  // transformation of linear form
  virtual void LinearFormNaturalToMy(DenseVector& v) const;

  // transformation of vector
  virtual void ChangeBasisMyToNatural(DenseVector& v) const;
  virtual void ChangeBasisNaturalToMy(DenseVector& v) const;

  // Recover polynomial in the natural basis
  virtual Polynomial CalculatePolynomial(const Teuchos::RCP<const AmanziMesh::MeshLight>& mesh,
                                         int c,
                                         int order,
                                         DenseVector& coefs) const;

  // access
  const Polynomial& monomial_scales() const { return monomial_scales_; }
  const Polynomial& monomial_ortho() const { return monomial_ortho_; }

 private:
  using Basis::id_;
  Polynomial monomial_scales_, monomial_ortho_;
};

} // namespace WhetStone
} // namespace Amanzi

#endif
