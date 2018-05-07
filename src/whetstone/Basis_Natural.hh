/*
  WhetStone, version 2.1
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  The natural basis for dG methods: (x - x0)^m / h^m.
*/

#ifndef AMANZI_DG_BASIS_NATURAL_HH_
#define AMANZI_DG_BASIS_NATURAL_HH_

#include <vector>

#include "Basis.hh"
#include "WhetStoneDefs.hh"

namespace Amanzi {
namespace WhetStone {

class Basis_Natural : public Basis { 
 public:
  Basis_Natural() { id_ = TAYLOR_BASIS_NATURAL; }
  ~Basis_Natural() {};

  // initialization
  virtual void Init(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh, int c, int order);

  // transformation from natural basis to owned basis (nothing to do)
  virtual void ChangeBasisMatrix(DenseMatrix& A) const {};
  virtual void ChangeBasisVector(DenseVector& v) const {};

  virtual void ChangeBasisMatrix(Basis* bl, Basis* br, DenseMatrix& A) const {};
  virtual void ChangeBasisMatrix(std::shared_ptr<Basis> bl, std::shared_ptr<Basis> br, DenseMatrix& A) const {};

  // Recover polynomial in regular basis
  virtual Polynomial CalculatePolynomial(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
                                         int c, int order, DenseVector& coefs) const;

  // access
  const std::vector<double>& monomial_scales() { return monomial_scales_; }

 private:
  using Basis::id_;
  std::vector<double> monomial_scales_;
};

}  // namespace WhetStone
}  // namespace Amanzi

#endif

