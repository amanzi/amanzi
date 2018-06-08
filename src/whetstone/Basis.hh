/*
  WhetStone, version 2.1
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  The base class for dG basis.
  We highlight two special bases: natural, with basis monomials 
  of form (x-x0)^k, and regularized, with basis monomials of 
  form (x-x0)^k / h^k, where h is cell measure. Matrix and vector
  Tranformations are available to and from these bases.
*/

#ifndef AMANZI_DG_BASIS_HH_
#define AMANZI_DG_BASIS_HH_

#include "Teuchos_RCP.hpp"

#include "Mesh.hh"

#include "DenseMatrix.hh"
#include "DenseVector.hh"
#include "Polynomial.hh"

namespace Amanzi {
namespace WhetStone {

class Basis { 
 public:
  Basis() {};
  virtual ~Basis() {};

  // initialization
  virtual void Init(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh, int c, int order) = 0;

  // transformation from natural basis to owned basis
  virtual void ChangeBasisMatrix(DenseMatrix& A) const = 0;
  virtual void ChangeBasisVector(DenseVector& v) const = 0;

  virtual void ChangeBasisMatrix(std::shared_ptr<Basis> bl, std::shared_ptr<Basis> br, DenseMatrix& A) const = 0;

  // recover polynomial in natural basis
  virtual Polynomial CalculatePolynomial(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
                                         int c, int order, DenseVector& coefs) const = 0;

  // assess 
  int id() { return id_; };

 protected:
  int id_;
};

}  // namespace WhetStone
}  // namespace Amanzi

#endif

