/*
  WhetStone, version 2.1
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Base class for mimetic inner products.
*/

#ifndef AMANZI_INNER_PRODUCT_HH_
#define AMANZI_INNER_PRODUCT_HH_

#include "Teuchos_RCP.hpp"

#include "Mesh.hh"

#include "DenseMatrix.hh"
#include "Tensor.hh"
#include "WhetStoneDefs.hh"

namespace Amanzi {
namespace WhetStone {

class InnerProduct { 
 public:
  explicit InnerProduct(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh) : mesh_(mesh) {};
  ~InnerProduct() {};

 protected:
  void StabilityScalar_(DenseMatrix& N, DenseMatrix& M);
  double CalculateStabilityScalar_(DenseMatrix& Mc);
  void GrammSchmidt_(DenseMatrix& N);

 protected:
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  double scaling_factor_;
};

}  // namespace WhetStone
}  // namespace Amanzi

#endif

