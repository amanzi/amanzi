/*
  WhetStone, Version 2.2
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

#include <vector>

#include "Teuchos_RCP.hpp"

#include "Mesh.hh"

#include "DenseMatrix.hh"
#include "Tensor.hh"
#include "WhetStoneDefs.hh"

namespace Amanzi {
namespace WhetStone {

class InnerProduct {
 public:
  InnerProduct(){};
  InnerProduct(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh)
    : mesh_(mesh), stability_method_(WHETSTONE_STABILITY_GENERIC)
  {
    d_ = mesh_->space_dimension();
  }
  ~InnerProduct(){};

  // access
  double scalar_stability() { return scalar_stability_; }
  double scaling_factor() { return scaling_factor_; }

 protected:
  // supporting stability methods
  void StabilityScalar_(DenseMatrix& N, DenseMatrix& M);
  int StabilityOptimized_(const Tensor& T, DenseMatrix& N, DenseMatrix& M);

  double CalculateStabilityScalar_(DenseMatrix& Mc);
  void GrammSchmidt_(DenseMatrix& N);

 protected:
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  int d_;

  int stability_method_; // stability parameters
  double scalar_stability_, scaling_factor_;
};

} // namespace WhetStone
} // namespace Amanzi

#endif
