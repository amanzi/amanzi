/*
  WhetStone, Version 2.2
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Derham complex: mimetic inner products at nodes.
*/

#ifndef AMANZI_DERHAM_NODE_HH_
#define AMANZI_DERHAM_NODE_HH_

#include "Teuchos_RCP.hpp"

#include "Mesh.hh"

#include "DenseMatrix.hh"
#include "InnerProductL2.hh"
#include "Tensor.hh"

namespace Amanzi {
namespace WhetStone {

class DeRham_Node : virtual public InnerProductL2 {
 public:
  DeRham_Node(){};
  DeRham_Node(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh)
    : InnerProduct(mesh){};
  ~DeRham_Node(){};

  virtual int L2consistency(int c, const Tensor& T, DenseMatrix& N,
                            DenseMatrix& Mc, bool symmetry);
  virtual int MassMatrix(int c, const Tensor& T, DenseMatrix& M);
};

} // namespace WhetStone
} // namespace Amanzi

#endif
