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

  Derham complex: mimetic inner products at nodes.
*/

#ifndef AMANZI_DERHAM_NODE_HH_
#define AMANZI_DERHAM_NODE_HH_

#include "Teuchos_RCP.hpp"

#include "Mesh.hh"

#include "BilinearForm.hh"
#include "DenseMatrix.hh"
#include "MFD3D.hh"
#include "Tensor.hh"

namespace Amanzi {
namespace WhetStone {

class DeRham_Node : public MFD3D {
 public:
  DeRham_Node(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh) : MFD3D(mesh){};

  virtual std::vector<SchemaItem> schema() const override
  {
    return std::vector<SchemaItem>(
      1, std::make_tuple(AmanziMesh::Entity_kind::NODE, DOF_Type::SCALAR, 1));
  }

  int L2consistency(int c, const Tensor& T, DenseMatrix& N, DenseMatrix& Mc);

  virtual int MassMatrix(int c, const Tensor& T, DenseMatrix& M) override;
};

} // namespace WhetStone
} // namespace Amanzi

#endif
