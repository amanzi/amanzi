/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon (coonet@ornl.gov)
*/

//! <MISSING_ONELINE_DOCSTRING>

#ifndef AMANZI_OP_NODE_NODE_HH_
#define AMANZI_OP_NODE_NODE_HH_

#include <vector>
#include "Operator.hh"
#include "Op.hh"

namespace Amanzi {
namespace Operators {

class Op_Node_Node : public Op {
 public:
  Op_Node_Node(const std::string& name,
               const Teuchos::RCP<const AmanziMesh::Mesh> mesh, int nvec)
    : Op(OPERATOR_SCHEMA_BASE_NODE | OPERATOR_SCHEMA_DOFS_NODE, name, mesh)
  {
    // TODO add diag?
    //diag = Teuchos::rcp(new Epetra_MultiVector(mesh->node_map(false), nvec));
    //diag_shadow =
    //  Teuchos::rcp(new Epetra_MultiVector(mesh->node_map(false), nvec));
    Kokkos::resize(diags,mesh->node_map(false)->getNodeNumElements (), nvec);
    Kokkos::resize(diags_shadow,mesh->node_map(false)->getNodeNumElements (), nvec);  
  }

  virtual void
  ApplyMatrixFreeOp(const Operator* assembler, const CompositeVector& X,
                    CompositeVector& Y) const
  {
    assembler->ApplyMatrixFreeOp(*this, X, Y);
  }

  virtual void
  ApplyTransposeMatrixFreeOp(const Operator* assembler,
                             const CompositeVector& X, CompositeVector& Y) const
  {
    assembler->ApplyMatrixFreeOp(*this, X, Y);
  }

  virtual void
  SymbolicAssembleMatrixOp(const Operator* assembler, const SuperMap& map,
                           GraphFE& graph, int my_block_row,
                           int my_block_col) const
  {
    assembler->SymbolicAssembleMatrixOp(
      *this, map, graph, my_block_row, my_block_col);
  }

  virtual void
  AssembleMatrixOp(const Operator* assembler, const SuperMap& map,
                   MatrixFE& mat, int my_block_row, int my_block_col) const
  {
    assembler->AssembleMatrixOp(*this, map, mat, my_block_row, my_block_col);
  }

  virtual void Rescale(const CompositeVector& scaling)
  {
    if (scaling.HasComponent("node")) {
      const auto& s_v = scaling.ViewComponent("node", false);
      for (int k = 0; k != s_v.extent(0); ++k) {
        for (int i = 0; i != s_v.extent(1); ++i) {
          diags(k,i) *= s_v(0,i);
        }
      }
    }
  }
};

} // namespace Operators
} // namespace Amanzi


#endif
