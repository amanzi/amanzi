/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  Operators

*/

#ifndef AMANZI_OP_CELL_NODE_HH_
#define AMANZI_OP_CELL_NODE_HH_

#include <vector>
#include "DenseMatrix.hh"
#include "Operator.hh"
#include "Op.hh"

namespace Amanzi {
namespace Operators {

class Op_Cell_Node : public Op {
 public:
  Op_Cell_Node(const std::string& name, const Teuchos::RCP<const AmanziMesh::Mesh> mesh)
    : Op(OPERATOR_SCHEMA_BASE_CELL | OPERATOR_SCHEMA_DOFS_NODE, name, mesh)
  {
    WhetStone::DenseMatrix null_matrix;
    matrices.resize(mesh->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED),
                    null_matrix);
    matrices_shadow = matrices;
  }

  virtual void
  ApplyMatrixFreeOp(const Operator* assembler, const CompositeVector& X, CompositeVector& Y) const
  {
    assembler->ApplyMatrixFreeOp(*this, X, Y);
  }

  virtual void SymbolicAssembleMatrixOp(const Operator* assembler,
                                        const SuperMap& map,
                                        GraphFE& graph,
                                        int my_block_row,
                                        int my_block_col) const
  {
    assembler->SymbolicAssembleMatrixOp(*this, map, graph, my_block_row, my_block_col);
  }

  virtual void AssembleMatrixOp(const Operator* assembler,
                                const SuperMap& map,
                                MatrixFE& mat,
                                int my_block_row,
                                int my_block_col) const
  {
    assembler->AssembleMatrixOp(*this, map, mat, my_block_row, my_block_col);
  }

  // rescaling columns of local matrices
  virtual void Rescale(const CompositeVector& scaling)
  {
    if (scaling.HasComponent("node")) {
      const Epetra_MultiVector& s_n = *scaling.ViewComponent("node", true);
      AmanziMesh::Entity_ID_View nodes;

      for (int c = 0; c != matrices.size(); ++c) {
        WhetStone::DenseMatrix& Acell = matrices[c];
        nodes = mesh_->getCellNodes(c);

        for (int n = 0; n != nodes.size(); ++n) {
          for (int m = 0; m != nodes.size(); ++m) { Acell(n, m) *= s_n[0][nodes[n]]; }
        }
      }
    }
  }
};

} // namespace Operators
} // namespace Amanzi


#endif
