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

#ifndef AMANZI_OP_CELL_FACECELL_HH_
#define AMANZI_OP_CELL_FACECELL_HH_

#include <vector>
#include "DenseMatrix.hh"
#include "Operator.hh"
#include "Op.hh"

namespace Amanzi {
namespace Operators {

class Op_Cell_FaceCell : public Op {
 public:
  Op_Cell_FaceCell(const std::string& name, const Teuchos::RCP<const AmanziMesh::Mesh> mesh)
    : Op(OPERATOR_SCHEMA_BASE_CELL | OPERATOR_SCHEMA_DOFS_FACE | OPERATOR_SCHEMA_DOFS_CELL,
         name,
         mesh)
  {
    WhetStone::DenseMatrix null_matrix;
    matrices.resize(
      mesh->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED),
      null_matrix);
    matrices_shadow = matrices;
  }

  virtual Teuchos::RCP<Op> DeepClone() const override
  {
    auto op = Teuchos::rcp(new Op_Cell_FaceCell(*this));
    *op->diag = *diag;
    *op->diag_shadow = *diag_shadow;
    return op;
  }

  virtual void ApplyMatrixFreeOp(const Operator* assembler,
                                 const CompositeVector& X,
                                 CompositeVector& Y) const override
  {
    assembler->ApplyMatrixFreeOp(*this, X, Y);
  }

  virtual void SymbolicAssembleMatrixOp(const Operator* assembler,
                                        const SuperMap& map,
                                        GraphFE& graph,
                                        int my_block_row,
                                        int my_block_col) const override
  {
    assembler->SymbolicAssembleMatrixOp(*this, map, graph, my_block_row, my_block_col);
  }

  virtual void AssembleMatrixOp(const Operator* assembler,
                                const SuperMap& map,
                                MatrixFE& mat,
                                int my_block_row,
                                int my_block_col) const override
  {
    assembler->AssembleMatrixOp(*this, map, mat, my_block_row, my_block_col);
  }

  virtual void Rescale(const CompositeVector& scaling) override
  {
    if (scaling.HasComponent("face")) {
      const Epetra_MultiVector& s_c = *scaling.ViewComponent("face", true);
      for (int c = 0; c != matrices.size(); ++c) {
        const auto& faces = mesh_->getCellFaces(c);
        for (int n = 0; n != faces.size(); ++n) {
          for (int m = 0; m != faces.size(); ++m) {
            matrices[c](n, m) *= s_c[0][faces[n]];
          }
          matrices[c](n, faces.size()) *= s_c[0][faces[n]];
        }
      }
    }

    if (scaling.HasComponent("cell")) {
      const Epetra_MultiVector& s_c = *scaling.ViewComponent("cell", true);
      AMANZI_ASSERT(s_c.MyLength() == matrices.size());

      AmanziMesh::Entity_ID_View face;
      for (int c = 0; c != matrices.size(); ++c) {
        int nfaces = mesh_->getCellNumFaces(c);
        for (int m = 0; m != nfaces; ++m) {
          matrices[c](nfaces, m) *= s_c[0][c];
        }
        matrices[c](nfaces, nfaces) *= s_c[0][c];
      }
    }
  }
};

} // namespace Operators
} // namespace Amanzi


#endif
