/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_OP_CELL_FACE_HH_
#define AMANZI_OP_CELL_FACE_HH_

#include <vector>
#include "DenseMatrix.hh"
#include "Operator.hh"
#include "Op.hh"

namespace Amanzi {
namespace Operators {

class Op_Cell_Face : public Op {
 public:
  Op_Cell_Face(const std::string& name,
               const Teuchos::RCP<const AmanziMesh::Mesh> mesh) :
      Op(OPERATOR_SCHEMA_BASE_CELL |
         OPERATOR_SCHEMA_DOFS_FACE, name, mesh) {
    WhetStone::DenseMatrix null_matrix;
    matrices.resize(mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED), null_matrix);
    matrices_shadow = matrices;
  }

  virtual void ApplyMatrixFreeOp(const Operator* assembler,
          const CompositeVector& X, CompositeVector& Y) const {
    assembler->ApplyMatrixFreeOp(*this, X, Y);
  }

  virtual void SymbolicAssembleMatrixOp(const Operator* assembler,
          const SuperMap& map, GraphFE& graph,
          int my_block_row, int my_block_col) const {
    assembler->SymbolicAssembleMatrixOp(*this,
            map, graph, my_block_row, my_block_col);
  }

  virtual void AssembleMatrixOp(const Operator* assembler,
          const SuperMap& map, MatrixFE& mat,
          int my_block_row, int my_block_col) const {
    assembler->AssembleMatrixOp(*this, map, mat, my_block_row, my_block_col);
  }

  virtual void Rescale(const CompositeVector& scaling) {
    if (scaling.HasComponent("face")) {
      const Epetra_MultiVector& s_f = *scaling.ViewComponent("face",true);
      for (int c = 0; c != matrices.size(); ++c) {
        const auto& faces = mesh_->cell_get_faces(c);
        for (int n = 0; n != faces.size(); ++n) {
          for (int m = 0; m != faces.size(); ++m) {
            matrices[c](n,m) *= s_f[0][faces[n]];
          }
        }
      }
    }
  }
};

}  // namespace Operators
}  // namespace Amanzi


#endif


