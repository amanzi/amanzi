/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon (coonet@ornl.gov)
*/

//! <MISSING_ONELINE_DOCSTRING>

#ifndef AMANZI_OP_FACE_CELL_HH_
#define AMANZI_OP_FACE_CELL_HH_

#include <vector>
#include "DenseMatrix.hh"
#include "Operator.hh"
#include "Op.hh"

/*
  Op classes are small structs that play two roles:

  1. They provide a class name to the schema, enabling visitor patterns.
  2. They are a container for local matrices.

  This Op class is for storing local matrices of length nfaces and with dofs
  on cells, i.e. for Advection or for TPFA.
*/

namespace Amanzi {
namespace Operators {

class Op_Face_Cell : public Op {
 public:
  Op_Face_Cell(const std::string& name,
               const Teuchos::RCP<const AmanziMesh::Mesh> mesh)
    : Op(OPERATOR_SCHEMA_BASE_FACE | OPERATOR_SCHEMA_DOFS_CELL, name, mesh)
  {
    WhetStone::DenseMatrix null_matrix;
    nfaces_owned =
      mesh->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED);
    matrices.resize(nfaces_owned, null_matrix);
    matrices_shadow = matrices;
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
    if (scaling.HasComponent("cell")) {
      const Epetra_MultiVector& s_c = *scaling.ViewComponent("cell", true);
      AmanziMesh::Entity_ID_List cells;
      for (int f = 0; f != matrices.size(); ++f) {
        mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
        matrices[f](0, 0) *= s_c[0][cells[0]];
        if (cells.size() > 1) {
          matrices[f](0, 1) *= s_c[0][cells[1]];
          matrices[f](1, 0) *= s_c[0][cells[0]];
          matrices[f](1, 1) *= s_c[0][cells[1]];
        }
      }
    }
  }

 protected:
  int nfaces_owned;
};

} // namespace Operators
} // namespace Amanzi


#endif
