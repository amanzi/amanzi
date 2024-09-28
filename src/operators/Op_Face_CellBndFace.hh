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

#ifndef AMANZI_OP_FACE_CELLFACE_HH_
#define AMANZI_OP_FACE_CELLFACE_HH_

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

class Op_Face_CellBndFace : public Op {
 public:
  Op_Face_CellBndFace(std::string& name, const Teuchos::RCP<const AmanziMesh::Mesh> mesh)
    : Op(OPERATOR_SCHEMA_BASE_FACE | OPERATOR_SCHEMA_DOFS_CELL | OPERATOR_SCHEMA_DOFS_FACE |
           OPERATOR_SCHEMA_DOFS_BNDFACE,
         name,
         mesh)
  {
    WhetStone::DenseMatrix null_matrix;
    nfaces_owned =
      mesh->getNumEntities(AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_kind::OWNED);
    matrices.resize(nfaces_owned, null_matrix);
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

  virtual void Rescale(const CompositeVector& scaling)
  {
    if ((scaling.HasComponent("cell")) && (scaling.HasComponent("boundary_face"))) {
      const Epetra_MultiVector& s_c = *scaling.ViewComponent("cell", true);
      const Epetra_MultiVector& s_bnd = *scaling.ViewComponent("boundary_face");
      for (int f = 0; f != matrices.size(); ++f) {
        auto cells = mesh_->getFaceCells(f);
        if (cells.size() > 1) {
          matrices[f](0, 0) *= s_c[0][cells[0]];
          matrices[f](0, 1) *= s_c[0][cells[1]];
          matrices[f](1, 0) *= s_c[0][cells[0]];
          matrices[f](1, 1) *= s_c[0][cells[1]];
        } else if (cells.size() == 1) {
          int bf = mesh_->getMap(AmanziMesh::Entity_kind::BOUNDARY_FACE, false)
                     .LID(mesh_->getMap(AmanziMesh::Entity_kind::FACE, false).GID(f));
          matrices[f](0, 0) *= s_c[0][cells[0]];
          matrices[f](1, 0) *= s_c[0][cells[0]];
          matrices[f](0, 1) *= s_bnd[0][bf];
          matrices[f](1, 1) *= s_bnd[0][bf];
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
