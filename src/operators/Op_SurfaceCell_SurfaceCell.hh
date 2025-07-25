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

#ifndef AMANZI_OP_SURFACECELL_SURFACECELL_HH_
#define AMANZI_OP_SURFACECELL_SURFACECELL_HH_

#include <vector>
#include "DenseMatrix.hh"
#include "Operator.hh"
#include "Op_Cell_Cell.hh"

namespace Amanzi {
namespace Operators {

class Op_SurfaceCell_SurfaceCell : public Op_Cell_Cell {
 public:
  Op_SurfaceCell_SurfaceCell(const std::string& name,
                             const Teuchos::RCP<const AmanziMesh::Mesh> surf_mesh_)
    : Op_Cell_Cell(name, surf_mesh_, 1), surf_mesh(surf_mesh_)
  {}

  virtual void ApplyMatrixFreeOp(const Operator* assembler,
                                 const CompositeVector& X,
                                 CompositeVector& Y) const
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
    if (scaling.HasComponent("cell") &&
        scaling.ViewComponent("cell", false)->MyLength() ==
          surf_mesh->getNumEntities(AmanziMesh::Entity_kind::CELL,
                                    AmanziMesh::Parallel_kind::OWNED)) {
      Op_Cell_Cell::Rescale(scaling);
    }
    if (scaling.HasComponent("face") &&
        scaling.ViewComponent("face", false)->MyLength() ==
          mesh_->getNumEntities(AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_kind::OWNED)) {
      const Epetra_MultiVector& s_f = *scaling.ViewComponent("face", false);
      for (int k = 0; k != s_f.NumVectors(); ++k) {
        for (int sc = 0; sc != diag->MyLength(); ++sc) {
          auto f = surf_mesh->getEntityParent(AmanziMesh::Entity_kind::CELL, sc);
          (*diag)[k][sc] *= s_f[0][f];
        }
      }
    }
  }

 public:
  Teuchos::RCP<const AmanziMesh::Mesh> surf_mesh;
};

} // namespace Operators
} // namespace Amanzi


#endif
