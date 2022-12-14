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

#ifndef AMANZI_OP_SURFACEFACE_SURFACECELL_HH_
#define AMANZI_OP_SURFACEFACE_SURFACECELL_HH_

#include <vector>
#include "DenseMatrix.hh"
#include "Operator.hh"
#include "Op_Face_Cell.hh"

namespace Amanzi {
namespace Operators {

class Op_SurfaceFace_SurfaceCell : public Op_Face_Cell {
 public:
  Op_SurfaceFace_SurfaceCell(const std::string& name,
                             const Teuchos::RCP<const AmanziMesh::Mesh> surf_mesh_)
    : Op_Face_Cell(name, surf_mesh_), surf_mesh(surf_mesh_){};

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
    if (scaling.HasComponent("cell") &&
        scaling.ViewComponent("cell", false)->MyLength() ==
          surf_mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED)) {
      // scaling's cell entry is defined on the surface mesh
      Op_Face_Cell::Rescale(scaling);
    }

    if (scaling.HasComponent("face") &&
        scaling.ViewComponent("face", false)->MyLength() ==
          mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED)) {
      AMANZI_ASSERT(mesh_ != surf_mesh);
      Exceptions::amanzi_throw(
        "Scaling surface cell entities with subsurface face vector not yet implemented");
    }
  }


 public:
  Teuchos::RCP<const AmanziMesh::Mesh> surf_mesh;
};

} // namespace Operators
} // namespace Amanzi


#endif
