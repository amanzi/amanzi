/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon (coonet@ornl.gov)
*/

//! <MISSING_ONELINE_DOCSTRING>

#ifndef AMANZI_OP_SURFACECELL_SURFACECELL_HH_
#define AMANZI_OP_SURFACECELL_SURFACECELL_HH_

#include "Operator.hh"
#include "Op_Cell_Cell.hh"

namespace Amanzi {
namespace Operators {

class Op_SurfaceCell_SurfaceCell : public Op_Cell_Cell {
 public:
  Op_SurfaceCell_SurfaceCell(const std::string& name,
                             const Teuchos::RCP<const AmanziMesh::Mesh> surf_mesh_)
    : Op_Cell_Cell(name, surf_mesh_),
      surf_mesh(surf_mesh_)
  {}

  virtual void
  getLocalDiagCopy(CompositeVector& X) const
  {
    auto Xv = X.ViewComponent("face", false);
    // entity_get_parent is not available on device
    for(int sc = 0 ; sc < data.extent(0); ++sc){
      auto f = surf_mesh->entity_get_parent(AmanziMesh::CELL, sc);
      Xv(f,0) += data(sc,0);
    }
  }
  
  virtual void
  ApplyMatrixFreeOp(const Operator* assembler, const CompositeVector& X,
                    CompositeVector& Y) const
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
    if (scaling.HasComponent("cell") &&
        scaling.GetComponent("cell", false)->getLocalLength() ==
          surf_mesh->num_entities(AmanziMesh::CELL,
                                  AmanziMesh::Parallel_type::OWNED)) {
      Op_Cell_Cell::Rescale(scaling);
    }
    if (scaling.HasComponent("face") &&
        scaling.GetComponent("face", false)->getLocalLength() ==
          mesh->num_entities(AmanziMesh::FACE,
                              AmanziMesh::Parallel_type::OWNED)) {

      auto scaling_v = scaling.ViewComponent("face", false);

      for(int sc = 0 ; sc < scaling_v.extent(0); ++sc){
        auto f = surf_mesh->entity_get_parent(AmanziMesh::CELL, sc);
        data(sc, 0) *= scaling_v(f, 0);
      }
    }
  }

 public:
  Teuchos::RCP<const AmanziMesh::Mesh> surf_mesh;
};

} // namespace Operators
} // namespace Amanzi


#endif
