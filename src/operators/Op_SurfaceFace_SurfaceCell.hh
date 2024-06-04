/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@lanl.gov)
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
    : Op_Face_Cell(name, surf_mesh_){};

  virtual void
  ApplyMatrixFreeOp(const Operator* assembler, const CompositeVector& X, CompositeVector& Y) const
  {
    assembler->ApplyMatrixFreeOp(*this, X, Y);
  }

  virtual void SumLocalDiag(CompositeVector& X) const
  {
    AMANZI_ASSERT(false); // this is incorrectly implemented... --ETC
    // AmanziMesh::Mesh const* mesh_ = mesh.get();
    // AmanziMesh::Mesh const* surf_mesh_ = surf_mesh.get();
    // auto Xv = X.viewComponent<Amanzi::MirrorHost>("face", true);

    // Kokkos::parallel_for(
    //   "Op_SurfaceFace_SurfaceCell::SumLocalDiag", A.size(), KOKKOS_LAMBDA(const int& sf) {
    //     auto cells = mesh_->getFaceCells(sf);

    //     auto lm = A[sf];
    //     auto f0 = surf_mesh_->getEntityParent(AmanziMesh::CELL, cells(0));
    //     Kokkos::atomic_add(&Xv(f0, 0), lm(0, 0));
    //     if (cells.extent(0) > 1) {
    //       auto f1 = surf_mesh_->getEntityParent(AmanziMesh::CELL, cells(1));
    //       Kokkos::atomic_add(&Xv(f1, 0), lm(1, 1));
    //     }
    //   });
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
    if (scaling.hasComponent("cell") &&
        scaling.getComponent("cell", false)->getLocalLength() ==
          mesh->getNumEntities(AmanziMesh::CELL, AmanziMesh::Parallel_kind::OWNED)) {
      // scaling's cell entry is defined on the surface mesh
      Op_Face_Cell::Rescale(scaling);
    }

    if (scaling.hasComponent("face") && scaling.getComponent("face", false)->getLocalLength() ==
                                          mesh->getParentMesh()->getNumEntities(
                                            AmanziMesh::FACE, AmanziMesh::Parallel_kind::OWNED)) {
      Exceptions::amanzi_throw(
        "Scaling surface cell entities with subsurface face vector not yet implemented");
    }
  }
};

} // namespace Operators
} // namespace Amanzi


#endif
