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
               const Teuchos::RCP<const AmanziMesh::Mesh> surf_mesh_) :
      Op_Face_Cell(name, surf_mesh_)
  {};

  virtual void ApplyMatrixFreeOp(const Operator* assembler,
          const CompositeVector& X, CompositeVector& Y) const {
    assembler->ApplyMatrixFreeOp(*this, X, Y);
  }

  virtual void
  GetLocalDiagCopy(CompositeVector& X) const
  {
    auto Xv = X.ViewComponent("face", true);
    
    AmanziMesh::Mesh const * surf_mesh = mesh.get();
    // entity_get_parent is not available on device
    for(int sf = 0 ; sf < csr.size(); ++sf){
      AmanziMesh::Entity_ID_View cells;
      surf_mesh->face_get_cells(sf, AmanziMesh::Parallel_type::ALL, cells);
      auto f0 = surf_mesh->entity_get_parent(AmanziMesh::CELL, cells(0));
      WhetStone::DenseMatrix lm(csr.at(sf),csr.size(sf,0),csr.size(sf,1)); 
      Kokkos::atomic_add(&Xv(f0,0), lm(0,0));
      if (cells.extent(0) > 1) {
        auto f1 = surf_mesh->entity_get_parent(AmanziMesh::CELL, cells(1));
        Kokkos::atomic_add(&Xv(f1,0), lm(1,1));
      }
    }
  }
  
  //virtual void SymbolicAssembleMatrixOp(const Operator* assembler,
  //        const SuperMap& map, GraphFE& graph,
  //        int my_block_row, int my_block_col) const {
  //  assembler->SymbolicAssembleMatrixOp(*this, map, graph, my_block_row, my_block_col);
  //}

  //virtual void AssembleMatrixOp(const Operator* assembler,
  //        const SuperMap& map, MatrixFE& mat,
  //        int my_block_row, int my_block_col) const {
  //  assembler->AssembleMatrixOp(*this, map, mat, my_block_row, my_block_col);
  //}
  
  virtual void Rescale(const CompositeVector& scaling) {
    if (scaling.HasComponent("cell") &&
        scaling.ViewComponent("cell", false).extent(1) ==
        mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED)) {
      // scaling's cell entry is defined on the surface mesh
      Op_Face_Cell::Rescale(scaling);
    } 

    if (scaling.HasComponent("face") &&
        scaling.GetComponent("face", false)->getLocalLength() ==
        mesh->parent()->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED)) {
      Exceptions::amanzi_throw("Scaling surface cell entities with subsurface face vector not yet implemented");
    }
  }


 public:
  Teuchos::RCP<const AmanziMesh::Mesh> surf_mesh;
};

}  // namespace Operators
}  // namespace Amanzi


#endif

