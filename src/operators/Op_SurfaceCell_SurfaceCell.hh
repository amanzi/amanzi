/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@lanl.gov)
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
               const Teuchos::RCP<const AmanziMesh::Mesh> surf_mesh_) :
      Op_Cell_Cell(name, surf_mesh_)
  {}

  virtual void
  SumLocalDiag(CompositeVector& X) const
  {
    auto Xv = X.ViewComponent("face", false);
    auto diag_v = diag->getLocalViewHost();

    const AmanziMesh::Mesh* mesh_ = mesh.get();
    Kokkos::parallel_for(
        "Op_Surfacecell_Surfacecell::SumLocalDiag",
        diag_v.extent(0),
        KOKKOS_LAMBDA(const int& sc) {
          auto f = mesh_->entity_get_parent(AmanziMesh::CELL, sc);
          Xv(f,0) += diag_v(sc,0);
        });
  }

  virtual void ApplyMatrixFreeOp(const Operator* assembler,
          const CompositeVector& X, CompositeVector& Y) const {
    assembler->ApplyMatrixFreeOp(*this, X, Y);
  }

  virtual void SymbolicAssembleMatrixOp(const Operator* assembler,
         const SuperMap& map, GraphFE& graph,
         int my_block_row, int my_block_col) const {
   assembler->SymbolicAssembleMatrixOp(*this, map, graph, my_block_row, my_block_col);
  }

  virtual void AssembleMatrixOp(const Operator* assembler,
         const SuperMap& map, MatrixFE& mat,
         int my_block_row, int my_block_col) const {
   assembler->AssembleMatrixOp(*this, map, mat, my_block_row, my_block_col);
  }
  
  virtual void Rescale(const CompositeVector& scaling) {
    if (scaling.HasComponent("cell") &&
        scaling.GetComponent("cell", false)->getLocalLength() ==
        mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED)) {
      Op_Cell_Cell::Rescale(scaling);
    }

    if (scaling.HasComponent("face") &&
        scaling.GetComponent("face", false)->getLocalLength() ==
        mesh->parent()->num_entities(AmanziMesh::FACE,
                AmanziMesh::Parallel_type::OWNED)) {

      const auto s_f = scaling.ViewComponent("face", false);
      auto diag_v = diag->getLocalViewDevice();
      const AmanziMesh::Mesh* mesh_ = mesh.get();

      Kokkos::parallel_for(
          "Op_SurfaceCell_SurfaceCell::Rescale",
          diag_v.extent(0),
          KOKKOS_LAMBDA(const int& sc) {
            auto f = mesh_->entity_get_parent(AmanziMesh::CELL, sc);
            diag_v(sc,0) *= s_f(f,0);
          });
    }
  }
};

}  // namespace Operators
}  // namespace Amanzi


#endif

