/*
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@lanl.gov)
*/

//! Op on all FACES with matrices for CELL-adjacencies.

#ifndef AMANZI_OP_FACE_CELL_HH_
#define AMANZI_OP_FACE_CELL_HH_

#include <vector>
#include "DenseMatrix.hh"
#include "Operator.hh"
#include "Op.hh"

namespace Amanzi {
namespace Operators {

class Op_Face_Cell : public Op {
 public:
  Op_Face_Cell(const std::string& name,
               const Teuchos::RCP<const AmanziMesh::Mesh> mesh) :
      Op(OPERATOR_SCHEMA_BASE_FACE |
         OPERATOR_SCHEMA_DOFS_CELL, name, mesh) {
    int nfaces_owned = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED);
    matrices.resize(nfaces_owned);

    for (int f=0; f!=nfaces_owned; ++f) {
      AmanziMesh::Entity_ID_View cells;
      mesh->face_get_cells(f, AmanziMesh::Parallel_type::ALL, cells);
      matrices[f].reshape(cells.extent(0), cells.extent(0));
    }    
  }

  virtual void
  GetLocalDiagCopy(CompositeVector& X) const
  {
    auto Xv = X.ViewComponent("cell", true);

    AmanziMesh::Mesh const * mesh_ = mesh.get();
    Kokkos::parallel_for(
        matrices.size(),
        KOKKOS_LAMBDA(const int f) {
          AmanziMesh::Entity_ID_View cells;
          mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, cells);
          Kokkos::atomic_add(&Xv(cells(0), 0), matrices[f](0,0));
          if (cells.extent(0) > 1) {
            Kokkos::atomic_add(&Xv(cells(1),0), matrices[f](1,1));
          }
        });
  }

  
  virtual void ApplyMatrixFreeOp(const Operator* assembler,
          const CompositeVector& X, CompositeVector& Y) const {
    assembler->ApplyMatrixFreeOp(*this, X, Y);
  }

  //virtual void SymbolicAssembleMatrixOp(const Operator* assembler,
  //        const SuperMap& map, GraphFE& graph,
  //       int my_block_row, int my_block_col) const {
  //  assembler->SymbolicAssembleMatrixOp(*this, map, graph, my_block_row, my_block_col);
  //}

  //virtual void AssembleMatrixOp(const Operator* assembler,
  //        const SuperMap& map, MatrixFE& mat,
  //        int my_block_row, int my_block_col) const {
  //  assembler->AssembleMatrixOp(*this, map, mat, my_block_row, my_block_col);
  //}
  
  virtual void Rescale(const CompositeVector& scaling) {
    if (scaling.HasComponent("cell")) {
      const auto s_c = scaling.ViewComponent("cell",true);
      const Amanzi::AmanziMesh::Mesh* m = mesh.get(); 
      Kokkos::parallel_for(
          matrices.size(), 
          KOKKOS_LAMBDA(const int& f) {
            AmanziMesh::Entity_ID_View cells;
            m->face_get_cells(f, AmanziMesh::Parallel_type::ALL, cells);
            matrices[f](0,0) *= s_c(0, cells(0));
            if (cells.size() > 1) {
              matrices[f](0,1) *= s_c(0, cells(1));          
              matrices[f](1,0) *= s_c(0, cells(0));          
              matrices[f](1,1) *= s_c(0, cells(1));
            }          
          });
    }
  }
};

}  // namespace Operators
}  // namespace Amanzi


#endif

