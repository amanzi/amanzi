/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_OP_CELL_FACECELL_HH_
#define AMANZI_OP_CELL_FACECELL_HH_

#include <vector>
#include "DenseMatrix.hh"
#include "Operator.hh"
#include "Op.hh"

namespace Amanzi {
namespace Operators {

class Op_Cell_FaceCell : public Op {
 public:
  Op_Cell_FaceCell(const std::string& name,
                   const Teuchos::RCP<const AmanziMesh::Mesh> mesh) :
      Op(OPERATOR_SCHEMA_BASE_CELL |
         OPERATOR_SCHEMA_DOFS_FACE |
         OPERATOR_SCHEMA_DOFS_CELL, name, mesh) {
    std::cout<<"Matrix init"<<std::endl;
    int ncells_owned = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
    A = CSR_Matrix(ncells_owned); 

    for (int c=0; c!=ncells_owned; ++c) {
      AmanziMesh::Entity_ID_View faces;
      mesh->cell_get_faces(c, faces);      // This perform the prefix_sum
      int nfaces = faces.extent(0); 
      int loc[2] = {nfaces,nfaces}; 
      A.set_shape_host(c, loc);
    }
    A.prefix_sum();
  }

  virtual void
  SumLocalDiag(CompositeVector& X) const
  {
    std::cout<<"Op_Cell_FaceCell.hh::SumLocalDiag"<<std::endl;
    
    AmanziMesh::Mesh const * mesh_ = mesh.get();
    auto Xv = X.ViewComponent("face", true);
    Kokkos::parallel_for(
        "Op_Cell_FaceCell::GetLocalDiagCopy",
        A.size(),
        KOKKOS_LAMBDA(const int f) {
          // Extract matrix 
          WhetStone::DenseMatrix<DeviceOnlyMemorySpace> lm(
            A.at(f),A.size(f,0),A.size(f,1)); 

          AmanziMesh::Entity_ID_View faces;
          mesh_->cell_get_faces(f, faces);
          Kokkos::atomic_add(&Xv(faces(0), 0), lm(0,0));
          if(faces.extent(0) > 1) {
            Kokkos::atomic_add(&Xv(faces(1),0), lm(1,1));
          }
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
    const Amanzi::AmanziMesh::Mesh* mesh_ = mesh.get(); 
    std::cout<<"Op_Cell_FaceCell::Rescale"<<std::endl;
    if (scaling.HasComponent("face")) {
      std::cout<<"Op_Cell_FaceCell::Rescale(hascomponent(face))"<<std::endl;
      const auto s_c = scaling.ViewComponent("face",true);
      Kokkos::parallel_for(
        "Op_Cell_FaceCell::Rescale",
        A.size(), 
        KOKKOS_LAMBDA(const int& c){
          AmanziMesh::Entity_ID_View faces;
          mesh_->cell_get_faces(c, faces);
          WhetStone::DenseMatrix<DeviceOnlyMemorySpace> lm = getFromCSR<WhetStone::DenseMatrix>(A,c); 
          for (int n = 0; n != faces.size(); ++n) {
            for (int m = 0; m != faces.size(); ++m) {
              lm(n,m) *= s_c(0,faces[n]);
            }
            lm(n,faces.size()) *= s_c(0,faces[n]);          
          }
        }); 
    }
    if (scaling.HasComponent("cell")) {
      std::cout<<"Op_Cell_FaceCell::Rescale(hascomponent(cell))"<<std::endl;
      const auto s_c = scaling.ViewComponent("cell",true);
      Kokkos::parallel_for(
        "Op_Cell_FaceCell::Rescale",
        A.size(), 
        KOKKOS_LAMBDA(const int& c){
          WhetStone::DenseMatrix<DeviceOnlyMemorySpace> lm = getFromCSR<WhetStone::DenseMatrix>(A,c); 
          int nfaces = mesh_->cell_get_num_faces(c);
          for (int m = 0; m != nfaces; ++m) {
            lm(nfaces,m) *= s_c(0,c);
          }
          lm(nfaces,nfaces) *= s_c(0,c);
        });
    }
  }
};

}  // namespace Operators
}  // namespace Amanzi


#endif


