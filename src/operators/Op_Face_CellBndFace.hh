/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@lanl.gov)
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
  Op_Face_CellBndFace(std::string& name,
               const Teuchos::RCP<const AmanziMesh::Mesh> mesh) :
      Op(OPERATOR_SCHEMA_BASE_FACE |
         OPERATOR_SCHEMA_DOFS_CELL | OPERATOR_SCHEMA_DOFS_FACE | OPERATOR_SCHEMA_DOFS_BNDFACE,
         name, mesh) {

    int nfaces_owned = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED);
    WhetStone::DenseMatrix null_matrix(2,2);
    matrices.resize(nfaces_owned, null_matrix);
    // for (int f=0; f!=nfaces_owned; ++f) {
    //   matrices[f].reshape(2,2);
    // }    
  }

  virtual void
  GetLocalDiagCopy(CompositeVector& X) const
  {
    AmanziMesh::Mesh const* mesh_ = mesh.get();
    auto Xc = X.ViewComponent("cell", true);
    Kokkos::parallel_for(
        matrices.size(),
        KOKKOS_LAMBDA(const int f) {
          AmanziMesh::Entity_ID_View cells;
          mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, cells);
          Xc(cells(0), 0) += matrices[f](0,0);
          if (cells.extent(0) > 1) {
            Kokkos::atomic_add(&Xc(cells(1),0), matrices[f](1,1));
          }
        });

    auto Xbf = X.ViewComponent("boundary_face", true);
    auto import_same_face = mesh->exterior_face_importer()->getNumSameIDs();
    Kokkos::parallel_for(
        import_same_face,
        KOKKOS_LAMBDA(const int bf) {
          // atomic not needed, each bf touched once
          Xbf(bf,0) += matrices[bf](1,1);
        });

    auto import_permute_face = mesh->exterior_face_importer()
        ->getPermuteFromLIDs_dv().view<AmanziDefaultDevice>();
    Kokkos::parallel_for(
        import_permute_face.extent(0),
        KOKKOS_LAMBDA(const int bf_offset) {
          // atomic not needed, each bf touched once
          Xbf(bf_offset+import_same_face,0) += matrices[import_permute_face(bf_offset)](1,1);
        });
  }


  virtual void ApplyMatrixFreeOp(const Operator* assembler,
          const CompositeVector& X, CompositeVector& Y) const {
    assembler->ApplyMatrixFreeOp(*this, X, Y);
  }

  // virtual void SymbolicAssembleMatrixOp(const Operator* assembler,
  //         const SuperMap& map, GraphFE& graph,
  //         int my_block_row, int my_block_col) const {
  //   assembler->SymbolicAssembleMatrixOp(*this, map, graph, my_block_row, my_block_col);
  // }

  // virtual void AssembleMatrixOp(const Operator* assembler,
  //         const SuperMap& map, MatrixFE& mat,
  //         int my_block_row, int my_block_col) const {
  //   assembler->AssembleMatrixOp(*this, map, mat, my_block_row, my_block_col);
  // }
  
  virtual void Rescale(const CompositeVector& scaling)
  {
    if ((scaling.HasComponent("cell")) &&
        (scaling.HasComponent("boundary_face"))) {

      { // context for views
        auto s_c = scaling.ViewComponent("cell", true);
        AmanziMesh::Mesh const * mesh_ = mesh.get();

        Kokkos::parallel_for(
            matrices.size(),
            KOKKOS_LAMBDA(const int f) {
              AmanziMesh::Entity_ID_View cells;
              mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, cells);
              matrices[f](0,0) *= s_c(cells(0),0);
              matrices[f](1,0) *= s_c(cells(0),0);
              if (cells.extent(0) > 1) {
                matrices[f](0,1) *= s_c(cells(1),0);
                matrices[f](1,1) *= s_c(cells(1),0);
              }
            });
      }

      {
        auto s_bnd = scaling.ViewComponent("boundary_face", true);

        // The importer gets us the face LID of a given boundary face LID
        // through two parts -- the first import_same_face are bf = f, then the
        // remainder are permutations.  Note this is guaranteed because we
        // know, by definition, that the faces and boundary faces are owned by
        // the same process.        
        auto import_same_face = mesh->exterior_face_importer()->getNumSameIDs();
        auto import_permute_face =
            mesh->exterior_face_importer()
            ->getPermuteFromLIDs_dv().view<AmanziDefaultDevice>();
      
        Kokkos::parallel_for(import_same_face,
                             KOKKOS_LAMBDA(const int bf) {
                               matrices[bf](0,1) *= s_bnd(bf,0);
                               matrices[bf](1,1) *= s_bnd(bf,0);
                             });

        Kokkos::parallel_for(import_permute_face.extent(0),
                             KOKKOS_LAMBDA(const int bf_offset) {
                               matrices[import_permute_face(bf_offset)](0,1) *=
                                   s_bnd(bf_offset+import_same_face,0);
                               matrices[import_permute_face(bf_offset)](1,1) *=
                                   s_bnd(bf_offset+import_same_face,0);
                             });
      }
    }
  }
};

}  // namespace Operators
}  // namespace Amanzi


#endif

