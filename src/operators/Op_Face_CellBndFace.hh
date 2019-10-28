/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon (coonet@ornl.gov)
*/

//! An Op for face-based, cell + boundary face entities

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
                      const Teuchos::RCP<const AmanziMesh::Mesh> mesh)
    : Op(OPERATOR_SCHEMA_BASE_FACE | OPERATOR_SCHEMA_DOFS_CELL |
           OPERATOR_SCHEMA_DOFS_FACE | OPERATOR_SCHEMA_DOFS_BNDFACE,
         name, mesh)
  {
    // 2 cells per face, or a cell plus a boundary face, squared
    data = Kokkos::View<double**>(name, mesh->face_map(false)->getNodeNumElements(), 4);
  }

  virtual void
  getLocalDiagCopy(CompositeVector& X) const
  {
    AmanziMesh::Mesh const * mesh_ = mesh.get();

    auto Xc = X.ViewComponent("cell", true);
    Kokkos::parallel_for(
        data.extent(0),
        KOKKOS_LAMBDA(const int f) {
          AmanziMesh::Entity_ID_View cells;
          mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, cells);
          Xc(cells(0), 0) += data(f,0);
          if (cells.extent(0) > 1) {
            Xc(cells(1),0) += data(f,3);
          }
        });

    auto Xbf = X.ViewComponent("boundary_face", true);
    auto import_same_face = mesh->exterior_face_importer()->getNumSameIDs();
    Kokkos::parallel_for(import_same_face,
                         KOKKOS_LAMBDA(const int bf) {
                           Xbf(bf,0) += data(bf,3);
                         });

    auto import_permute_face = mesh->exterior_face_importer()
        ->getPermuteFromLIDs_dv().view<AmanziDefaultDevice>();
    Kokkos::parallel_for(
        import_permute_face.extent(0),
        KOKKOS_LAMBDA(const int bf_offset) {
          Xbf(bf_offset+import_same_face,0) += data(import_permute_face(bf_offset),3);
        });
  }
  
  virtual void
  ApplyMatrixFreeOp(const Operator* assembler, const CompositeVector& X,
                    CompositeVector& Y) const
  {
    assembler->ApplyMatrixFreeOp(*this, X, Y);
  }

  virtual void
  ApplyTransposeMatrixFreeOp(const Operator* assembler,
                             const CompositeVector& X, CompositeVector& Y) const
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
    if ((scaling.HasComponent("cell")) &&
        (scaling.HasComponent("boundary_face"))) {

      {
        auto s_c = scaling.ViewComponent("cell", true);
        AmanziMesh::Mesh const * mesh_ = mesh.get();

        Kokkos::parallel_for(
            data.extent(0),
            KOKKOS_LAMBDA(const int f) {
              AmanziMesh::Entity_ID_View cells;
              mesh->face_get_cells(f, AmanziMesh::Parallel_type::ALL, cells);
              data(f,0) *= s_c(cells(0),0);
              data(f,1) *= s_c(cells(0),0);
              if (cells.extent(0) > 1) {
                data(f,2) *= s_c(cells(1),0);
                data(f,3) *= s_c(cells(1),0);
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
                               data(bf,2) *= s_bnd(bf,0);
                               data(bf,3) *= s_bnd(bf,0);
                             });

        Kokkos::parallel_for(import_permute_face.extent(0),
                             KOKKOS_LAMBDA(const int bf_offset) {
                               data(import_permute_face(bf_offset),2) *=
                                   s_bnd(bf_offset+import_same_face,0);
                               data(import_permute_face(bf_offset),3) *=
                                   s_bnd(bf_offset+import_same_face,0);
                             });
      }
    }
  }

 protected:
  int nfaces_owned;
};

} // namespace Operators
} // namespace Amanzi


#endif
