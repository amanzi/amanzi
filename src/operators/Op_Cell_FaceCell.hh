/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon (coonet@ornl.gov)
*/

//! <MISSING_ONELINE_DOCSTRING>

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
                   const Teuchos::RCP<const AmanziMesh::Mesh> mesh)
    : Op(OPERATOR_SCHEMA_BASE_CELL | OPERATOR_SCHEMA_DOFS_FACE |
           OPERATOR_SCHEMA_DOFS_CELL, name, mesh)
  {
    //WhetStone::DenseMatrix null_matrix;
    //matrices.resize(
    //  mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED),
    //  null_matrix);
    //matrices_shadow = matrices;
    data = Kokkos::View<double**>(name, 
      mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED), 1);
    shadow = Kokkos::View<double**>(name, 
      mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED), 1);
  }

  virtual void
  getLocalDiagCopy(CompositeVector& X) const
  {
    auto Xv = X.ViewComponent("cell", false);
    Kokkos::parallel_for(
        "Op_Cell_FaceCell::getLocalDiagCopy",
        Xv.extent(0),
        KOKKOS_LAMBDA(const int i) {
          Xv(i,0) += data(i,0);
        });
  }

  virtual void
  ApplyMatrixFreeOp(const Operator* assembler, 
                    const CompositeVector& X,
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

    const Amanzi::AmanziMesh::Mesh* m = mesh.get(); 
    if (scaling.HasComponent("face")) {
      auto s_c = scaling.ViewComponent("face", true);
      Kokkos::parallel_for(data.extent(0),
        KOKKOS_LAMBDA(const int c){
          Kokkos::View<AmanziMesh::Entity_ID*> faces;
          m->cell_get_faces(c, faces);
            for (int n = 0; n != faces.extent(0); ++n) {
              for (int m = 0; m != faces.extent(0); ++m) {
                data(n, m) *= s_c(0,faces[n]);
              }
            data(n, faces.extent(0)) *= s_c(0,faces[n]);
          }
        }); 
    }

    if (scaling.HasComponent("cell")) {
      auto s_c = scaling.ViewComponent("cell", true);
      AMANZI_ASSERT(s_c.extent(0) == data.extent(0));

      AmanziMesh::Entity_ID_List face;
      Kokkos::parallel_for(data.extent(0),
        KOKKOS_LAMBDA(const int c){
          // TODO ! not callable on device 
          int nfaces = m->cell_get_num_faces(c);
          for (int m = 0; m != nfaces; ++m) {
            data(nfaces, m) *= s_c(0,c);
          }
          data(nfaces, nfaces) *= s_c(0,c);
        }); 
    }
  }
};

} // namespace Operators
} // namespace Amanzi


#endif
