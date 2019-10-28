/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon (coonet@ornl.gov)
*/

//! An Op for face-based, cell entities.

#ifndef AMANZI_OP_FACE_CELL_HH_
#define AMANZI_OP_FACE_CELL_HH_

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

class Op_Face_Cell : public Op {
 public:
  Op_Face_Cell(const std::string& name,
               const Teuchos::RCP<const AmanziMesh::Mesh> mesh)
    : Op(OPERATOR_SCHEMA_BASE_FACE | OPERATOR_SCHEMA_DOFS_CELL, name, mesh)
  {
    // 2 cells per face, squared
    data = Kokkos::View<double**>(name, mesh->face_map(false)->getNodeNumElements(), 4);
  }

  virtual void
  getLocalDiagCopy(CompositeVector& X) const
  {
    auto Xv = X.ViewComponent("cell", true);

    AmanziMesh::Mesh const * mesh_ = mesh.get();
    Kokkos::parallel_for(
        data.extent(0),
        KOKKOS_LAMBDA(const int f) {
          AmanziMesh::Entity_ID_View cells;
          mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, cells);
          Xv(cells(0), 0) += data(f,0);
          if (cells.extent(0) > 1) {
            Xv(cells(1),0) += data(f,3);
          }
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
    if (scaling.HasComponent("cell")) {
      auto scaling_v = scaling.ViewComponent("cell", true);

      AmanziMesh::Mesh const * mesh_ = mesh.get();
      Kokkos::parallel_for(
          data.extent(0),
          KOKKOS_LAMBDA(const int f) {
            AmanziMesh::Entity_ID_View cells;
            mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, cells);
            data(f, 0) *= scaling_v(cells(0),0);

            if (cells.extent(0) > 1) {
              data(f, 1) *= scaling_v(cells(0),0);
              data(f, 2) *= scaling_v(cells(1),0);
              data(f, 3) *= scaling_v(cells(1),0);
            }
          });
    }
  }

 protected:
  int nfaces_owned;
};

} // namespace Operators
} // namespace Amanzi


#endif
