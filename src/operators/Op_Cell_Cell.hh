/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon (coonet@ornl.gov)
*/

//! An Op for diagonal Cell-based entries.

#ifndef AMANZI_OP_CELL_CELL_HH_
#define AMANZI_OP_CELL_CELL_HH_

#include <vector>
#include "Operator.hh"
#include "Op.hh"

namespace Amanzi {
namespace Operators {

class Op_Cell_Cell : public Op {
 public:
  Op_Cell_Cell(const std::string& name,
               const Teuchos::RCP<const AmanziMesh::Mesh> mesh)
    : Op(OPERATOR_SCHEMA_BASE_CELL | OPERATOR_SCHEMA_DOFS_CELL, name, mesh)
  {
    data = Kokkos::View<double**>(name, mesh->cell_map(false)->getNodeNumElements(), 1);
  }

  virtual void
  getLocalDiagCopy(CompositeVector& X) const
  {
    auto Xv = X.ViewComponent("cell", false);
    Kokkos::parallel_for(
        "Op_Cell_Cell::getLocalDiagCopy",
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
    if (scaling.HasComponent("cell")) {
      auto scaling_v = scaling.ViewComponent("cell", false);
      AMANZI_ASSERT(scaling_v.extent(0) == data.extent(0));
      AMANZI_ASSERT(scaling_v.extent(1) == 1);
      
      Kokkos::parallel_for(scaling_v.extent(0),
                           KOKKOS_LAMBDA(const int i) {
                             data(i,0) *= scaling_v(i,0);
                           });
    }
  }
};

} // namespace Operators
} // namespace Amanzi


#endif
