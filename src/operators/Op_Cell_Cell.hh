/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_OP_CELL_CELL_HH_
#define AMANZI_OP_CELL_CELL_HH_

#include <vector>
#include "Operator.hh"
#include "Op.hh"

namespace Amanzi {
namespace Operators {

class Op_Cell_Cell : public Op {
 public:
  Op_Cell_Cell(const std::string& name, const Teuchos::RCP<const AmanziMesh::Mesh> mesh)
    : Op(OPERATOR_SCHEMA_BASE_CELL | OPERATOR_SCHEMA_DOFS_CELL, name, mesh)
  {
    diag =
      Teuchos::rcp(new MultiVector_type(mesh->getMap(AmanziMesh::Entity_kind::CELL, false), 1));
  }

  virtual void SumLocalDiag(CompositeVector& X) const
  {
    X.getComponent("cell", false)->update(1., *diag, 1.);
  }

  virtual void
  ApplyMatrixFreeOp(const Operator* assembler, const CompositeVector& X, CompositeVector& Y) const
  {
    assembler->ApplyMatrixFreeOp(*this, X, Y);
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
    if (scaling.hasComponent("cell")) {
      auto scaling_v = scaling.viewComponent("cell", false);
      auto diag_v = diag->getLocalView<DefaultDevice>(Tpetra::Access::ReadWrite);
      AMANZI_ASSERT(scaling_v.extent(0) == diag_v.extent(0));
      AMANZI_ASSERT(scaling_v.extent(1) == diag_v.extent(1));
      Kokkos::MDRangePolicy<Kokkos::Rank<2>> policy({ 0, 0 },
                                                    { diag_v.extent(0), diag_v.extent(1) });
      Kokkos::parallel_for(
        "Op_Cell_Cell::Rescale", policy, KOKKOS_LAMBDA(const int& i, const int& j) {
          diag_v(i, j) *= scaling_v(i, j);
        });
    }
  }
};

} // namespace Operators
} // namespace Amanzi


#endif
