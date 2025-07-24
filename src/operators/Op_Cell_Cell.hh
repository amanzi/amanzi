/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  Operators

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
  Op_Cell_Cell(const std::string& name, const Teuchos::RCP<const AmanziMesh::Mesh> mesh, int nvec)
    : Op(OPERATOR_SCHEMA_BASE_CELL | OPERATOR_SCHEMA_DOFS_CELL, name, mesh)
  {
    diag = Teuchos::rcp(
      new Epetra_MultiVector(mesh->getMap(AmanziMesh::Entity_kind::CELL, false), nvec));
    diag_shadow = Teuchos::rcp(
      new Epetra_MultiVector(mesh->getMap(AmanziMesh::Entity_kind::CELL, false), nvec));
  }

  virtual Teuchos::RCP<Op> DeepClone() const override
  {
    auto op = Teuchos::rcp(new Op_Cell_Cell(*this));
    *op->diag = *diag;
    *op->diag_shadow = *diag_shadow;
    return op;
  }

  virtual void ApplyMatrixFreeOp(const Operator* assembler,
                                 const CompositeVector& X,
                                 CompositeVector& Y) const override
  {
    assembler->ApplyMatrixFreeOp(*this, X, Y);
  }

  virtual void SymbolicAssembleMatrixOp(const Operator* assembler,
                                        const SuperMap& map,
                                        GraphFE& graph,
                                        int my_block_row,
                                        int my_block_col) const override
  {
    assembler->SymbolicAssembleMatrixOp(*this, map, graph, my_block_row, my_block_col);
  }

  virtual void AssembleMatrixOp(const Operator* assembler,
                                const SuperMap& map,
                                MatrixFE& mat,
                                int my_block_row,
                                int my_block_col) const override
  {
    assembler->AssembleMatrixOp(*this, map, mat, my_block_row, my_block_col);
  }

  virtual void Rescale(const CompositeVector& scaling) override
  {
    if (scaling.HasComponent("cell")) {
      const Epetra_MultiVector& s_c = *scaling.ViewComponent("cell", false);
      AMANZI_ASSERT(s_c.MyLength() == diag->MyLength());
      for (int k = 0; k != s_c.NumVectors(); ++k) {
        for (int i = 0; i != s_c.MyLength(); ++i) {
          (*diag)[k][i] *= s_c[0][i];
        }
      }
    }
  }
};

} // namespace Operators
} // namespace Amanzi


#endif
