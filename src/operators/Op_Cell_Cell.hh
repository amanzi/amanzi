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
  Op_Cell_Cell(const std::string& name,
               const Teuchos::RCP<const AmanziMesh::Mesh> mesh) :
      Op(OPERATOR_SCHEMA_BASE_CELL | OPERATOR_SCHEMA_DOFS_CELL,
         name, mesh) {

    diag = Teuchos::rcp(new MultiVector_type(mesh->cell_map(false), 1));
    diag_shadow = Teuchos::rcp(new MultiVector_type(mesh->cell_map(false), 1));
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
  
  virtual void Rescale(const CompositeVector& scaling) {
    assert(false);
    #if 0 
    if (scaling.HasComponent("cell")) {
      const Epetra_MultiVector& s_c = *scaling.ViewComponent("cell", false);
      AMANZI_ASSERT(s_c.MyLength() == diag->MyLength());
      for (int k = 0; k != s_c.NumVectors(); ++k) {
        for (int i = 0; i != s_c.MyLength(); ++i) {
          (*diag)[k][i] *= s_c[0][i];
        }
      }
    }
    #endif 
  }
};

}  // namespace Operators
}  // namespace Amanzi


#endif

