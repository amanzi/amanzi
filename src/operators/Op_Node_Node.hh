/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_OP_NODE_NODE_HH_
#define AMANZI_OP_NODE_NODE_HH_

#include <vector>
#include "Operator.hh"
#include "Op.hh"

namespace Amanzi {
namespace Operators {

class Op_Node_Node : public Op {
 public:
  Op_Node_Node(std::string& name,
               const Teuchos::RCP<const AmanziMesh::Mesh> mesh) :
      Op(OPERATOR_SCHEMA_BASE_NODE |
         OPERATOR_SCHEMA_DOFS_NODE, name, mesh) {
    vals.resize(mesh->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED), 0.);
    vals_shadow = vals;
  }

  virtual void ApplyMatrixFreeOp(const Operator* assembler,
          const CompositeVector& X, CompositeVector& Y) const {
    assembler->ApplyMatrixFreeOp(*this, X, Y);
  }

  virtual void SymbolicAssembleMatrixOp(const Operator* assembler,
          const SuperMap& map, GraphFE& graph,
          int my_block_row, int my_block_col) const {
    assembler->SymbolicAssembleMatrixOp(*this,
            map, graph, my_block_row, my_block_col);
  }

  virtual void AssembleMatrixOp(const Operator* assembler,
          const SuperMap& map, MatrixFE& mat,
          int my_block_row, int my_block_col) const {
    assembler->AssembleMatrixOp(*this, map, mat,
            my_block_row, my_block_col);
  }
  
  virtual void Rescale(const CompositeVector& scaling) {
    if (scaling.HasComponent("node")) {
      const Epetra_MultiVector& s_c = *scaling.ViewComponent("node",false);
      for (int i = 0; i != vals.size(); ++i) {
        vals[i] *= s_c[0][i];
      }
    }
  }
};

}  // namespace Operators
}  // namespace Amanzi


#endif


