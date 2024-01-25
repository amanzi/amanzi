/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#ifndef AMANZI_OP_EDGE_EDGE_HH_
#define AMANZI_OP_EDGE_EDGE_HH_

#include <vector>
#include "Operator.hh"
#include "Op.hh"

namespace Amanzi {
namespace Operators {

class Op_Edge_Edge : public Op {
 public:
  Op_Edge_Edge(const std::string& name,
               const Teuchos::RCP<const AmanziMesh::Mesh> mesh) :
      Op(OPERATOR_SCHEMA_BASE_EDGE |
         OPERATOR_SCHEMA_DOFS_EDGE, name, mesh) {
    diag = Teuchos::rcp(new Epetra_MultiVector(mesh->edge_map(false), 1));
    diag_shadow = Teuchos::rcp(new Epetra_MultiVector(mesh->edge_map(false), 1));
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
    if (scaling.hasComponent("edge")) {
      const Epetra_MultiVector& s_e = *scaling.viewComponent("edge", false);
      for (int k = 0; k != s_e.NumVectors(); ++k) {
        for (int i = 0; i != s_e.MyLength(); ++i) {
          (*diag)[k][i] *= s_e[0][i];
        }
      }
    }
  }
};

}  // namespace Operators
}  // namespace Amanzi

#endif


