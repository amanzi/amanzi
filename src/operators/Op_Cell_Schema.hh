/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@lanl.gov)

  This operator is a container for local rectangular matrices of 
  length ncells and with dofs given by two schemas.
*/

#ifndef AMANZI_OP_CELL_SCHEMA_HH_
#define AMANZI_OP_CELL_SCHEMA_HH_

#include <vector>
#include "DenseMatrix.hh"
#include "Operator.hh"
#include "Op.hh"
#include "Schema.hh"

namespace Amanzi {
namespace Operators {

class Op_Cell_Schema : public Op {
 public:
  Op_Cell_Schema(const Schema& schema_row,
                 const Schema& schema_col,
                 const Teuchos::RCP<const AmanziMesh::Mesh> mesh)
    : Op(schema_row, schema_col, mesh)
  {
    WhetStone::DenseMatrix null_matrix;
    matrices.resize(mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED),
                    null_matrix);
    matrices_shadow = matrices;
  }
  ~Op_Cell_Schema(){};

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

  // incomplete members
  virtual void Rescale(const CompositeVector& scaling) { AMANZI_ASSERT(0); }
};

} // namespace Operators
} // namespace Amanzi

#endif
