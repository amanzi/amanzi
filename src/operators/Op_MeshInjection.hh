/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@lanl.gov)
*/
//! A local matrix from one schema to another, each with their own mesh.

#pragma once

#include "Epetra_Map.h"

#include "Op.hh"
#include "Mesh.hh"
#include "OperatorDefs.hh"
#include "Schema.hh"

namespace Amanzi {
namespace Operators {

class Op_MeshInjection : public Op {
 public:
  Op_MeshInjection(const Schema& schema_row,
                   const Teuchos::RCP<const AmanziMesh::Mesh>& mesh_row,
                   const Schema& schema_col,
                   const Teuchos::RCP<const AmanziMesh::Mesh>& mesh_col,
                   const Teuchos::RCP<Epetra_Map>& injection_,
                   bool transpose_ = false)
    : Op(schema_row, schema_col, mesh_row),
      transpose(transpose_),
      injection(injection_),
      mesh_col_(mesh_col)
  {
    AMANZI_ASSERT(schema_row_.size() == 1);
    AMANZI_ASSERT(schema_col_.size() == 1);
    schema_string = schema_row_.CreateUniqueName() + '+' + schema_col_.CreateUniqueName();
    diag = Teuchos::rcp(new Epetra_MultiVector(*injection, 1));
    diag_shadow = Teuchos::rcp(new Epetra_MultiVector(*injection, 1));
  }

  // linear operator functionality.
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

  virtual void Rescale(const CompositeVector& scaling) override { AMANZI_ASSERT(false); }


  const AmanziMesh::Mesh& get_row_mesh() const { return *mesh_; }
  const AmanziMesh::Mesh& get_col_mesh() const { return *mesh_col_; }

 public:
  bool transpose;
  Teuchos::RCP<Epetra_Map> injection;

 protected:
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_col_;
};


} // namespace Operators
} // namespace Amanzi
