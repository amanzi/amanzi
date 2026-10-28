/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

//! A local matrix from one schema to another, each with their own mesh.
/*
  Operators

  A genuinely rectangular, cross-mesh local matrix block: row entity space is
  one mesh, column entity space is another, with exactly one nonzero per row
  (a scalar coefficient, `diag`), mapped row-to-column via an injection --
  for row-local index i, injection->getGlobalElement(i) gives the
  corresponding GLOBAL index of the matching column-side entity.

  Used by Operator_Schema as the off-diagonal coupling block for a 2x2
  TreeOperator coupling two PKs defined on different meshes (e.g. a
  surface/subsurface water MPC's rectangular Jacobian cross-terms), where
  each surface cell maps to exactly one subsurface top cell (or vice versa).
*/

#ifndef AMANZI_OP_MESH_INJECTION_HH_
#define AMANZI_OP_MESH_INJECTION_HH_

#include "AmanziTypes.hh"
#include "Op.hh"
#include "Operator.hh"
#include "Mesh.hh"
#include "OperatorDefs.hh"
#include "Schema.hh"

namespace Amanzi {
namespace Operators {

class Op_MeshInjection : public Op {
 public:
  // injection: built from and indexed by the "short" side's local CELL
  // index -- i.e. when transpose=false, that's mesh_row (row-local index
  // i maps via injection->getGlobalElement(i) to a GID on mesh_col's CELL
  // map); when transpose=true, the SAME kind of map is instead indexed by
  // mesh_col's local CELL index, mapping to a GID on mesh_row's CELL map
  // (so a single injection map, with transpose flipped, drives both
  // directions of a coupling -- e.g. the (0,1) and (1,0) blocks of a 2x2
  // TreeOperator -- without needing to build two separate maps). diag is
  // always sized/indexed to match injection's own domain (the short side),
  // regardless of transpose.
  Op_MeshInjection(const Schema& schema_row,
                   const Teuchos::RCP<const AmanziMesh::Mesh>& mesh_row,
                   const Schema& schema_col,
                   const Teuchos::RCP<const AmanziMesh::Mesh>& mesh_col,
                   const Map_ptr_type& injection,
                   bool transpose = false)
    : Op(schema_row, schema_col, mesh_row), transpose(transpose), injection(injection),
      mesh_col_(mesh_col)
  {
    AMANZI_ASSERT(schema_row.size() == 1);
    AMANZI_ASSERT(schema_col.size() == 1);
    schema_string = schema_row.CreateUniqueName() + '+' + schema_col.CreateUniqueName();
    diag = Teuchos::rcp(new MultiVector_type(injection, 1));
  }

  // linear operator functionality.
  virtual void SumLocalDiag(CompositeVector& X) const override
  {
    // Off-diagonal blocks are never separately diagonal-summed (they are
    // not square); this is only ever called on a TreeOperator's diagonal
    // blocks.
    AMANZI_ASSERT(false);
  }

  virtual void
  ApplyMatrixFreeOp(const Operator* assembler, const CompositeVector& X, CompositeVector& Y) const override
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

  const AmanziMesh::Mesh& get_row_mesh() const { return *mesh; }
  const AmanziMesh::Mesh& get_col_mesh() const { return *mesh_col_; }

 public:
  bool transpose;
  Map_ptr_type injection;
  Teuchos::RCP<MultiVector_type> diag;

 protected:
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_col_;
};


} // namespace Operators
} // namespace Amanzi

#endif
