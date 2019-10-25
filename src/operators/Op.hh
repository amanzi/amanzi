/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon (coonet@ornl.gov)
*/

//! Container for local matrices.

/*

  Op classes are a small container that generically stores local matrices and
  provides metadata about the structure of those local matrices via visitor
  pattern.
 
*/

#ifndef AMANZI_OP_HH_
#define AMANZI_OP_HH_

#include "Teuchos_RCP.hpp"
#include "Kokkos_View.hpp"

#include "OperatorDefs.hh"
#include "Schema.hh"


namespace Amanzi {

namespace AmanziMesh {
class Mesh;
}

namespace Operators {

class SuperMap;
class GraphFE;
class MatrixFE;
class Operator;

struct Op {
 public:
  Op(const Schema& schema_row_,
     const Schema& schema_col_,
     const Teuchos::RCP<const AmanziMesh::Mesh> mesh_)
      : schema_row(schema_row_),
        schema_col(schema_col_),
        mesh(mesh_),
        schema_string(schema_row.CreateUniqueName()+'+'+schema_col.CreateUniqueName())
  {
    AMANZI_ASSERT(schema_row_.base() == schema_col.base());
  }

  // Clean the operator without destroying memory
  void Zero();

  // Restore pristine value of the matrices, i.e. before BCs.
  virtual int CopyShadowToMaster();

  // Matching rules for schemas.
  virtual bool Matches(int match_schema, int matching_rule);

  // linear operator functionality.
  virtual void
  ApplyMatrixFreeOp(const Operator* assembler,
                    const CompositeVector& X,
                    CompositeVector& Y) const = 0;

  virtual void ApplyTransposeMatrixFreeOp(const Operator* assembler,
                                          const CompositeVector& X,
                                          CompositeVector& Y) const = 0;

  virtual void
  SymbolicAssembleMatrixOp(const Operator* assembler,
                           const SuperMap& map,
                           GraphFE& graph,
                           int my_block_row,
                           int my_block_col) const = 0;

  virtual void
  AssembleMatrixOp(const Operator* assembler,
                   const SuperMap& map,
                   MatrixFE& mat,
                   int my_block_row,
                   int my_block_col) const = 0;

  // Mutators of local matrices.
  // -- rescale local matrices in the container using a CV
  virtual void Rescale(const CompositeVector& scaling) = 0;

  // -- rescale local matrices in the container using a double
  virtual void Rescale(double scaling);

  // access
  const Schema& schema_row() const { return schema_row_; }

 public:
  // diagonal matrix
  std::string schema_string;

  Kokkos::View<double**> data;
  Kokkos::View<double**> shadow;
  Kokkos::View<int*> shadow_indices;

 protected:
  int schema_old_;
  Schema schema_row_, schema_col_;

  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
};

} // namespace Operators
} // namespace Amanzi


#endif
