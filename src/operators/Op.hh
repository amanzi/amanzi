/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_OP_HH_
#define AMANZI_OP_HH_

#include "Epetra_MultiVector.h"
#include "Teuchos_RCP.hpp"

#include "DenseMatrix.hh"
#include "Mesh.hh"
#include "OperatorDefs.hh"
#include "Schema.hh"

/*
  Op classes are small structs that play two roles:

  1. They provide a class name to the schema, enabling visitor patterns.
  2. They are a container for local matrices.

  This Op class is a container for storing local matrices that spans
  the whole mesh. The dofs vary and defined by operator's schema.
*/

namespace Amanzi {

namespace AmanziMesh { class Mesh; }
class CompositeVector;

namespace Operators {

class SuperMap;
class GraphFE;
class MatrixFE;
class Operator;

class Op {
 public:
  Op(int schema, const std::string& schema_string_,
     const Teuchos::RCP<const AmanziMesh::Mesh> mesh) :
      schema_string(schema_string_),
      schema_old_(schema),
      schema_row_(schema),
      schema_col_(schema),
      mesh_(mesh)
  {};

  Op(const Schema& schema_row, const Schema& schema_col,
     const Teuchos::RCP<const AmanziMesh::Mesh> mesh) :
      schema_row_(schema_row),
      schema_col_(schema_col),
      mesh_(mesh) {
    schema_string = schema_row.CreateUniqueName() + '+' + schema_col.CreateUniqueName();
  }

  Op(int schema, const std::string& schema_string_) :
      schema_string(schema_string_),
      schema_old_(schema),
      mesh_(Teuchos::null)
  {};

  virtual ~Op() = default;

  // Clean the operator without destroying memory
  void Init() {
    if (diag != Teuchos::null) {
      diag->PutScalar(0.0);
      diag_shadow->PutScalar(0.0);
    }

    WhetStone::DenseMatrix null_mat;
    for (int i = 0; i < matrices.size(); ++i) {
      matrices[i] = 0.0;
      matrices_shadow[i] = null_mat;
    }
  }

  // Restore pristine value of the matrices, i.e. before BCs.
  virtual int CopyShadowToMaster() {
    for (int i = 0; i != matrices.size(); ++i) {
      if (matrices_shadow[i].NumRows() != 0) {
        matrices[i] = matrices_shadow[i];
      }
    }
    *diag = *diag_shadow;
    return 0;
  }

  // For backward compatibility... must go away
  virtual void RestoreCheckPoint() {
    for (int i = 0; i != matrices.size(); ++i) {
      if (matrices_shadow[i].NumRows() != 0) {
        matrices[i] = matrices_shadow[i];
      }
    }
    *diag = *diag_shadow;
  }

  // Matching rules for schemas.
  virtual bool Matches(int match_schema, int matching_rule) {
    if (matching_rule == OPERATOR_SCHEMA_RULE_EXACT) {
      if (match_schema == schema_old_) return true;
    } else if (matching_rule == OPERATOR_SCHEMA_RULE_SUPERSET) {
      if ((match_schema & schema_old_) == schema_old_) return true;
    } else if (matching_rule == OPERATOR_SCHEMA_RULE_SUBSET) {
      if ((match_schema & schema_old_) == match_schema) return true;
    }
    return false;
  }

  // linear operator functionality.
  virtual void ApplyMatrixFreeOp(const Operator* assembler,
          const CompositeVector& X, CompositeVector& Y) const = 0;

  virtual void SymbolicAssembleMatrixOp(const Operator* assembler,
          const SuperMap& map, GraphFE& graph,
          int my_block_row, int my_block_col) const = 0;

  virtual void AssembleMatrixOp(const Operator* assembler,
          const SuperMap& map, MatrixFE& mat,
          int my_block_row, int my_block_col) const = 0;

  // Mutators of local matrices.
  // -- rescale local matrices in the container using a CV
  virtual void Rescale(const CompositeVector& scaling) = 0;
  // -- rescale local matrices in the container using a double
  virtual void Rescale(double scaling);

  // access
  const Schema& schema_row() const { return schema_row_; }
  const Schema& schema_col() const { return schema_col_; }
  int schema_old() const { return schema_old_; }

  // quality control (more work is needed)
  // -- veify that containers are not empty
  void Verify() const; 

 public:
  std::string schema_string;

  // diagonal matrix
  Teuchos::RCP<Epetra_MultiVector> diag;
  Teuchos::RCP<Epetra_MultiVector> diag_shadow;

  // collection of local matrices
  std::vector<WhetStone::DenseMatrix> matrices;
  std::vector<WhetStone::DenseMatrix> matrices_shadow;

 protected:
  int schema_old_;
  Schema schema_row_, schema_col_;
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
};


/* ******************************************************************
* Default implementation
****************************************************************** */
inline
void Op::Rescale(double scaling)
{
  if (scaling != 1.0) {
    for (int i = 0; i != matrices.size(); ++i) {
      matrices[i] *= scaling;
    }
    if (diag.get()) diag->Scale(scaling);
  }
}


/* ******************************************************************
* Rudimentary verification of container quality.
****************************************************************** */
inline
void Op::Verify() const
{
 int nmatrices = matrices.size();
 for (int i = 0; i < nmatrices; ++i) {
   AMANZI_ASSERT(matrices[i].NumRows() > 0 && matrices[i].NumCols() > 0);
 }
}

}  // namespace Operators
}  // namespace Amanzi


#endif


