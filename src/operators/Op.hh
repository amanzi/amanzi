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

#include "Teuchos_RCP.hpp"

#include "DenseMatrix.hh"
#include "BCs.hh"
#include "Mesh.hh"
#include "OperatorDefs.hh"

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
  Op(int schema_, std::string& schema_string_,
     const Teuchos::RCP<const AmanziMesh::Mesh> mesh) :
      schema(schema_),
      schema_string(schema_string_),
      mesh_(mesh)
  {};

  // Clean the operator without destroying memory
  void Init() {
    if (vals.size() > 0) {
      for (int i = 0; i != vals.size(); ++i) {
        vals[i] = 0.0;
        vals_shadow[i] = 0.0;
      }
    }

    WhetStone::DenseMatrix null_mat;
    if (matrices.size() > 0) {
      for (int i = 0; i != matrices.size(); ++i) {
        matrices[i] = 0.0;
        matrices_shadow[i] = null_mat;
      }
    }
  }
    
  // Restore pristine value of the matrices, i.e. before BCs.
  virtual int CopyShadowToMaster() {
    for (int i = 0; i != matrices.size(); ++i) {
      if (matrices_shadow[i].NumRows() != 0) {
        matrices[i] = matrices_shadow[i];
      }
    }
    vals = vals_shadow;
    return 0;
  }

  // For backward compatibility... must go away
  virtual void RestoreCheckPoint() {
    for (int i = 0; i != matrices.size(); ++i) {
      if (matrices_shadow[i].NumRows() != 0) {
        matrices[i] = matrices_shadow[i];
      }
    }
    vals = vals_shadow;
  }

  // Matching rules for schemas.
  virtual bool Matches(int match_schema, int matching_rule) {
    if (matching_rule == OPERATOR_SCHEMA_RULE_EXACT) {
      if ((match_schema & schema) == schema) return true;
    } else if (matching_rule == OPERATOR_SCHEMA_RULE_SUBSET) {
      if (match_schema & schema) return true;
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
  virtual void Rescale(double scaling) {
    for (int i = 0; i != matrices.size(); ++i) {
      matrices[i] *= scaling;
    }
  }

 public:
  int schema;
  std::string schema_string;

  std::vector<double> vals;
  std::vector<double> vals_shadow;  
  std::vector<WhetStone::DenseMatrix> matrices;
  std::vector<WhetStone::DenseMatrix> matrices_shadow;

 protected:
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
};

}  // namespace Operators
}  // namespace Amanzi


#endif


