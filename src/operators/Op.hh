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

#include <Kokkos_Core.hpp>
#include <Kokkos_Vector.hpp>
#include "Teuchos_RCP.hpp"

#include "DenseMatrix.hh"
#include "Mesh.hh"
#include "OperatorDefs.hh"
#include "Schema.hh"

#include "AmanziTypes.hh"
#include "AmanziVector.hh"

#include "CSR.hh"
#include "DenseMatrix_Vector.hh"
#include "DenseVector_Vector.hh"
#include "Mesh.hh"

namespace Amanzi {

class GraphFE;
class MatrixFE;
class SuperMap;

namespace Operators {

class Operator;

class Op {
 public:
  Op(const Schema& schema_row_,
     const Schema& schema_col_,
     const Teuchos::RCP<const AmanziMesh::Mesh> mesh_)
      : schema_row(schema_row_),
        schema_col(schema_col_),
        schema_old(-1),
        mesh(mesh_),
        schema_string(schema_row_.CreateUniqueName()+'+'+schema_col_.CreateUniqueName())
  {
    AMANZI_ASSERT(schema_row.base() == schema_col.base()); 
  }

  Op(int schema_old_,
     const std::string& name,     
     const Teuchos::RCP<const AmanziMesh::Mesh> mesh_)
      : schema_old(schema_old_),
        schema_col(schema_old_),
        schema_row(schema_old_),
        mesh(mesh_),
        schema_string(name) {}     

  virtual ~Op() = default;
  
  // Clean the operator without destroying memory
  void Zero();

  KOKKOS_INLINE_FUNCTION
  void Zero(const int i) {
    auto lm = A[i]; 
    lm.putScalar(0.); 
    // See PDE_DiffusionFV::ApplyBCs for canonical usage example. --etc
    //matrices(i).putScalar(0.);
  }

  // Matching rules for schemas.
  virtual bool Matches(int match_schema, int matching_rule);

  // linear operator functionality.
  virtual void
  SumLocalDiag(CompositeVector& X) const = 0;

  virtual void
  ApplyMatrixFreeOp(const Operator* assembler, const CompositeVector& X,
                    CompositeVector& Y) const = 0;

  void PreallocateWorkVectors() const;
  
  //virtual void ApplyTransposeMatrixFreeOp(const Operator* assembler,
  //                                        const CompositeVector& X,
  //                                        CompositeVector& Y) const = 0;

  virtual void
  SymbolicAssembleMatrixOp(const Operator* assembler, const SuperMap& map,
                          GraphFE& graph, int my_block_row,
                          int my_block_col) const = 0;

  virtual void
  AssembleMatrixOp(const Operator* assembler, const SuperMap& map,
                  MatrixFE& mat, int my_block_row, int my_block_col) const = 0;

  // Mutators of local matrices.
  // -- rescale local matrices in the container using a CV
  virtual void Rescale(const CompositeVector& scaling) = 0;

  // -- rescale local matrices in the container using a double
  virtual void Rescale(double scaling);

  // access
 public:
  int schema_old;
  Schema schema_row, schema_col;
  std::string schema_string;

  // diagonal matrix
  MultiVector_ptr_type diag;

  // collection of local matrices
  mutable DenseMatrix_Vector A;

  // work space for local vectors, domain and range
  mutable DenseVector_Vector v; // work vector in domain
  mutable DenseVector_Vector Av; // work vector in range

  Teuchos::RCP<const AmanziMesh::Mesh> mesh;
};

} // namespace Operators
} // namespace Amanzi


#endif
