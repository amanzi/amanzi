/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon (coonet@ornl.gov)
*/

//! Container for local matrices.

#ifndef AMANZI_OP_HH_
#define AMANZI_OP_HH_

#include "Teuchos_RCP.hpp"
#include "Kokkos_StaticCrsGraph.hpp"

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

namespace AmanziMesh {
class Mesh;
}
class CompositeVector;

namespace Operators {

class SuperMap;
class GraphFE;
class MatrixFE;
class Operator;

class Op {
 public:
  Op(int schema, const std::string& schema_string_,
     const Teuchos::RCP<const AmanziMesh::Mesh> mesh)
    : schema_old_(schema),
      schema_row_(schema),
      schema_col_(schema),
      schema_string(schema_string_),
      mesh_(mesh){};

  Op(const Schema& schema_row, const Schema& schema_col,
     const Teuchos::RCP<const AmanziMesh::Mesh> mesh)
    : schema_row_(schema_row), schema_col_(schema_col), mesh_(mesh)
  {
    schema_string =
      schema_row.CreateUniqueName() + '+' + schema_col.CreateUniqueName();
  }

  Op(int schema, const std::string& schema_string_)
    : schema_old_(schema),
      schema_string(schema_string_),
      mesh_(Teuchos::null){};

  virtual ~Op() = default;

  // Clean the operator without destroying memory
  void Init()
  {
    Kokkos::MDRangePolicy<Kokkos::Rank<2>> policy({0,0},
            {data_.extent(0), data_.extent(1)});
    Kokkos::parallel_for(policy,
                         KOKKOS_LAMBDA(const int i, const int j) {
                           data_(i,j) = 0.0;
                         });
  }

  // Restore pristine value of the matrices, i.e. before BCs.
  virtual int CopyShadowToMaster()
  {
    // FIXME --etc
    for (int i = 0; i != matrices.size(); ++i) {
      if (matrices_shadow[i].NumRows() != 0) {
        matrices[i] = matrices_shadow[i];
      }
    }
    *diag = *diag_shadow;
    return 0;
  }

  // For backward compatibility... must go away
  virtual void RestoreCheckPoint()
  {
    // FIXME --etc
    for (int i = 0; i != matrices.size(); ++i) {
      if (matrices_shadow[i].NumRows() != 0) {
        matrices[i] = matrices_shadow[i];
      }
    }
  }

  // Matching rules for schemas.
  virtual bool Matches(int match_schema, int matching_rule)
  {
    if (matching_rule == OPERATOR_SCHEMA_RULE_EXACT) {
      if ((match_schema & schema_old_) == schema_old_) return true;
    } else if (matching_rule == OPERATOR_SCHEMA_RULE_SUBSET) {
      if (match_schema & schema_old_) return true;
    }
    return false;
  }

  // linear operator functionality.
  virtual void
  ApplyMatrixFreeOp(const Operator* assembler, const CompositeVector& X,
                    CompositeVector& Y) const = 0;

  virtual void ApplyTransposeMatrixFreeOp(const Operator* assembler,
                                          const CompositeVector& X,
                                          CompositeVector& Y) const = 0;

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
  virtual void Rescale(double scaling)
  {
    Kokkos::MDRangePolicy<Kokkos::Rank<2>> policy({0,0},
            {data_.extent(0), data_.extent(1)});
    Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const int i,const int j) {
        data_(i,j) *= scaling;
        });
  }

  // access
  const Schema& schema_row() const { return schema_row_; }

 public:
  int schema_old_;
  Schema schema_row_, schema_col_;
  std::string schema_string;

  // diagonal matrix
  Kokkos::View<double**> data_;
  Kokkos::Crs<double> sparse_shadow_;
  Kokkos::View<double**> dense_shadow_;

 protected:
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
};

} // namespace Operators
} // namespace Amanzi


#endif
