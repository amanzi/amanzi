/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/
/*!

  TreeOperators are the block analogue of Operators -- they provide
  a linear operator acting on a TreeVectorSpace.

  Note that these are really intended for enabling preconditioners -- it is
  unlikely that these need assembled for the operator itself, and therefore no
  ComputeResidual() methods are provided.

*/

#pragma once
#include <vector>

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Epetra_CrsMatrix.h"

#include "TreeVector.hh"
#include "TreeVectorSpace.hh"
#include "VerboseObject.hh"

namespace Amanzi {
namespace Operators {

// forward declarations
class SuperMap;
class MatrixFE;
class Operator;
namespace Impl {
class TreeOperator_BlockDiagonalPreconditioner;
}

class TreeOperator : public Matrix<TreeVector,TreeVectorSpace> {
public:
  using Vector_t = TreeVector;
  using VectorSpace_t = TreeVector::VectorSpace_t;

  TreeOperator();
  TreeOperator(Teuchos::ParameterList& plist);
  TreeOperator(const Teuchos::RCP<const TreeVectorSpace>& row_map,
               const Teuchos::RCP<const TreeVectorSpace>& col_map,
               Teuchos::ParameterList& plist);
  TreeOperator(const Teuchos::RCP<const TreeVectorSpace>& row_map,
               const Teuchos::RCP<const TreeVectorSpace>& col_map);
  TreeOperator(const Teuchos::RCP<const TreeVectorSpace>& tvs,
               Teuchos::ParameterList& plist);
  TreeOperator(const Teuchos::RCP<const TreeVectorSpace>& tvs);

  //
  // NOTE: this is not a shallow copy, but not quite a deep copy either!
  //
  // It is deep in the sense that it calls Clone on all leaf Operator objects,
  // but those Clone() methods themselves are shallow!
  TreeOperator(const TreeOperator& other);
  Teuchos::RCP<TreeOperator> Clone() const {
    return Teuchos::rcp(new TreeOperator(*this));
  }

  // accessors
  virtual const TreeVectorSpace& DomainMap() const override { return *col_map_; }
  virtual const TreeVectorSpace& RangeMap() const override { return *row_map_; }

  Teuchos::RCP<const TreeVectorSpace> get_domain_map() const { return col_map_; }
  Teuchos::RCP<const TreeVectorSpace> get_range_map() const { return row_map_; }
  Teuchos::RCP<const TreeVectorSpace> get_col_map() const { return col_map_; }
  Teuchos::RCP<const TreeVectorSpace> get_row_map() const { return row_map_; }
  Teuchos::RCP<const SuperMap> get_col_supermap() const { return col_supermap_; }
  Teuchos::RCP<const SuperMap> get_row_supermap() const { return row_supermap_; }

  // here for use by clients that still assume square
  Teuchos::RCP<const SuperMap> get_supermap() const { return get_row_supermap(); }

  Teuchos::RCP<TreeOperator> get_block(std::size_t i, std::size_t j) { return blocks_[i][j]; }
  Teuchos::RCP<const TreeOperator> get_block(std::size_t i, std::size_t j) const { return blocks_[i][j]; }
  Teuchos::RCP<Operator> get_operator() { return data_; }
  Teuchos::RCP<const Operator> get_operator() const { return data_; }
  Teuchos::RCP<Operator> get_operator_block(std::size_t i, std::size_t j) {
    if (get_block(i,j) != Teuchos::null) return get_block(i,j)->get_operator();
    return Teuchos::null;
  }
  Teuchos::RCP<const Operator> get_operator_block(std::size_t i, std::size_t j) const {
    if (get_block(i,j) != Teuchos::null) return get_block(i,j)->get_operator();
    return Teuchos::null;
  }

  Teuchos::RCP<Epetra_CrsMatrix> A() { return A_; }
  Teuchos::RCP<const Epetra_CrsMatrix> A() const { return A_; }

  void set_block(std::size_t i, std::size_t j, const Teuchos::RCP<TreeOperator>& op);
  void set_operator(const Teuchos::RCP<Operator>& op);
  void set_operator_block(std::size_t i, std::size_t j, const Teuchos::RCP<Operator>& op);

  bool IsSquare() const { return get_col_size() == get_row_size(); }
  std::size_t get_col_size() const { return col_size_; }
  std::size_t get_row_size() const { return row_size_; }

  void set_coloring(int num_colors, const Teuchos::RCP<std::vector<int>>& coloring) {
    num_colors_ = num_colors;
    coloring_ = coloring;
  }

  // i/o
  std::string PrintDiagnostics() const { return std::string(); }

  // forward operator
  virtual int Apply(const TreeVector& X, TreeVector& Y) const override {
    return Apply(X,Y,0.0);
  }
  int Apply(const TreeVector& X, TreeVector& Y, double scalar) const;
  int ApplyFlattened(const TreeVector& X, TreeVector& Y) const;
  int ApplyAssembled(const TreeVector& X, TreeVector& Y) const;

  void SymbolicAssembleMatrix();
  void AssembleMatrix();

  // inverse operator
  virtual int ApplyInverse(const TreeVector& X, TreeVector& Y) const override;

  void set_inverse_parameters(const std::string& prec_name,
                              const Teuchos::ParameterList& plist);

  virtual void set_inverse_parameters(Teuchos::ParameterList& plist) override;
  virtual void InitializeInverse() override;
  virtual void ComputeInverse() override;

  // Inverse diagnostics... these may change
  virtual double residual() const override {
    if (preconditioner_.get()) return preconditioner_->residual();
    return 0.;
  }
  virtual int num_itrs() const override {
    if (preconditioner_.get()) return preconditioner_->num_itrs();
    return 0;
  }
  virtual int returned_code() const override {
    if (preconditioner_.get()) return preconditioner_->returned_code();
    return 0;
  }
  virtual std::string returned_code_string() const override {
    if (preconditioner_.get()) return preconditioner_->returned_code_string();
    return "success";
  }
  virtual std::string name() const override {
    if (preconditioner_.get()) return std::string("TreeOperator: ")+preconditioner_->name();
    return "TreeOperator: block diagonal";
  }


 protected:
int ApplyInverseBlockDiagonal_(const TreeVector& X, TreeVector& Y) const;

 private:
  friend Impl::TreeOperator_BlockDiagonalPreconditioner;

  Teuchos::RCP<const TreeVectorSpace> row_map_, col_map_;
  std::size_t row_size_, col_size_;
  Teuchos::Array<Teuchos::Array<Teuchos::RCP<TreeOperator> > > blocks_;
  std::vector<std::vector<Teuchos::RCP<Operator>>> leaves_;
  Teuchos::RCP<Operator> data_;

  Teuchos::RCP<Epetra_CrsMatrix> A_;
  Teuchos::RCP<MatrixFE> Amat_;
  Teuchos::RCP<SuperMap> row_supermap_, col_supermap_;

  Teuchos::RCP<Matrix<TreeVector>> preconditioner_;
  bool block_diagonal_;

  int num_colors_;
  Teuchos::RCP<std::vector<int>> coloring_;
  Teuchos::ParameterList inv_plist_;
  bool inited_, updated_, computed_;
  Teuchos::RCP<VerboseObject> vo_;
};


namespace Impl {
//
// This class simply wraps TreeOperator with ApplyInverse() calling
// ApplyInverseBlockDiagonal() so that it can be used as a preconditioner.  To
// provide the common interface to the client, TreeOperator's ApplyInverse must
// call the linear solver's ApplyInverse(), which wants a preconditioner, which
// might here be the block diagonal case.  This would be circular without
// introducing a separate object.
//
class TreeOperator_BlockDiagonalPreconditioner {
 public:
  using Vector_t = TreeOperator::Vector_t;
  using VectorSpace_t = TreeOperator::VectorSpace_t;

  TreeOperator_BlockDiagonalPreconditioner(TreeOperator& op) :
      op_(op) {}

  int Apply(const TreeVector& X, TreeVector& Y) const {
    Exceptions::amanzi_throw("TreeOperator Preconditioner does not implement Apply()");
    return 1;
  }

  int ApplyInverse(const TreeVector& X, TreeVector& Y) const {
    return op_.ApplyInverseBlockDiagonal_(X,Y);
  }
  void InitializeInverse() {
    for (int n=0; n!=op_.get_col_size(); ++n) {
      op_.get_block(n,n)->InitializeInverse();
    }
  }
  void ComputeInverse() {
    for (int n=0; n!=op_.get_col_size(); ++n) {
      op_.get_block(n,n)->ComputeInverse();
    }
  }

 private:
  TreeOperator& op_;
};

std::pair<int,int>
collectTreeOperatorLeaves(TreeOperator& tm, std::vector<std::vector<Teuchos::RCP<Operator>>>& leaves, std::size_t i, std::size_t j);

} // namespace Impl

}  // namespace Operators
}  // namespace Amanzi



