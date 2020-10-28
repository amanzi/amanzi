/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Konstantin Lipnikov (lipnikov@lanl.gov)
*/

//! <MISSING_ONELINE_DOCSTRING>

#ifndef AMANZI_TREE_OPERATOR_HH_
#define AMANZI_TREE_OPERATOR_HH_

#include <vector>

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "TreeVector.hh"
#include "TreeVectorSpace.hh"
#include "VerboseObject.hh"

#include "Operator.hh"

/* ******************************************************************
  TreeOperators are the block analogue of Operators -- they provide
  a linear operator acting on a TreeVectorSpace.  They are currently
  assumed R^n -> R^n, and furthermore each block is currently assumed
  to be from R^m --> R^m for n = i*m where i is an integer (every
  block's space is the same).

  Note that these are really intended for preconditioners -- it is
  unlikely that these need assembled for the operator itself, and
  therefore no ComputeResidual() methods are provided.  It would be
  difficult to manage a RHS for these systems.

  Future work will relax this constraint, but currently this can be
  used for things like multi-phased flow, thermal Richards, etc.
****************************************************************** */

namespace Amanzi {
namespace Operators {

class SuperMap;
class MatrixFE;

class TreeOperator : public Matrix<TreeVector,TreeVectorSpace> {
 public:
  using Vector_t = TreeVector;
  using VectorSpace_t = TreeVectorSpace;

  TreeOperator() : block_diagonal_(false){};
  TreeOperator(Teuchos::RCP<const TreeVectorSpace> tvs);

  Teuchos::RCP<TreeOperator> get_block(std::size_t i, std::size_t j) { return blocks_[i][j]; }
  Teuchos::RCP<const TreeOperator> get_block(std::size_t i, std::size_t j) const { return blocks_[i][j]; }
  Teuchos::RCP<Operator> get_operator() { return data_; }
  Teuchos::RCP<const Operator> get_operator() const { return data_; }
  Teuchos::RCP<Operator> get_operator_block(std::size_t i, std::size_t j) {
    if (get_block(i,j) != Teuchos::null) return get_block(i,j)->get_operator();
    return Teuchos::null;
  }

  // main members
  void SetOperatorBlock(int i, int j, const Teuchos::RCP<const Operator>& op);

  void getLocalDiagCopy(TreeVector& tv) const;

  int apply(const TreeVector& X, TreeVector& Y) const;
  int applyAssembled(const TreeVector& X, TreeVector& Y) const;
  int applyInverse(const TreeVector& X, TreeVector& Y) const;

  void SymbolicAssembleMatrix();
  void AssembleMatrix();

  const Teuchos::RCP<const TreeVectorSpace> getDomainMap() const { return tvs_; }
  const Teuchos::RCP<const TreeVectorSpace> getRangeMap() const { return tvs_; }
  const Teuchos::RCP<const TreeVectorSpace> getRowMap() const { return tvs_; }  

  void set_inverse_parameters(const std::string& prec_name,
                              const Teuchos::ParameterList& plist);
  void set_inverse_parameters(Teuchos::ParameterList& plist) override final;
  void initializeInverse() override final;
  void computeInverse() override final;

  bool IsSquare() const { return get_col_size() == get_row_size(); }
  std::size_t get_col_size() const { return col_size_; }
  std::size_t get_row_size() const { return row_size_; }

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

  // access
  Teuchos::RCP<Matrix_type> A() { return A_; }
  Teuchos::RCP<const Matrix_type> A() const { return A_; }
  Teuchos::RCP<const SuperMap> getSuperMap() const { return smap_; }

 private:
  Teuchos::RCP<const TreeVectorSpace> tvs_;

  Teuchos::RCP<const TreeVectorSpace> row_map_, col_map_;
  std::size_t row_size_, col_size_;
  Teuchos::Array<Teuchos::Array<Teuchos::RCP<TreeOperator> > > blocks_;
  std::vector<std::vector<Teuchos::RCP<Operator>>> leaves_;
  Teuchos::RCP<Operator> data_;

  Teuchos::RCP<Matrix_type> A_;
  Teuchos::RCP<MatrixFE> Amat_;
  Teuchos::RCP<SuperMap> smap_;

  Teuchos::RCP<Matrix<TreeVector,TreeVectorSpace>> preconditioner_;
  bool block_diagonal_;

  int num_colors_;
  Teuchos::RCP<std::vector<int>> coloring_;
  Teuchos::ParameterList inv_plist_;
  bool inited_, updated_, computed_;
  Teuchos::RCP<VerboseObject> vo_;
};

} // namespace Operators
} // namespace Amanzi


#endif
