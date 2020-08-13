/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#ifndef AMANZI_TREE_OPERATOR_HH_
#define AMANZI_TREE_OPERATOR_HH_

#include <vector>

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "TreeVector.hh"
#include "TreeVectorSpace.hh"
#include "VerboseObject.hh"

#include "InverseFactory.hh"
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

namespace Impl {
class TreeOperator_BlockPreconditioner;
}

class TreeOperator {
 public:
  using Vector_t = TreeVector;
  using VectorSpace_t = TreeVector::VectorSpace_t;

  TreeOperator() : block_diagonal_(false) {};
  TreeOperator(Teuchos::RCP<const TreeVectorSpace> tvs);
  // TreeOperator(Teuchos::RCP<const TreeVectorSpace> tvs, int nblocks);
  virtual ~TreeOperator() = default;

  // main members
  void SetOperatorBlock(int i, int j, const Teuchos::RCP<Operator>& op);
  
  virtual int Apply(const TreeVector& X, TreeVector& Y) const;
  virtual int ApplyAssembled(const TreeVector& X, TreeVector& Y) const;
  virtual int ApplyInverse(const TreeVector& X, TreeVector& Y) const;
  int ApplyFlattened(const TreeVector& X, TreeVector& Y) const;

  const TreeVectorSpace& DomainMap() const { return *tvs_; }
  const TreeVectorSpace& RangeMap() const { return *tvs_; }
  Teuchos::RCP<SuperMap> getSuperMap() const { return smap_; }
  
  void SymbolicAssembleMatrix();
  void AssembleMatrix();

  void InitializeInverse(const std::string& prec_name,
                         const Teuchos::ParameterList& plist);
  void InitializeInverse(Teuchos::ParameterList& plist) {
    inv_plist_ = plist;
    inited_ = true;
  }
  void UpdateInverse();
  void ComputeInverse();

  // Inverse diagnostics... these may change
  double residual() const {
    AMANZI_ASSERT(preconditioner_.get());
    return preconditioner_->residual();
  }
  int num_itrs() const {
    AMANZI_ASSERT(preconditioner_.get());
    return preconditioner_->num_itrs();
  }
  int returned_code() const { 
    AMANZI_ASSERT(preconditioner_.get());
    return preconditioner_->returned_code();
  }
  
  // access
  Teuchos::RCP<Epetra_CrsMatrix> A() { return A_; } 
  Teuchos::RCP<const Epetra_CrsMatrix> A() const { return A_; } 
  Teuchos::RCP<SuperMap> smap() const { return smap_; }

  Teuchos::RCP<const Operator> GetOperatorBlock(int i, int j) const { return blocks_[i][j];}
  int GetNumberBlocks() const {return blocks_.size();}
  void set_coloring(int num_colors, const Teuchos::RCP<std::vector<int>>& coloring) {
    num_colors_ = num_colors;
    coloring_ = coloring;
  }

  // i/o
  std::string PrintDiagnostics() const;

 protected:
  int ApplyInverseBlockDiagonal_(const TreeVector& X, TreeVector& Y) const;
  

 private:
  friend Impl::TreeOperator_BlockPreconditioner;

  Teuchos::RCP<const TreeVectorSpace> tvs_;
  Teuchos::Array<Teuchos::Array<Teuchos::RCP<Operator> > > blocks_;
  
  Teuchos::RCP<Epetra_CrsMatrix> A_;
  Teuchos::RCP<MatrixFE> Amat_;
  Teuchos::RCP<SuperMap> smap_;

  Teuchos::RCP<Matrix<TreeVector>> preconditioner_;
  bool block_diagonal_;

  int num_colors_;
  Teuchos::RCP<std::vector<int>> coloring_;
  Teuchos::ParameterList inv_plist_;
  bool inited_;

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
class TreeOperator_BlockPreconditioner {
 public:
  using Vector_t = TreeOperator::Vector_t;
  using VectorSpace_t = TreeOperator::VectorSpace_t;
  
  TreeOperator_BlockPreconditioner(TreeOperator& op) :
      op_(op) {}

  int Apply(const TreeVector& X, TreeVector& Y) const {
    Exceptions::amanzi_throw("TreeOperator Preconditioner does not implement Apply()");
    return 1;
  }
  
  int ApplyInverse(const TreeVector& X, TreeVector& Y) const {
    return op_.ApplyInverseBlockDiagonal_(X,Y);
  }
  void UpdateInverse() {
    for (int n=0; n!=op_.tvs_->size(); ++n) {
      op_.blocks_[n][n]->UpdateInverse();
    }
  }
  void ComputeInverse() {
    for (int n=0; n!=op_.tvs_->size(); ++n) {
      op_.blocks_[n][n]->ComputeInverse();
    }
  }
  
 private:
  TreeOperator& op_;
};

} // namespace Impl
  
}  // namespace Operators
}  // namespace Amanzi

#endif
