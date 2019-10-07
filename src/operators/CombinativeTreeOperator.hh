/*
  This is the Operator component of the Amanzi code.

  Copyright 2010-2013 held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#ifndef AMANZI_COMBINATIVE_TREE_OPERATOR_HH_
#define AMANZI_COMBINATIVE_TREE_OPERATOR_HH_

#include <vector>

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "TreeVector.hh"
#include "TreeVectorSpace.hh"
#include "TreeOperator.hh"
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

class CombinativeTreeOperator {
 public:
  CombinativeTreeOperator() {};
  CombinativeTreeOperator(const Teuchos::RCP<TreeOperator>& tree_op, Teuchos::RCP<Teuchos::ParameterList> plist, bool use_asc);
  ~CombinativeTreeOperator() {};

  virtual int Apply(const TreeVector& X, TreeVector& Y) const;
  virtual int ApplyInverse(const TreeVector& X, TreeVector& Y) const;

  void SymbolicAssembleMatrix();
  void AssembleMatrix();

  const TreeVectorSpace& DomainMap() const { return tree_op_->DomainMap(); }
  const TreeVectorSpace& RangeMap() const { return tree_op_->RangeMap(); }

  // preconditioners
  void InitPreconditionerGlobal(const std::string& prec_name, const Teuchos::ParameterList& plist);
  void InitPreconditionerBlock(const std::vector<std::string>& prec_names, const Teuchos::ParameterList& plist);

  // access
  Teuchos::RCP<Epetra_CrsMatrix> A() const { return tree_op_->A(); } 
  Teuchos::RCP<const TreeOperator> Op() const { return tree_op_; } 
  Teuchos::RCP<TreeOperator> Op() { return tree_op_; } 

  // deep copy for building interfaces to TPLs, mainly to solvers
  void CopyVectorToSuperVector(const TreeVector& cv, Epetra_Vector& sv) const {
    tree_op_->CopyVectorToSuperVector(cv, sv);
  }
  void CopySuperVectorToVector(const Epetra_Vector& sv, TreeVector& cv) const {
    tree_op_->CopySuperVectorToVector(sv, cv);
  }

 private:
  Teuchos::RCP<TreeOperator> tree_op_;
  bool global_solve_;
  std::vector<int> block_id_;
  bool use_asc_;

};

}  // namespace Operators
}  // namespace Amanzi


#endif


