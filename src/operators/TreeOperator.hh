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

class TreeOperator {
 public:
  TreeOperator() : block_diagonal_(false){};
  TreeOperator(Teuchos::RCP<const TreeVectorSpace> tvs);

  // main members
  void SetOperatorBlock(int i, int j, const Teuchos::RCP<const Operator>& op);

  void getLocalDiagCopy(TreeVector& tv) const;

  int apply(const TreeVector& X, TreeVector& Y) const;
  int applyAssembled(const TreeVector& X, TreeVector& Y) const;
  int applyInverse(const TreeVector& X, TreeVector& Y) const;

  void SymbolicAssembleMatrix();
  void AssembleMatrix();

  Teuchos::RCP<const TreeVectorSpace> getDomainMap() const { return tvs_; }
  Teuchos::RCP<const TreeVectorSpace> getRangeMap() const { return tvs_; }
  Teuchos::RCP<const TreeVectorSpace> getRowMap() const { return tvs_; }  

  // preconditioners (deprecate)
  void InitPreconditioner(const std::string& prec_name,
                          const Teuchos::ParameterList& plist);
  void InitPreconditioner(Teuchos::ParameterList& plist);
  void InitBlockDiagonalPreconditioner() { block_diagonal_ = true; }


  // two-stage initializeation (preferred)
  void InitializePreconditioner(Teuchos::ParameterList& plist);
  void UpdatePreconditioner();

  // access
  Teuchos::RCP<Matrix_type> A() { return A_; }
  Teuchos::RCP<const Matrix_type> A() const { return A_; }
  Teuchos::RCP<const SuperMap> getSuperMap() const { return smap_; }

 private:
  Teuchos::RCP<const TreeVectorSpace> tvs_;
  Teuchos::Array<Teuchos::Array<Teuchos::RCP<const Operator>>> blocks_;

  Teuchos::RCP<Matrix_type> A_;
  Teuchos::RCP<MatrixFE> Amat_;
  Teuchos::RCP<SuperMap> smap_;

  Teuchos::RCP<AmanziPreconditioners::Preconditioner<TreeOperator,TreeVector> > preconditioner_;
  bool block_diagonal_;

  Teuchos::RCP<VerboseObject> vo_;
};

} // namespace Operators
} // namespace Amanzi


#endif
