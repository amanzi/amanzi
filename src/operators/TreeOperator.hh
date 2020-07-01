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
  TreeOperator() : block_diagonal_(false) {};
  TreeOperator(Teuchos::RCP<const TreeVectorSpace> tvs);
  TreeOperator(Teuchos::RCP<const TreeVectorSpace> tvs, int nblocks);
  virtual ~TreeOperator() = default;

  // main members
  void SetOperatorBlock(int i, int j, const Teuchos::RCP<const Operator>& op);
  void SetTreeOperatorBlock(int i, int j, const Teuchos::RCP<const TreeOperator>& op_tree);  
  
  virtual int Apply(const TreeVector& X, TreeVector& Y) const;
  virtual int ApplyAssembled(const TreeVector& X, TreeVector& Y) const;
  virtual int ApplyInverse(const TreeVector& X, TreeVector& Y) const;
  int ApplyFlattened(const TreeVector& X, TreeVector& Y) const;

  void SymbolicAssembleMatrix();
  void AssembleMatrix();

  const TreeVectorSpace& DomainMap() const { return *tvs_; }
  const TreeVectorSpace& RangeMap() const { return *tvs_; }

  // preconditioners
  void InitPreconditioner(const std::string& prec_name, const Teuchos::ParameterList& plist);
  void InitPreconditioner(Teuchos::ParameterList& plist);
  void InitBlockDiagonalPreconditioner();

  void InitializePreconditioner(Teuchos::ParameterList& plist);
  void UpdatePreconditioner();

  // access
  Teuchos::RCP<Epetra_CrsMatrix> A() { return A_; } 
  Teuchos::RCP<const Epetra_CrsMatrix> A() const { return A_; } 

  // deep copy for building interfaces to TPLs, mainly to solvers
  void CopyVectorToSuperVector(const TreeVector& cv, Epetra_Vector& sv) const;
  void CopySuperVectorToVector(const Epetra_Vector& sv, TreeVector& cv) const;
  Teuchos::RCP<const Operator> GetOperatorBlock(int i, int j) const { return blocks_[i][j];}
  int GetNumberBlocks() const {return blocks_.size();}

  // i/o
  std::string PrintDiagnostics() const;

 private:
  Teuchos::RCP<const TreeVectorSpace> tvs_;
  Teuchos::Array<Teuchos::Array<Teuchos::RCP<const Operator> > > blocks_;
  
  Teuchos::RCP<Epetra_CrsMatrix> A_;
  Teuchos::RCP<MatrixFE> Amat_;
  Teuchos::RCP<SuperMap> smap_;

  Teuchos::RCP<AmanziPreconditioners::Preconditioner> preconditioner_;
  bool block_diagonal_;

  Teuchos::RCP<VerboseObject> vo_;
};

}  // namespace Operators
}  // namespace Amanzi


#endif


