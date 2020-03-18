/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@lanl.gov)
*/

#include <vector>

// TPLs
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

// Amanzi
#include "GraphFE.hh"
#include "MatrixFE.hh"
#include "PreconditionerFactory.hh"
#include "SuperMap.hh"
#include "VerboseObject.hh"

// Operators
#include "Operator.hh"
#include "OperatorUtils.hh"
#include "TreeOperator.hh"

namespace Amanzi {
namespace Operators {

/* ******************************************************************
* Constructor from a tree vector.
****************************************************************** */
TreeOperator::TreeOperator(Teuchos::RCP<const TreeVectorSpace> tvs) :
    tvs_(tvs),
    block_diagonal_(false)
{
  // make sure we have the right kind of TreeVectorSpace -- it should be
  // one parent node with all leaf node children.
  AMANZI_ASSERT(tvs_->Data() == Teuchos::null);
  for (TreeVectorSpace::const_iterator it = tvs_->begin(); it != tvs_->end(); ++it) {
    AMANZI_ASSERT((*it)->Data() != Teuchos::null);
  }

  // resize the blocks
  int n_blocks = tvs_->size();
  blocks_.resize(n_blocks, Teuchos::Array<Teuchos::RCP<const Operator> >(n_blocks, Teuchos::null));
}


/* ******************************************************************
* Populate block matrix with pointers to operators.
****************************************************************** */
void TreeOperator::SetOperatorBlock(int i, int j, const Teuchos::RCP<const Operator>& op) {
  blocks_[i][j] = op;
}

/* ******************************************************************
* Populate block matrix with pointers to operators.
****************************************************************** */
void TreeOperator::SetTreeOperatorBlock(int i, int j, const Teuchos::RCP<const TreeOperator>& op_tree) {
  int n_blocks = op_tree->GetNumberBlocks();
  for (int l=0; l<n_blocks; l++) {
    for (int m=0; m<n_blocks; m++) {
      blocks_[i+l][j+m] = op_tree->GetOperatorBlock(l,m);
    }
  }
}

  

/* ******************************************************************
* Calculate Y = A * X using matrix-free matvec on blocks of operators.
****************************************************************** */
int TreeOperator::Apply(const TreeVector& X, TreeVector& Y) const
{
  Y.PutScalar(0.0);

  int ierr(0), n(0);
  for (TreeVector::iterator yN_tv = Y.begin(); yN_tv != Y.end(); ++yN_tv, ++n) {
    CompositeVector& yN = *(*yN_tv)->Data();
    int m(0);
    for (TreeVector::const_iterator xM_tv = X.begin(); xM_tv != X.end(); ++xM_tv, ++m) {
      if (blocks_[n][m] != Teuchos::null) {
        ierr |= blocks_[n][m]->Apply(*(*xM_tv)->Data(), yN, 1.0);
      }
    }
  }
  return ierr;
}


/* ******************************************************************
* Calculate Y = A * X using matrix-free matvec on blocks of operators.
****************************************************************** */
int TreeOperator::ApplyAssembled(const TreeVector& X, TreeVector& Y) const
{
  Y.PutScalar(0.0);
  Epetra_Vector Xcopy(A_->RowMap());
  Epetra_Vector Ycopy(A_->RowMap());
  double x_norm, y_norm;

  int ierr = CopyTreeVectorToSuperVector(*smap_, X, Xcopy);

  ierr |= A_->Apply(Xcopy, Ycopy);

  ierr |= CopySuperVectorToTreeVector(*smap_, Ycopy, Y);
  AMANZI_ASSERT(!ierr);

  return ierr;
}


/* ******************************************************************
* Calculate Y = inv(A) * X using global matrix.
****************************************************************** */
int TreeOperator::ApplyInverse(const TreeVector& X, TreeVector& Y) const
{
  int code(0);
  if (!block_diagonal_) {
    Epetra_Vector Xcopy(A_->RowMap());
    Epetra_Vector Ycopy(A_->RowMap());

    int ierr = CopyTreeVectorToSuperVector(*smap_, X, Xcopy);
    code = preconditioner_->ApplyInverse(Xcopy, Ycopy);
    ierr |= CopySuperVectorToTreeVector(*smap_, Ycopy, Y);

    AMANZI_ASSERT(!ierr);
  } else {
    for (int n = 0; n < tvs_->size(); ++n) {
      const CompositeVector& Xn = *X.SubVector(n)->Data();
      CompositeVector& Yn = *Y.SubVector(n)->Data();
      code |= blocks_[n][n]->ApplyInverse(Xn, Yn);
    }
  }

  return code;
}

    
/* ******************************************************************
* Symbolic assemble global matrix from elemental matrices of block 
* operators. 
****************************************************************** */
void TreeOperator::SymbolicAssembleMatrix()
{
  int n_blocks = blocks_.size();

  // Currently we assume all diagonal schema are the same and well defined.
  // May be ways to relax this a bit in the future, but it currently covers
  // all uses.
  int schema = 0;
  std::vector<CompositeVectorSpace> cvs_vec;
  std::vector<std::string> cvs_names;

  // Check that each row has at least one non-null operator block
  // and save the position of this block, preferably diagonal.
  Teuchos::RCP<const Operator> an_op;
  for (int lcv_row = 0; lcv_row != n_blocks; ++lcv_row) {
    int block_col(-1);
    for (int lcv_col = 0; lcv_col != n_blocks; ++lcv_col) {
      if (blocks_[lcv_row][lcv_col] != Teuchos::null) {
        an_op = blocks_[lcv_row][lcv_col];
        if (block_col != lcv_row) block_col = lcv_col;
      }
    }
    AMANZI_ASSERT(block_col >= 0);

    cvs_vec.push_back(blocks_[lcv_row][block_col]->RangeMap());
    cvs_names.push_back(std::to_string(lcv_row));
  }

  // create the supermap and graph
  smap_ = createSuperMap(DomainMap());

  // NOTE: this probably needs to be fixed for differing meshes. -etc
  int row_size = MaxRowSize(*an_op->DomainMap().Mesh(), schema, n_blocks);
  auto graph = Teuchos::rcp(new GraphFE(smap_->Map(), 
      smap_->GhostedMap(), smap_->GhostedMap(), row_size));

  // fill the graph
  for (int lcv_row = 0; lcv_row != n_blocks; ++lcv_row) {
    for (int lcv_col = 0; lcv_col != n_blocks; ++lcv_col) {
      Teuchos::RCP<const Operator> block = blocks_[lcv_row][lcv_col];
      if (block != Teuchos::null) {
        block->SymbolicAssembleMatrix(*smap_, *graph, lcv_row, lcv_col);
      }
    }
  }

  // assemble the graph
  int ierr = graph->FillComplete(smap_->Map(), smap_->Map());
  AMANZI_ASSERT(!ierr);

  // create the matrix
  Amat_ = Teuchos::rcp(new MatrixFE(graph));
  A_ = Amat_->Matrix();
}


/* ******************************************************************
* Assemble global matrix from elemental matrices of block operators.
****************************************************************** */
void TreeOperator::AssembleMatrix() {
  int n_blocks = blocks_.size();
  Amat_->Zero();

  // check that each row has at least one non-null operator block
  for (int lcv_row = 0; lcv_row != n_blocks; ++lcv_row) {
    for (int lcv_col = 0; lcv_col != n_blocks; ++lcv_col) {
      Teuchos::RCP<const Operator> block = blocks_[lcv_row][lcv_col];
      if (block != Teuchos::null) {
        block->AssembleMatrix(*smap_, *Amat_, lcv_row, lcv_col);
      }
    }
  }

  int ierr = Amat_->FillComplete();
  AMANZI_ASSERT(!ierr);

  // std::stringstream filename_s2;
  // filename_s2 << "assembled_matrix" << 0 << ".txt";
  // EpetraExt::RowMatrixToMatlabFile(filename_s2.str().c_str(), *Amat_->Matrix());
}


/* ******************************************************************
* Create preconditioner using name and a factory.
****************************************************************** */
void TreeOperator::InitPreconditioner(
    const std::string& prec_name, const Teuchos::ParameterList& plist)
{
  AmanziPreconditioners::PreconditionerFactory factory;
  preconditioner_ = factory.Create(prec_name, plist);
  preconditioner_->Update(A_);
}


/* ******************************************************************
* Create preconditioner using name and a factory.
****************************************************************** */
void TreeOperator::InitPreconditioner(Teuchos::ParameterList& plist)
{
  // provide block ids for block strategies.
  if (plist.isParameter("preconditioner type") &&
      plist.get<std::string>("preconditioner type") == "boomer amg" &&
      plist.isSublist("boomer amg parameters")) {

    auto block_ids = smap_->BlockIndices();
    plist.sublist("boomer amg parameters").set("number of unique block indices", block_ids.first);
    plist.sublist("boomer amg parameters").set("block indices", block_ids.second);
  }

  AmanziPreconditioners::PreconditionerFactory factory;
  preconditioner_ = factory.Create(plist);
  preconditioner_->Update(A_);
}


/* ******************************************************************
* Two-stage initialization of preconditioner, part 1.
* Create the PC and set options.  SymbolicAssemble() must have been called.
****************************************************************** */
void TreeOperator::InitializePreconditioner(Teuchos::ParameterList& plist)
{
  AMANZI_ASSERT(A_.get());
  AMANZI_ASSERT(smap_.get());

  // provide block ids for block strategies.
  if (plist.isParameter("preconditioner type") &&
      plist.get<std::string>("preconditioner type") == "boomer amg" &&
      plist.isSublist("boomer amg parameters")) {

    // NOTE: Hypre frees this
    auto block_ids = smap_->BlockIndices();

    plist.sublist("boomer amg parameters").set("number of unique block indices", block_ids.first);

    // Note, this passes a raw pointer through a ParameterList.  I was surprised
    // this worked too, but ParameterList is a boost::any at heart... --etc
    plist.sublist("boomer amg parameters").set("block indices", block_ids.second);
  }

  AmanziPreconditioners::PreconditionerFactory factory;
  preconditioner_ = factory.Create(plist);
}


/* ******************************************************************
* Two-stage initialization of preconditioner, part 2.
* Set the matrix in the preconditioner.  Assemble() must have been called.
****************************************************************** */
void TreeOperator::UpdatePreconditioner()
{
  AMANZI_ASSERT(preconditioner_.get());
  AMANZI_ASSERT(A_.get());
  preconditioner_->Update(A_);
}


/* ******************************************************************
* Init block-diagonal preconditioner
****************************************************************** */
void TreeOperator::InitBlockDiagonalPreconditioner()
{
  block_diagonal_ = true;
}
 

/* ******************************************************************
* Copies to/from SuperVector for use by Amesos.
****************************************************************** */
void TreeOperator::CopyVectorToSuperVector(const TreeVector& cv, Epetra_Vector& sv) const
{
  CopyTreeVectorToSuperVector(*smap_, cv, sv);
}


void TreeOperator::CopySuperVectorToVector(const Epetra_Vector& sv, TreeVector& cv) const
{
  CopySuperVectorToTreeVector(*smap_, sv, cv);
}

}  // namespace Operators
}  // namespace Amanzi


