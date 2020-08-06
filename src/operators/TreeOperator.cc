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
#include "SuperMap.hh"
#include "VerboseObject.hh"

// Operators
#include "Op.hh"
#include "Operator.hh"
#include "OperatorUtils.hh"
#include "TreeOperator.hh"

#define TEST_MAPS 0

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
  blocks_.resize(n_blocks, Teuchos::Array<Teuchos::RCP<Operator> >(n_blocks, Teuchos::null));
}


/* ******************************************************************
* Populate block matrix with pointers to operators.
****************************************************************** */
void TreeOperator::SetOperatorBlock(int i, int j, const Teuchos::RCP<Operator>& op) {
  blocks_[i][j] = op;
}


/* ******************************************************************
* Calculate Y = A * X using matrix-free matvec on blocks of operators.
****************************************************************** */
int TreeOperator::Apply(const TreeVector& X, TreeVector& Y) const
{
#if TEST_MAPS
  AMANZI_ASSERT(DomainMap().SubsetOf(X.Map()));
  AMANZI_ASSERT(RangeMap().SubsetOf(Y.Map()));
#endif  
  Y.PutScalar(0.0);
  
  int ierr(0);
  for (int n=0; n!=Y.size(); ++n) {
    AMANZI_ASSERT(Y.SubVector(n)->Data() != Teuchos::null); // only works on 1-level heirarchy currently, see #453
    CompositeVector& yN = *Y.SubVector(n)->Data();
    for (int m=0; m!=X.size(); ++m) {
      if (blocks_[n][m] != Teuchos::null) {
	AMANZI_ASSERT(X.SubVector(m)->Data() != Teuchos::null); // only works on 1-level heirarchy currently, see #453
	const CompositeVector& xM = *X.SubVector(m)->Data();
        ierr |= blocks_[n][m]->Apply(xM, yN, 1.0);
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

  int ierr = copyToSuperVector(*smap_, X, Xcopy);

  ierr |= A_->Apply(Xcopy, Ycopy);

  ierr |= copyFromSuperVector(*smap_, Ycopy, Y);
  AMANZI_ASSERT(!ierr);

  return ierr;
}


/* ******************************************************************
* Calculate Y = inv(A) * X using global matrix.
****************************************************************** */
int TreeOperator::ApplyInverse(const TreeVector& Y, TreeVector& X) const
{
#if TEST_MAPS
  AMANZI_ASSERT(DomainMap().SubsetOf(X.Map()));
  AMANZI_ASSERT(RangeMap().SubsetOf(Y.Map()));
#endif
  if (preconditioner_.get()) {
    return preconditioner_->ApplyInverse(Y, X);
  } else {
    AMANZI_ASSERT(block_diagonal_);  // this assertion shouldn't be possible --
                                     // in all cases where block_diagonal_
                                     // isn't true, a preconditioner_ should
                                     // have been created.
    return ApplyInverseBlockDiagonal_(Y,X);
  }
}


/* ******************************************************************
* Calculate Y = inv(A) * X using the block diagonal
****************************************************************** */
int TreeOperator::ApplyInverseBlockDiagonal_(const TreeVector& Y, TreeVector& X) const
{
  int code = 0;
  for (int n = 0; n < tvs_->size(); ++n) {
    const CompositeVector& Yn = *Y.SubVector(n)->Data();
    CompositeVector& Xn = *X.SubVector(n)->Data();
    code |= blocks_[n][n]->ApplyInverse(Yn, Xn);
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


void TreeOperator:: InitializeInverse(const std::string& prec_name,
        const Teuchos::ParameterList& plist) {
  Teuchos::ParameterList inner_plist(plist.sublist(prec_name));
  InitializeInverse(inner_plist);
}  

/* ******************************************************************
* Two-stage initialization of preconditioner, part 1.
* Create the PC and set options.  SymbolicAssemble() must have been called.
****************************************************************** */
void TreeOperator::InitializeInverse(Teuchos::ParameterList& plist)
{
  // // provide block ids for block strategies.
  // if (plist.isParameter("preconditioning method") &&
  //     plist.get<std::string>("preconditioning method") == "boomer amg" &&
  //     plist.isSublist("boomer amg parameters")) {

  //   // NOTE: Hypre frees this
  //   auto block_ids = smap_->BlockIndices();

  //   plist.sublist("boomer amg parameters").set("number of unique block indices", block_ids.first);

  //   // Note, this passes a raw pointer through a ParameterList.  I was surprised
  //   // this worked too, but ParameterList is a boost::any at heart... --etc
  //   plist.sublist("boomer amg parameters").set("block indices", block_ids.second);
  // }

  // create the inverse
  if (plist.isParameter("preconditioning method")) {
    if (plist.get<std::string>("preconditioning method") == "block diagonal") {
      // are we using block diagonal preconditioning?  If so, this provides the
      // preconditioner...
      if (plist.isParameter("iterative method")) {
        // are we wrapping it in an iterative method?
        //
        // If so, we have to wrap this so that there are two accessible
        // ApplyInverse() methods -- one, provided by this, which gives the
        // full iterative method inverse to operator clients.  The second,
        // which is used by the iterative method, is the wrapped one which just
        // calls this's ApplyInverseBlockDiagional_ method.
        auto pc = Teuchos::rcp(new Impl::TreeOperator_BlockPreconditioner(*this));
        preconditioner_ = AmanziSolvers::createIterativeMethod(plist, Teuchos::rcpFromRef(*this), pc);
      } else {
        block_diagonal_ = true;
        preconditioner_ = Teuchos::null;
      }
    } else {
      // probably an assembled inverse method, with or without an iterative
      // method on top.  Call the factory, which assumes this is an assembler
      // factory.
      preconditioner_ = AmanziSolvers::createInverse(plist, Teuchos::rcpFromRef(*this));
    }
  } else {
    // probably a direct method...
    preconditioner_ = AmanziSolvers::createInverse(plist, Teuchos::rcpFromRef(*this));
  }
}


/* ******************************************************************
* Two-stage initialization of preconditioner, part 2.
* Set the matrix in the preconditioner.  Assemble() must have been called.
****************************************************************** */
void TreeOperator::UpdateInverse()
{
  if (preconditioner_.get()) {
    preconditioner_->UpdateInverse(); // calls SymbolicAssemble if needed
  } else if (block_diagonal_) {
    for (int n=0; n!=tvs_->size(); ++n) {
      blocks_[n][n]->UpdateInverse();
    }
  }
    
}

void TreeOperator::ComputeInverse()
{
  if (preconditioner_.get()) {
    preconditioner_->ComputeInverse(); // calls SymbolicAssemble if needed
  } else if (block_diagonal_) {
    for (int n=0; n!=tvs_->size(); ++n) {
      blocks_[n][n]->ComputeInverse();
    }
  }
}


/* ******************************************************************
* Populates matrix entries.
****************************************************************** */
std::string TreeOperator::PrintDiagnostics() const
{
  std::stringstream msg;
  int n_blocks = blocks_.size();
  for (int i = 0; i < n_blocks; ++i) {
    for (int j = 0; j < n_blocks; ++j) {
      auto block = blocks_[i][j];
      if (block != Teuchos::null) {
        msg << " block " << i << " " << j << ": ";
        for (auto it = block->begin(); it != block->end(); ++it) {
          msg << "<" << (*it)->schema_string << "> ";
        }
        msg << "\n";
      }
    }
  }
  return msg.str();
}

}  // namespace Operators
}  // namespace Amanzi
