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
#include "InverseFactory.hh"
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
QW}


/* ******************************************************************
* Populate block matrix with pointers to operators.
****************************************************************** */
void TreeOperator::SetOperatorBlock(int i, int j, const Teuchos::RCP<const Operator>& op) {
  blocks_[i][j] = op;
}

void TreeOperator::getLocalDiagCopy(TreeVector& tv) const {
  int i=0;
  for (const auto& tv_b : tv) {
    AMANZI_ASSERT(blocks_[i][i] != Teuchos::null);
    blocks_[i][i]->getLocalDiagCopy(*tv_b->Data());
  }    
}


/* ******************************************************************
* Calculate Y = A * X using matrix-free matvec on blocks of operators.
****************************************************************** */
int TreeOperator::apply(const TreeVector& X, TreeVector& Y) const
{
  Y.putScalar(0.0);

  int ierr(0), n(0);
  for (TreeVector::iterator yN_tv = Y.begin(); yN_tv != Y.end(); ++yN_tv, ++n) {
    CompositeVector& yN = *(*yN_tv)->Data();
    int m(0);
    for (TreeVector::const_iterator xM_tv = X.begin(); xM_tv != X.end(); ++xM_tv, ++m) {
      if (blocks_[n][m] != Teuchos::null) {
        ierr |= blocks_[n][m]->apply(*(*xM_tv)->Data(), yN, 1.0);
      }
    }
  }
  return ierr;
}


/* ******************************************************************
* Calculate Y = A * X using matrix-free matvec on blocks of operators.
****************************************************************** */
int TreeOperator::applyAssembled(const TreeVector& X, TreeVector& Y) const
{
  Y.putScalar(0.0);
  Vector_type Xcopy(A_->getRowMap());
  Vector_type Ycopy(A_->getRowMap());
  double x_norm, y_norm;

  int ierr = copyToSuperVector(*smap_, X, Xcopy);
  A_->apply(Xcopy, Ycopy);
  ierr |= copyFromSuperVector(*smap_, Ycopy, Y);

  if (ierr) {
    Errors::Message msg;
    msg << "TreeOperator: ApplyAssemble failed.\n";
    Exceptions::amanzi_throw(msg);
  }

  return ierr;
}


/* ******************************************************************
* Calculate Y = inv(A) * X using global matrix.
****************************************************************** */
int TreeOperator::applyInverse(const TreeVector& X, TreeVector& Y) const
{
  int code(0);
  if (!block_diagonal_) {
    if (preconditioner_.get() == nullptr) {
      Errors::Message msg("TreeOperator did not initialize a preconditioner.\n");
      Exceptions::amanzi_throw(msg);
    }
    int ierr = preconditioner_->applyInverse(X, Y);
    if (ierr) {
      Errors::Message msg("TreeOperator: applyInverse failed.\n");
      Exceptions::amanzi_throw(msg);
    }
    return ierr;

  } else {
    for (int n = 0; n < tvs_->size(); ++n) {
      const CompositeVector& Xn = *X.SubVector(n)->Data();
      CompositeVector& Yn = *Y.SubVector(n)->Data();
      code |= blocks_[n][n]->applyInverse(Xn, Yn);
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
  // NOTE: not yet implemented, as we have no preconditioner, and the user
  // shouldn't call this anyway (only the PC should!)
  //
  // What still needs to be implemented is GraphFE and MatrixFE (or ditch if no
  // longer necessary) and the SymbolicAssemble sub-calls.
  AMANZI_ASSERT(false);

  
  // int n_blocks = blocks_.size();

  // // Currently we assume all diagonal schema are the same and well defined.
  // // May be ways to relax this a bit in the future, but it currently covers
  // // all uses.
  // int schema = 0;
  // std::vector<CompositeVectorSpace> cvs_vec;
  // std::vector<std::string> cvs_names;

  // // Check that each row has at least one non-null operator block
  // // and save the position of this block, preferably diagonal.
  // Teuchos::RCP<const Operator> an_op;
  // for (int lcv_row = 0; lcv_row != n_blocks; ++lcv_row) {
  //   int block_col(-1);
  //   for (int lcv_col = 0; lcv_col != n_blocks; ++lcv_col) {
  //     if (blocks_[lcv_row][lcv_col] != Teuchos::null) {
  //       an_op = blocks_[lcv_row][lcv_col];
  //       if (block_col != lcv_row) block_col = lcv_col;
  //     }
  //   }
  //   AMANZI_ASSERT(block_col >= 0);

  //   cvs_vec.push_back(blocks_[lcv_row][block_col]->RangeMap());
  //   cvs_names.push_back(std::to_string(lcv_row));
  // }

  // // create the supermap and graph
  // smap_ = createSuperMap(*getDomainMap());

  // // NOTE: this probably needs to be fixed for differing meshes. -etc
  // int row_size = MaxRowSize(*an_op->DomainMap().Mesh(), schema, n_blocks);
  // auto graph = Teuchos::rcp(new GraphFE(smap_->getMap(), 
  //     smap_->getGhostedMap(), smap_->getGhostedMap(), row_size));

  // // fill the graph
  // for (int lcv_row = 0; lcv_row != n_blocks; ++lcv_row) {
  //   for (int lcv_col = 0; lcv_col != n_blocks; ++lcv_col) {
  //     Teuchos::RCP<const Operator> block = blocks_[lcv_row][lcv_col];
  //     if (block != Teuchos::null) {
  //       block->SymbolicAssembleMatrix(*smap_, *graph, lcv_row, lcv_col);
  //     }
  //   }
  // }

  // // assemble the graph
  // int ierr = graph->FillComplete(smap_->Map(), smap_->Map());
  // AMANZI_ASSERT(!ierr);

  // // create the matrix
  // Amat_ = Teuchos::rcp(new MatrixFE(graph));
  // A_ = Amat_->Matrix();
}


/* ******************************************************************
* Assemble global matrix from elemental matrices of block operators.
****************************************************************** */
void TreeOperator::AssembleMatrix() {
  // NOTE: not yet implemented, as we have no preconditioner, and the user
  // shouldn't call this anyway (only the PC should!)
  //
  // What still needs to be implemented is GraphFE and MatrixFE (or ditch if no
  // longer necessary) and the SymbolicAssemble sub-calls.
  AMANZI_ASSERT(false);

  // int n_blocks = blocks_.size();
  // Amat_->Zero();

  // // check that each row has at least one non-null operator block
  // for (int lcv_row = 0; lcv_row != n_blocks; ++lcv_row) {
  //   for (int lcv_col = 0; lcv_col != n_blocks; ++lcv_col) {
  //     Teuchos::RCP<const Operator> block = blocks_[lcv_row][lcv_col];
  //     if (block != Teuchos::null) {
  //       block->AssembleMatrix(*smap_, *Amat_, lcv_row, lcv_col);
  //     }
  //   }
  // }

  // int ierr = Amat_->FillComplete();
  // AMANZI_ASSERT(!ierr);
}



/* ******************************************************************
* Methods require to enable an Inverse
****************************************************************** */
void TreeOperator:: set_inverse_parameters(const std::string& prec_name,
        const Teuchos::ParameterList& plist) {
  Teuchos::ParameterList inner_plist(plist.sublist(prec_name));
  set_inverse_parameters(inner_plist);
}


void TreeOperator:: set_inverse_parameters(Teuchos::ParameterList& inv_plist) {
  inv_plist_ = inv_plist;
  inited_ = true;
}


void TreeOperator::initializeInverse()
{
  // if this is a leaf, just initialize the Operator
  if (get_operator() != Teuchos::null) {
    get_operator()->initializeInverse();
    return;
  }

  // otherwise this must be square
  AMANZI_ASSERT(IsSquare());
#if TEST_MAPS
  AMANZI_ASSERT(get_row_map()->SameAs(*get_col_map()));
#endif

  // and inited
  if (!inited_) {
    Errors::Message msg("Developer error: set_inverse_parameters() has not been called.");
    msg << " ref: " << typeid(*this).name() << "\n";
    Exceptions::amanzi_throw(msg);
  }

  // deal with preconditioner options that are local to this
  if (inv_plist_.isParameter("preconditioning method") &&
      inv_plist_.get<std::string>("preconditioning method") == "boomer amg" &&
      inv_plist_.isSublist("boomer amg parameters") &&
      inv_plist_.sublist("boomer amg parameters").get<bool>("use block indices", false)) {
    // Must do initialization here, because the parameter list carries the block
    // indices, which need structure.  Since not guaranteed structure until Initialize,
    // is called, we cannot set block indicies until now.
    // provide block ids for block strategies.
    if (!row_supermap_.get()) row_supermap_ = createSuperMap(*get_row_map());

    if (coloring_ == Teuchos::null || num_colors_ == 0) {
      auto block_ids = get_row_supermap()->BlockIndices();
      set_coloring(block_ids.first, block_ids.second);
    }
    inv_plist_.sublist("boomer amg parameters").set("number of unique block indices", num_colors_);
    inv_plist_.sublist("boomer amg parameters").set("block indices", coloring_);
  }

  // create the inverse
  if (inv_plist_.isParameter("preconditioning method") &&
      inv_plist_.get<std::string>("preconditioning method") == "block diagonal") {
    // are we using block diagonal preconditioning?  If so, this provides the
    // preconditioner...
    // check for an iterative method on top of the block diagional PC
    if (inv_plist_.isParameter("iterative method")) {
      // are we wrapping it in an iterative method?
      //
      // If so, we have to wrap this so that there are two accessible
      // ApplyInverse() methods -- one, provided by this, which gives the
      // full iterative method inverse to operator clients.  The second,
      // which is used by the iterative method, is the wrapped one which just
      // calls this's ApplyInverseBlockDiagional_ method.
      auto pc = Teuchos::rcp(new Impl::TreeOperator_BlockDiagonalPreconditioner(*this));
      preconditioner_ = AmanziSolvers::createIterativeMethod(inv_plist_, Teuchos::rcpFromRef(*this), pc);
    } else {
      block_diagonal_ = true;
      preconditioner_ = Teuchos::null;
    }
  } else {
    // Probably an assembled inverse method, with or without an iterative
    // method on top, or a direct method.  Call the factory, which assumes this is an assembler
    // factory.
    preconditioner_ = AmanziSolvers::createInverse(inv_plist_, Teuchos::rcpFromRef(*this));
  }

  // call the Initialize method
  if (preconditioner_.get()) {
    preconditioner_->InitializeInverse(); // calls SymbolicAssemble if needed
  } else if (block_diagonal_) {
    // initialize the diagonal
    for (std::size_t n=0; n!=get_row_map()->size(); ++n) {
      Teuchos::RCP<TreeOperator> block = get_block(n,n);
      if (block == Teuchos::null) {
        Errors::Message msg("TreeOperator: block diagonal preconditioner requested, but block (");
        msg << n << "," << n << ") has not been provided.";
        Exceptions::amanzi_throw(msg);
      }
      block->InitializeInverse();
    }
  }
  updated_ = true;
  computed_ = false;
}


void TreeOperator::computeInverse()
{
  if (!updated_) InitializeInverse();
  if (get_operator() != Teuchos::null) {
    // leaf operator
    get_operator()->ComputeInverse();
  } else {
    if (preconditioner_.get()) {
      preconditioner_->ComputeInverse(); // calls SymbolicAssemble if needed
    } else if (block_diagonal_) {
      for (std::size_t n=0; n!=get_row_map()->size(); ++n) {
        blocks_[n][n]->ComputeInverse();
      }
    }
  }
  computed_ = true;
}


}  // namespace Operators
}  // namespace Amanzi


