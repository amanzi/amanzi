/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon (coonet@ornl.gov)
*/

//! <MISSING_ONELINE_DOCSTRING>

#define ASSEMBLY_DONE 0


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
TreeOperator::TreeOperator(Teuchos::RCP<const TreeVectorSpace> tvs)
  : tvs_(tvs), block_diagonal_(false)
{
  // make sure we have the right kind of TreeVectorSpace -- it should be
  // one parent node with all leaf node children.
  AMANZI_ASSERT(tvs_->Data() == Teuchos::null);
  for (TreeVectorSpace::const_iterator it = tvs_->begin(); it != tvs_->end();
       ++it) {
    AMANZI_ASSERT((*it)->Data() != Teuchos::null);
  }

  // resize the blocks
  int n_blocks = tvs_->size();
  blocks_.resize(
    n_blocks,
    Teuchos::Array<Teuchos::RCP<const Operator>>(n_blocks, Teuchos::null));
  transpose_.resize(n_blocks, Teuchos::Array<bool>(n_blocks, false));
}


/* ******************************************************************
 * Populate block matrix with pointers to operators.
 ****************************************************************** */
void
TreeOperator::SetOperatorBlock(int i, int j,
                               const Teuchos::RCP<const Operator>& op,
                               bool transpose)
{
  blocks_[i][j] = op;
  transpose_[i][j] = transpose;
}


void
TreeOperator::getLocalDiagCopy(TreeVector& tv) const
{
  int i=0;
  for (auto sv : tv) blocks_[i][i]->getLocalDiagCopy(*sv->Data());
}


/* ******************************************************************
 * Calculate Y = A * X using matrix-free matvec on blocks of operators.
 ****************************************************************** */
int
TreeOperator::apply(const TreeVector& X, TreeVector& Y) const
{
  Y.putScalar(0.0);

  int ierr(0), n(0);
  for (TreeVector::iterator yN_tv = Y.begin(); yN_tv != Y.end(); ++yN_tv, ++n) {
    CompositeVector& yN = *(*yN_tv)->Data();
    int m(0);
    for (TreeVector::const_iterator xM_tv = X.begin(); xM_tv != X.end();
         ++xM_tv, ++m) {
      if (blocks_[n][m] != Teuchos::null) {
        if (transpose_[n][m]) {
          ierr |= blocks_[n][m]->applyTranspose(*(*xM_tv)->Data(), yN, 1.0);
        } else {
          ierr |= blocks_[n][m]->apply(*(*xM_tv)->Data(), yN, 1.0);
        }
      }
    }
  }
  return ierr;
}


/* ******************************************************************
 * Calculate Y = A * X using matrix-free matvec on blocks of operators.
 ****************************************************************** */
int
TreeOperator::applyAssembled(const TreeVector& X, TreeVector& Y) const
{
#if ASSEMBLY_DONE
  Y.putScalar(0.0);
  Vector_type Xcopy(A_->getRowMap());
  Vector_type Ycopy(A_->getRowMap());
  double x_norm, y_norm;

  int ierr = copyToSuperVector(*smap_, X, Xcopy);
  A_->apply(Xcopy, Ycopy);
  ierr |= copyFromSuperVector(*smap_, Ycopy, Y);
  AMANZI_ASSERT(!ierr);
  return ierr;
#else 
  return -1; 
#endif
}


/* ******************************************************************
 * Calculate Y = inv(A) * X using global matrix.
 ****************************************************************** */
int
TreeOperator::applyInverse(const TreeVector& X, TreeVector& Y) const
{
  int code(0);
  if (!block_diagonal_) {
    code = preconditioner_->applyInverse(X, Y);
  } else {
    for (int n = 0; n < tvs_->size(); ++n) {
      const CompositeVector& Xn = *X.SubVector(n)->Data();
      CompositeVector& Yn = *Y.SubVector(n)->Data();
      code |= blocks_[n][n]->applyInverse(Xn, Yn);
    }
  }
  if (code) {
    Errors::Message msg("TreeOperator: ApplyInverse failed.\n");
    throw(msg);
  }
  return code;
}


/* ******************************************************************
 * Symbolic assemble global matrix from elemental matrices of block
 * operators.
 ****************************************************************** */
void
TreeOperator::SymbolicAssembleMatrix()
{
#if ASSEMBLY_DONE  
  int n_blocks = blocks_.size();

  // Currently we assume all diagonal schema are the same and well defined.
  // May be ways to relax this a bit in the future, but it currently covers
  // all uses.
  int schema = 0;
  std::vector<CompositeVectorSpace> cvs_vec;
  std::vector<std::string> cvs_names;

  // Check that each row has at least one non-null operator block
  Teuchos::RCP<const Operator> an_op;
  for (int lcv_row = 0; lcv_row != n_blocks; ++lcv_row) {
    bool is_block(false);
    for (int lcv_col = 0; lcv_col != n_blocks; ++lcv_col) {
      if (blocks_[lcv_row][lcv_col] != Teuchos::null) {
        an_op = blocks_[lcv_row][lcv_col];
        is_block = true;
      }

      if (lcv_row == lcv_col) {
        AMANZI_ASSERT(blocks_[lcv_row][lcv_col] != Teuchos::null);
        cvs_vec.push_back(blocks_[lcv_row][lcv_col]->getDomainMap());
        cvs_names.push_back(std::to_string(lcv_row));
      }
    }
    AMANZI_ASSERT(is_block);
  }

  // create the supermap and graph
  smap_ = createSuperMap(*getDomainMap());

  // NOTE: this probably needs to be fixed for differing meshes. -etc
  int row_size = MaxRowSize(*an_op->DomainMap().Mesh(), schema, n_blocks);
  auto graph = Teuchos::rcp(new GraphFE(
    smap_->getMap(), smap_->getGhostedMap(), smap_->getGhostedMap(), row_size));

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
  int ierr = graph->FillComplete(smap_->getMap(), smap_->getMap());
  AMANZI_ASSERT(!ierr);

  // create the matrix
  Amat_ = Teuchos::rcp(new MatrixFE(graph));
  A_ = Amat_->Matrix();
#endif
}


/* ******************************************************************
 * Assemble global matrix from elemental matrices of block operators.
 ****************************************************************** */
void
TreeOperator::AssembleMatrix()
{
#if ASSEMBLY_DONE
  if (!Amat_.get()) SymbolicAssembleMatrix();
  
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
#endif
}


/* ******************************************************************
 * Create preconditioner using name and a factory.
 ****************************************************************** */
void
TreeOperator::InitPreconditioner(const std::string& prec_name,
                                 const Teuchos::ParameterList& plist)
{
  AmanziPreconditioners::PreconditionerFactory<TreeOperator,TreeVector> factory;
  preconditioner_ = factory.Create(prec_name, plist);
  UpdatePreconditioner();
}


/* ******************************************************************
 * Create preconditioner using name and a factory.
 ****************************************************************** */
void
TreeOperator::InitPreconditioner(Teuchos::ParameterList& plist)
{
  AmanziPreconditioners::PreconditionerFactory<TreeOperator,TreeVector> factory;
  preconditioner_ = factory.Create(plist);
  UpdatePreconditioner();
}


/* ******************************************************************
 * Two-stage initialization of preconditioner, part 1.
 * Create the PC and set options.  SymbolicAssemble() must have been called.
 ****************************************************************** */
void
TreeOperator::InitializePreconditioner(Teuchos::ParameterList& plist)
{
  // // provide block ids for block strategies.
  // if (plist.isParameter("preconditioner type")) {
  //   std::string pc_type = plist.get<std::string>("preconditioner type");
  //   if (pc_type == "boomer amg") {
  //     // NOTE: Hypre frees this
  //     auto block_ids = smap_->BlockIndices();

  //     plist.sublist("boomer amg parameters")
  //         .set("number of unique block indices", block_ids.first);

  //     // Note, this passes a raw pointer through a ParameterList.  I was surprised
  //     // this worked too, but ParameterList is a boost::any at heart... --etc
  //     plist.sublist("boomer amg parameters")
  //         .set("block indices", block_ids.second);

  //   } else if (pc_type == "block diagonal") {
  //     InitBlockDiagonalPreconditioner();
  //     return;
  //   }
  // }

  AmanziPreconditioners::PreconditionerFactory<TreeOperator,TreeVector> factory;
  preconditioner_ = factory.Create(plist);
}


/* ******************************************************************
 * Two-stage initialization of preconditioner, part 2.
 * Set the matrix in the preconditioner.  Assemble() must have been called.
 ****************************************************************** */
void
TreeOperator::UpdatePreconditioner()
{
  if (!block_diagonal_) {
    AMANZI_ASSERT(preconditioner_.get());
    preconditioner_->Update(Teuchos::rcpFromRef(*this));
  }
}


/* ******************************************************************
 * Init block-diagonal preconditioner
 ****************************************************************** */
void
TreeOperator::InitBlockDiagonalPreconditioner()
{
  block_diagonal_ = true;
}


} // namespace Operators
} // namespace Amanzi
