/*
  This is the Operator component of the Amanzi code.

  Copyright 2010-2013 held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include <vector>

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "VerboseObject.hh"
#include "SuperMap.hh"
#include "MatrixFE.hh"
#include "GraphFE.hh"
#include "PreconditionerFactory.hh"

#include "Operator.hh"
#include "OperatorUtils.hh"
#include "TreeOperator.hh"

/* ******************************************************************
TreeOperators are the block analogue of Operators -- they provide a linear
operator acting on a TreeVectorSpace.  They are currently assumed R^n -> R^n,
and furthermore each block is currently assumed to be from R^m --> R^m for n =
i*m where i is an integer (every block's space is the same).

Future work will relax this constraint, but currently this can be used for
things like multi-phased flow, thermal+flow, etc.

****************************************************************** */ 

namespace Amanzi {
namespace Operators {

TreeOperator::TreeOperator(Teuchos::RCP<const TreeVectorSpace> tvs) :
    tvs_(tvs) {

  // make sure we have the right kind of TreeVectorSpace -- it should be
  // one parent node with all leaf node children.
  ASSERT(tvs_->Data() == Teuchos::null);
  for (std::vector<Teuchos::RCP<const TreeVectorSpace> >::const_iterator
           it = tvs_->SubVectors().begin();
       it != tvs_->SubVectors().end(); ++it) {
    ASSERT((*it)->Data() != Teuchos::null);
  }

  // resize the blocks
  int n_blocks = tvs_->SubVectors().size();
  blocks_.resize(n_blocks, Teuchos::Array<Teuchos::RCP<const Operator> >(n_blocks, Teuchos::null));
}


void
TreeOperator::SetOperatorBlock(int i, int j, const Teuchos::RCP<const Operator>& op) {
  blocks_[i][j] = op;
}


int
TreeOperator::Apply(const TreeVector& X, TreeVector& Y) const {
  // apply is done matrix-free
  Y.PutScalar(0.0);

  int ierr(0);
  for (int n = 0; n != blocks_.size(); ++n) {
    CompositeVector& yN = *Y.SubVectors()[n]->Data();
    for (int m = 0; m != blocks_.size(); ++m) {
      if (blocks_[n][m] != Teuchos::null) {
        ierr |= blocks_[n][m]->Apply(*X.SubVectors()[m]->Data(), yN, 1.0);
      }
    }
  }
  return ierr;
}


int
TreeOperator::ApplyInverse(const TreeVector& X, TreeVector& Y) const {
  Epetra_Vector Xcopy(A_->RowMap());
  Epetra_Vector Ycopy(A_->RowMap());
  int ierr = CopyTreeVectorToSuperVector(*smap_, X, Xcopy);

  // std::cout << "r - staggered" << std::endl;
  // Xcopy.Print(std::cout);
  
  ierr |= preconditioner_->ApplyInverse(Xcopy, Ycopy);
  ierr |= CopySuperVectorToTreeVector(*smap_, Ycopy, Y);
  ASSERT(!ierr);
  return ierr;
}

    
void
TreeOperator::SymbolicAssembleMatrix() {
  int n_blocks = blocks_.size();

  // Currently we assume all diagonal schema are the same and well defined.
  // May be ways to relax this a bit in the future, but it currently covers
  // all uses.
  int schema = 0;

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
        ASSERT(blocks_[lcv_row][lcv_col] != Teuchos::null);
        if (schema == 0) {
          schema = blocks_[lcv_row][lcv_col]->schema();
        } else {
          ASSERT(schema == blocks_[lcv_row][lcv_col]->schema());
        }
      }
    }
    ASSERT(is_block);
  }

  // create the supermap and graph
  smap_ = createSuperMap(an_op->DomainMap(), schema, n_blocks);
  int row_size = MaxRowSize(*an_op->DomainMap().Mesh(), schema, n_blocks);
  Teuchos::RCP<GraphFE> graph = Teuchos::rcp(new GraphFE(smap_->Map(),
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
  ASSERT(!ierr);

  // create the matrix
  Amat_ = Teuchos::rcp(new MatrixFE(graph));
  A_ = Amat_->Matrix();
}


void
TreeOperator::AssembleMatrix() {
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
  ASSERT(!ierr);
}


// preconditioners
void
TreeOperator::InitPreconditioner(const std::string& prec_name, const Teuchos::ParameterList& plist) {
  AmanziPreconditioners::PreconditionerFactory factory;
  preconditioner_ = factory.Create(prec_name, plist);
  preconditioner_->Update(A_);
}


// preconditioners
void
TreeOperator::InitPreconditioner(Teuchos::ParameterList& plist) {
  AmanziPreconditioners::PreconditionerFactory factory;
  preconditioner_ = factory.Create(plist);
  preconditioner_->Update(A_);
}

}  // namespace Operators
}  // namespace Amanzi


