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
#include "TreeVector_Utils.hh"
#include "InverseFactory.hh"

#define TEST_MAPS 0

namespace Amanzi {
namespace Operators {

/* ******************************************************************
* Constructors
****************************************************************** */
TreeOperator::TreeOperator()
  : inited_(false),
    updated_(false),
    computed_(false),
    block_diagonal_(false),
    num_colors_(0),
    coloring_(Teuchos::null)
{
  vo_ = Teuchos::rcp(new VerboseObject("TreeOperator", Teuchos::ParameterList()));
}


TreeOperator::TreeOperator(Teuchos::ParameterList& plist)
  : TreeOperator()
{
  vo_ = Teuchos::rcp(new VerboseObject("TreeOperator", plist));
  if (plist.isSublist("inverse")) set_inverse_parameters(plist.sublist("inverse"));
}


TreeOperator::TreeOperator(const Teuchos::RCP<const TreeVectorSpace>& row_map,
                           const Teuchos::RCP<const TreeVectorSpace>& col_map,
                           Teuchos::ParameterList& plist)
  : TreeOperator(plist)
{
  col_map_ = col_map;
  row_map_ = row_map;
  col_size_ = col_map_->size() + (col_map_->Data() == Teuchos::null ? 0 : 1);
  row_size_ = row_map_->size() + (row_map_->Data() == Teuchos::null ? 0 : 1);

  // resize the blocks
  blocks_.resize(row_size_, Teuchos::Array<Teuchos::RCP<TreeOperator> >(col_size_, Teuchos::null));
}


TreeOperator::TreeOperator(const Teuchos::RCP<const TreeVectorSpace>& row_map,
                           const Teuchos::RCP<const TreeVectorSpace>& col_map)
  : TreeOperator()
{
  col_map_ = col_map;
  row_map_ = row_map;
  col_size_ = col_map_->size() + (col_map_->Data() == Teuchos::null ? 0 : 1);
  row_size_ = row_map_->size() + (row_map_->Data() == Teuchos::null ? 0 : 1);

  // resize the blocks
  blocks_.resize(row_size_, Teuchos::Array<Teuchos::RCP<TreeOperator> >(col_size_, Teuchos::null));
}

TreeOperator::TreeOperator(const Teuchos::RCP<const TreeVectorSpace>& tvs,
                           Teuchos::ParameterList& plist)
  : TreeOperator(tvs, tvs, plist) {}

TreeOperator::TreeOperator(const Teuchos::RCP<const TreeVectorSpace>& tvs)
  : TreeOperator(tvs, tvs) {}


/* ******************************************************************
* Copy constructor does a deep copy.
****************************************************************** */
TreeOperator::TreeOperator(const TreeOperator& other)
  : TreeOperator(other.row_map_, other.col_map_)
{
  vo_ = other.vo_;
  inv_plist_ = other.inv_plist_;

  for (int i=0; i!=row_size_; ++i) {
    for (int j=0; j!=col_size_; ++j) {
      if (other.blocks_[i][j] != Teuchos::null) {
        blocks_[i][j] = other.blocks_[i][j]->Clone();
      }
    }
  }

  if (other.data_ != Teuchos::null) {
    data_ = other.data_->Clone();
  }
}


/* ******************************************************************
* block setters
****************************************************************** */
void TreeOperator::set_block(std::size_t i, std::size_t j, const Teuchos::RCP<TreeOperator>& op)
{
  blocks_[i][j] = op;
  updated_ = false;
}


void TreeOperator::set_operator(const Teuchos::RCP<Operator>& op)
{
  data_ = op;
  updated_ = false;
}


void TreeOperator::set_operator_block(std::size_t i, std::size_t j, const Teuchos::RCP<Operator>& op) {
  auto block_row_map = get_row_map()->SubVector(i);
  auto block_col_map = get_col_map()->SubVector(j);

  // can only call this on a leaf block
  AMANZI_ASSERT(block_row_map->Data() != Teuchos::null);
  AMANZI_ASSERT(block_col_map->Data() != Teuchos::null); // can only call this on a leaf block
#if TEST_MAPS
  AMANZI_ASSERT(block_col_map->Data()->SameAs(op->get_domain_map()));
  AMANZI_ASSERT(block_row_map->Data()->SameAs(op->get_range_map()));
#endif

  Teuchos::ParameterList plist;
  auto top = Teuchos::rcp(new TreeOperator(block_row_map, block_col_map, plist));
  top->set_operator(op);
  set_block(i,j,top);
}


/* ******************************************************************
* Calculate Y = A * X using matrix-free matvec on blocks of operators.
****************************************************************** */
int TreeOperator::Apply(const TreeVector& X, TreeVector& Y, double scalar) const
{
#if TEST_MAPS
  AMANZI_ASSERT(get_domain_map().SubsetOf(X.Map()));
  AMANZI_ASSERT(get_range_map().SubsetOf(Y.Map()));
#endif

  int ierr(0);
  if (scalar == 0.0) {
    Y.PutScalarMasterAndGhosted(0.0);
  } else if (scalar != 1.0) {
    Y.Scale(scalar);
  }

  if (Y.Data() != Teuchos::null) {
    if (X.Data() != Teuchos::null) {
      ierr = get_operator()->Apply(*X.Data(), *Y.Data(), 1.0);
      if (ierr) return ierr;
    } else {
      for (std::size_t j=0; j!=X.size(); ++j) {
        const TreeVector& Xj = *X.SubVector(j);
        ierr = get_block(0,j)->Apply(Xj, Y, 1.0);
        if (ierr) return ierr;
      }
    }
  } else {
    for (std::size_t i=0; i!=Y.size(); ++i) {
      TreeVector& yi = *Y.SubVector(i);

      if (X.Data() != Teuchos::null) {
        ierr = get_block(i,0)->Apply(X, yi, 1.0);
        if (ierr) return ierr;
      } else {
        for (std::size_t j=0; j!=X.size(); ++j) {
          const TreeVector& xj = *X.SubVector(j);
          auto block = get_block(i,j);
          if (block != Teuchos::null) {
            ierr = block->Apply(xj, yi, 1.0);
            if (ierr) return ierr;
          }
        }
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
  Epetra_Vector Xcopy(A_->RangeMap());
  Epetra_Vector Ycopy(A_->DomainMap());

  int ierr = copyToSuperVector(*get_col_supermap(), X, Xcopy);
  ierr |= A_->Apply(Xcopy, Ycopy);
  ierr |= copyFromSuperVector(*get_row_supermap(), Ycopy, Y);
  AMANZI_ASSERT(!ierr);
  return ierr;
}


/* ******************************************************************
* Calculate Y = inv(A) * X using global matrix.
****************************************************************** */
int TreeOperator::ApplyInverse(const TreeVector& X, TreeVector& Y) const
{
#if TEST_MAPS
  AMANZI_ASSERT(get_domain_map().SubsetOf(Y.Map()));
  AMANZI_ASSERT(get_range_map().SubsetOf(X.Map()));
#endif

  int ierr = 0;
  if (get_operator() != Teuchos::null) {
    // this is a leaf
    ierr = get_operator()->ApplyInverse(*X.Data(), *Y.Data());
  } else {
    if (!computed_) const_cast<TreeOperator*>(this)->ComputeInverse();
    if (preconditioner_.get()) {
      ierr = preconditioner_->ApplyInverse(X,Y);
    } else {
      AMANZI_ASSERT(block_diagonal_);  // this assertion shouldn't be possible --
                                     // in all cases where block_diagonal_
                                     // isn't true, a preconditioner_ should
                                     // have been created.
      ierr = ApplyInverseBlockDiagonal_(X,Y);
    }
  }
  return ierr;
}


/* ******************************************************************
* Calculate Y = inv(A) * X using the block diagonal
****************************************************************** */
int TreeOperator::ApplyInverseBlockDiagonal_(const TreeVector& X, TreeVector& Y) const
{
  AMANZI_ASSERT(IsSquare());
#if TEST_MAPS
  AMANZI_ASSERT(get_domain_map().SameAs(get_range_map()));
  AMANZI_ASSERT(get_domain_map().SubsetOf(Y.Map()));
  AMANZI_ASSERT(get_range_map().SubsetOf(X.Map()));
#endif
  int ierr = 0;
  if (get_col_size() == 1) {
    AMANZI_ASSERT(get_operator() != Teuchos::null);
    AMANZI_ASSERT(X.Data() != Teuchos::null);
    AMANZI_ASSERT(Y.Data() != Teuchos::null);
    return get_operator()->ApplyInverse(*X.Data(), *Y.Data());

  } else {
    for (std::size_t i=0; i!=get_col_size(); ++i) {
      ierr = get_block(i,i)->ApplyInverse(*X.SubVector(i), *Y.SubVector(i));
      if (ierr) return ierr;
    }
  }
  return ierr;
}



/* ******************************************************************
* Symbolic assemble global matrix from elemental matrices of block
* operators.
****************************************************************** */
void TreeOperator::SymbolicAssembleMatrix()
{
  // create the supermaps
  if (!row_supermap_.get()) row_supermap_ = createSuperMap(*get_row_map());
  if (!col_supermap_.get()) col_supermap_ = createSuperMap(*get_col_map());

  // NOTE: this can be an overshoot, as we do this once, then FillComplete()
  // and clean up extra space.  From then on it is a static graph.  So there
  // should be no issue with overshooting, as long as it fits in memory once.
  //
  // Count the global number of leaves
  std::size_t n_row_leaves = getNumTreeVectorLeaves(*get_row_map());
  std::size_t n_col_leaves = getNumTreeVectorLeaves(*get_col_map());

  // allocate space, graph
  leaves_.resize(n_row_leaves, std::vector<Teuchos::RCP<Operator>>(n_col_leaves, Teuchos::null));
  int cols_per_row = n_col_leaves * OPERATOR_MAX_NUM_FACES;
  auto graph = Teuchos::rcp(new GraphFE(row_supermap_->Map(),
      row_supermap_->GhostedMap(), col_supermap_->GhostedMap(), cols_per_row));

  // get a flattened array of blocks that are leaves
  auto shape = Impl::collectTreeOperatorLeaves(*this, leaves_, 0, 0);
  AMANZI_ASSERT(n_row_leaves == shape.first);
  AMANZI_ASSERT(n_col_leaves == shape.second);

  // check at every row and column has at least one non-zero entry
  std::vector<bool> cols_ok(n_col_leaves, false);
  for (std::size_t lcv_row = 0; lcv_row != n_row_leaves; ++lcv_row) {
    bool row_ok = false;
    for (std::size_t lcv_col = 0 ; lcv_col != n_col_leaves; ++lcv_col) {
      if (leaves_[lcv_row][lcv_col] != Teuchos::null) {
        row_ok = true;
        cols_ok[lcv_col] = true;
      }
    }
    if (!row_ok) {
      Errors::Message msg("TreeOperator::SymbolicAssemble: no nonempty block in row ");
      msg << lcv_row;
      Exceptions::amanzi_throw(msg);
    }
  }
  for (std::size_t lcv_col = 0; lcv_col != n_col_leaves; ++lcv_col) {
    if (!cols_ok[lcv_col]) {
      Errors::Message msg("TreeOperator::SymbolicAssemble: no nonempty block in column ");
      msg << lcv_col;
      Exceptions::amanzi_throw(msg);
    }
  }

  // assemble the graph structure
  for (std::size_t lcv_row = 0; lcv_row != n_row_leaves; ++lcv_row) {
    for (std::size_t lcv_col = 0 ; lcv_col != n_col_leaves; ++lcv_col) {
      Teuchos::RCP<const Operator> leaf = leaves_[lcv_row][lcv_col];
      if (leaf != Teuchos::null) {
        // NOTE: currently the operator interface keeps this from working
        // on non-square matrices... need to pass both row and col supermaps.
        leaf->SymbolicAssembleMatrix(*get_row_supermap(),
                                      *graph, lcv_row, lcv_col);
      }
    }
  }

  // assemble the graph
  int ierr = graph->FillComplete(get_row_supermap()->Map(),
                                 get_col_supermap()->Map());
  AMANZI_ASSERT(!ierr);

  // create the matrix
  Amat_ = Teuchos::rcp(new MatrixFE(graph));
  A_ = Amat_->Matrix();
}


/* ******************************************************************
* Assemble global matrix from elemental matrices of block operators.
****************************************************************** */
void TreeOperator::AssembleMatrix() {
  AMANZI_ASSERT(leaves_.size() != 0);
  computed_ = false;
  Amat_->Zero();

  // check that each row has at least one non-null operator block
  std::size_t n_row_leaves = leaves_.size();
  std::size_t n_col_leaves = leaves_[0].size();
  for (std::size_t lcv_row = 0; lcv_row != n_row_leaves; ++lcv_row) {
    for (std::size_t lcv_col = 0; lcv_col != n_col_leaves; ++lcv_col) {
      Teuchos::RCP<const Operator> leaf = leaves_[lcv_row][lcv_col];
      if (leaf != Teuchos::null) {
        // NOTE: currently the operator interface keeps this from working
        // on non-square matrices... need to pass both row and col supermaps.
        leaf->AssembleMatrix(*get_row_supermap(),
                              *Amat_, lcv_row, lcv_col);
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


void TreeOperator::InitializeInverse()
{
  // if this is a leaf, just initialize the Operator
  if (get_operator() != Teuchos::null) {
    get_operator()->InitializeInverse();
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


void TreeOperator::ComputeInverse()
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




std::pair<int,int>
Impl::collectTreeOperatorLeaves(TreeOperator& tm, std::vector<std::vector<Teuchos::RCP<Operator>>>& leaves,
                                std::size_t i, std::size_t j)
{
  if (tm.get_operator() != Teuchos::null) {
    leaves[i][j] = tm.get_operator();
    return std::make_pair<int,int>(1,1);
  } else {
    std::size_t i0 = i;
    std::size_t j0 = j;

    int nj = -1;
    for (std::size_t lcv_i = 0; lcv_i != tm.get_row_size(); ++lcv_i) {
      std::size_t ni = -1;
      for (std::size_t lcv_j = 0; lcv_j != tm.get_col_size(); ++lcv_j) {
        std::pair<int,int> delta;
        if (tm.get_block(lcv_i, lcv_j) != Teuchos::null) {
          delta = collectTreeOperatorLeaves(*tm.get_block(lcv_i,lcv_j), leaves, i, j);
        } else {
          delta = std::make_pair<int,int>(getNumTreeVectorLeaves(*tm.get_row_map()->SubVector(lcv_i)),
                                          getNumTreeVectorLeaves(*tm.get_col_map()->SubVector(lcv_j)));
        }
        if (ni == -1) {
          ni = delta.first;
        } else {
          AMANZI_ASSERT(delta.first == ni);
        }
        j += delta.second;
      }

      if (nj == -1) {
        nj = j;
      } else {
        AMANZI_ASSERT(nj == j);
      }
      j = j0;
      i += ni;
    }
    return std::make_pair<int,int>(i - i0, nj - j0);
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
        auto op = block->get_operator();
        if (op == Teuchos::null) {
           msg << "TreeOperator ";
        } else {
          for (auto it = op->begin(); it != op->end(); ++it) {
            msg << "<" << (*it)->schema_string << "> ";
          }
        }
        msg << "\n";
      }
    }
  }
  return msg.str();
}

}  // namespace Operators
}  // namespace Amanzi
