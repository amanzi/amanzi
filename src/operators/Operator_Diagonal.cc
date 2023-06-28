/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
           Ethan Coon (ecoon@lanl.gov)
*/

/*
  Operators

  Operator whose unknowns are gives by the list of
*/

#include "DenseMatrix.hh"
#include "Op_Diagonal.hh"

#include "SuperMap.hh"
#include "GraphFE.hh"
#include "MatrixFE.hh"

#include "OperatorDefs.hh"
#include "Operator_Diagonal.hh"

namespace Amanzi {
namespace Operators {

/* ******************************************************************
* Apply the local matrices directly as schemas match.
****************************************************************** */
int
Operator_Diagonal::ApplyMatrixFreeOp(const Op_Diagonal& op,
                                     const CompositeVector& X,
                                     CompositeVector& Y) const
{
  const Epetra_MultiVector& Xi = *X.ViewComponent(op.col_compname(), true);
  Epetra_MultiVector& Yi = *Y.ViewComponent(op.row_compname(), true);

  const auto& col_lids = op.col_inds();
  const auto& row_lids = op.row_inds();
  AMANZI_ASSERT(row_lids.size() == op.matrices.size());
  AMANZI_ASSERT(col_lids.size() == op.matrices.size());

  for (int n = 0; n != row_lids.size(); ++n) {
    const WhetStone::DenseMatrix& Acell = op.matrices[n];
    int nrows = Acell.NumRows();
    int ncols = Acell.NumCols();

    AMANZI_ASSERT(col_lids[n].size() == ncols);
    AMANZI_ASSERT(row_lids[n].size() == nrows);

    WhetStone::DenseVector v(ncols), av(nrows);
    for (int i = 0; i != ncols; ++i) {
      int col_lid = col_lids[n][i];
      AMANZI_ASSERT(col_lid >= 0 && col_lid < Xi.MyLength());
      v(i) = Xi[0][col_lid];
    }

    Acell.Multiply(v, av, false);

    for (int i = 0; i != nrows; ++i) {
      int row_lid = row_lids[n][i];
      AMANZI_ASSERT(row_lid >= 0 && row_lid < Yi.MyLength());
      Yi[0][row_lid] += av(i);
    }
  }

  return 0;
}


/* ******************************************************************
* Visit methods for symbolic assemble.
****************************************************************** */
void
Operator_Diagonal::SymbolicAssembleMatrixOp(const Op_Diagonal& op,
                                            const SuperMap& map,
                                            GraphFE& graph,
                                            int my_block_row,
                                            int my_block_col) const
{
  const std::vector<int>& row_gids = map.GhostIndices(my_block_row, op.row_compname(), 0);
  const std::vector<int>& col_gids = map.GhostIndices(my_block_col, op.col_compname(), 0);

  const auto& col_lids = op.col_inds();
  const auto& row_lids = op.row_inds();

  std::vector<int> lid_r, lid_c;

  int ierr(0);
  for (int n = 0; n != col_lids.size(); ++n) {
    int ncols = col_lids[n].size();
    int nrows = row_lids[n].size();

    lid_r.clear();
    lid_c.clear();

    for (int i = 0; i != ncols; ++i) lid_c.push_back(col_gids[col_lids[n][i]]);
    for (int i = 0; i != nrows; ++i) lid_r.push_back(row_gids[row_lids[n][i]]);

    ierr |= graph.InsertMyIndices(nrows, lid_r.data(), ncols, lid_c.data());
  }
  AMANZI_ASSERT(!ierr);
}


/* ******************************************************************
* Visit methods for assemble
****************************************************************** */
void
Operator_Diagonal::AssembleMatrixOp(const Op_Diagonal& op,
                                    const SuperMap& map,
                                    MatrixFE& mat,
                                    int my_block_row,
                                    int my_block_col) const
{
  const std::vector<int>& row_gids = map.GhostIndices(my_block_row, op.row_compname(), 0);
  const std::vector<int>& col_gids = map.GhostIndices(my_block_col, op.col_compname(), 0);

  const auto& col_lids = op.col_inds();
  const auto& row_lids = op.row_inds();

  std::vector<int> lid_r, lid_c;

  int ierr(0);
  for (int n = 0; n != col_lids.size(); ++n) {
    int ncols = col_lids[n].size();
    int nrows = row_lids[n].size();

    lid_r.clear();
    lid_c.clear();

    for (int i = 0; i != ncols; ++i) lid_c.push_back(col_gids[col_lids[n][i]]);
    for (int i = 0; i != nrows; ++i) lid_r.push_back(row_gids[row_lids[n][i]]);

    ierr |= mat.SumIntoMyValues(lid_r.data(), lid_c.data(), op.matrices[n]);
  }
  AMANZI_ASSERT(!ierr);
}


/* ******************************************************************
* Copy constructor.
****************************************************************** */
Teuchos::RCP<Operator>
Operator_Diagonal::Clone() const
{
  return Teuchos::rcp(new Operator_Diagonal(*this));
}

} // namespace Operators
} // namespace Amanzi
