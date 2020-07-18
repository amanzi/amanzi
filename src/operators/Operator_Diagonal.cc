/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
           Ethan Coon (ecoon@lanl.gov)

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
int Operator_Diagonal::ApplyMatrixFreeOp(
    const Op_Diagonal& op, const CompositeVector& X, CompositeVector& Y) const
{
X.Print(std::cout, false);
  const Epetra_MultiVector& Xi = *X.ViewComponent(col_compname_, true);
  Epetra_MultiVector& Yi = *Y.ViewComponent(row_compname_, true);
 
  const auto& col_lids = op.col_inds();
  const auto& row_lids = op.row_inds();

  for (int n = 0; n != row_lids.size(); ++n) {
    const WhetStone::DenseMatrix& Acell = op.matrices[n];
    int nrows = Acell.NumRows();
    int ncols = Acell.NumCols();

    WhetStone::DenseVector v(ncols), av(nrows);
    for (int i = 0; i != ncols; ++i) {
      v(i) = Xi[0][col_lids[n][i]];
    }

    Acell.Multiply(v, av, false);

    for (int i = 0; i != nrows; ++i) {
      Yi[0][row_lids[n][i]] += av(i);
std::cout << n << ":  " << row_lids[n][i] << " " << col_lids[n][i] << " a=" << av(i) << std::endl;
    }
  }

  return 0;
}


/* ******************************************************************
* Visit methods for symbolic assemble.
****************************************************************** */
void Operator_Diagonal::SymbolicAssembleMatrixOp(
    const Op_Diagonal& op, const SuperMap& map, GraphFE& graph,
    int my_block_row, int my_block_col) const
{
std::cout << "Component: " << row_compname_ << " " << col_compname_ << std::endl;
  const std::vector<int>& row_gids = map.GhostIndices(my_block_row, row_compname_, 0);
  const std::vector<int>& col_gids = map.GhostIndices(my_block_col, col_compname_, 0);

  const auto& col_lids = op.col_inds();
  const auto& row_lids = op.row_inds();

  std::vector<int> lid_r, lid_c;

  int ierr(0);
  for (int n = 0; n != col_lids.size(); ++n) {
    int ndofs = col_lids[n].size();

    lid_r.clear();
    lid_c.clear();

    for (int i = 0; i != ndofs; ++i) {
      lid_r.push_back(row_gids[row_lids[n][i]]);
      lid_c.push_back(col_gids[col_lids[n][i]]);
std::cout << n << ":  " << row_lids[n][i] << " " << col_lids[n][i]  << std::endl;
    }
    ierr |= graph.InsertMyIndices(ndofs, lid_r.data(), ndofs, lid_c.data());
  }
  AMANZI_ASSERT(!ierr);
}


/* ******************************************************************
* Visit methods for assemble
****************************************************************** */
void Operator_Diagonal::AssembleMatrixOp(
    const Op_Diagonal& op, const SuperMap& map, MatrixFE& mat,
    int my_block_row, int my_block_col) const
{
  const std::vector<int>& row_gids = map.GhostIndices(my_block_row, row_compname_, 0);
  const std::vector<int>& col_gids = map.GhostIndices(my_block_col, col_compname_, 0);

  const auto& col_lids = op.col_inds();
  const auto& row_lids = op.row_inds();

  std::vector<int> lid_r, lid_c;

  int ierr(0);
  for (int n = 0; n != col_lids.size(); ++n) {
    int ndofs = col_lids[n].size();

    lid_r.clear();
    lid_c.clear();

    for (int i = 0; i != ndofs; ++i) {
      lid_r.push_back(row_gids[row_lids[n][i]]);
      lid_c.push_back(col_gids[col_lids[n][i]]);
    }

    ierr |= mat.SumIntoMyValues(lid_r.data(), lid_c.data(), op.matrices[n]);
  }
  AMANZI_ASSERT(!ierr);
}

}  // namespace Operators
}  // namespace Amanzi



