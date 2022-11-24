/*
Author: Ethan Coon (ecoon@lanl.gov)

A plausibly scalable graph for use in FE-like systems, where assembly
must be done into rows of ghost entities as well as owned entities.
This graph is taken in the construction of Amanzi's MatrixFE.

This graph uses the "construct, insert, complete fill" paradigm of all
Epetra Graphs.  The only real difference is the use of
InserMyIndices(), which may now take local indices from the GHOSTED
map, not the true row map.

*/

#include <vector>

#include "AmanziComm.hh"
#include "Epetra_Map.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_Export.h"

#include "dbc.hh"
#include "GraphFE.hh"

namespace Amanzi {
namespace Operators {

// Constructor
GraphFE::GraphFE(const Teuchos::RCP<const Epetra_Map>& row_map,
                 const Teuchos::RCP<const Epetra_Map>& ghosted_row_map,
                 const Teuchos::RCP<const Epetra_Map>& col_map,
                 int max_nnz_per_row)
  : row_map_(row_map), ghosted_row_map_(ghosted_row_map), col_map_(col_map)
{
  // defaults to square matrix with elemental access patterns
  if (col_map_ == Teuchos::null) col_map_ = ghosted_row_map_;

  // create the offproc maps
  n_used_ = ghosted_row_map_->NumMyElements();
  n_owned_ = row_map_->NumMyElements();

  // offproc graph is not empty when at least one processor has a mesh
  int tmp2, tmp1(n_used_ - n_owned_);
  ghosted_row_map_->Comm().MaxAll(&tmp1, &tmp2, 1);
  includes_ghosted_ = (tmp2 > 0);

  // create the graphs
  graph_ = Teuchos::rcp(new Epetra_CrsGraph(Copy, *row_map_, *col_map_, max_nnz_per_row, false));

  // potentially create the offproc graphs
  if (includes_ghosted_) {
    std::vector<int> offproc_gids(n_used_ - n_owned_);
    for (int i = n_owned_; i != n_used_; ++i) offproc_gids[i - n_owned_] = ghosted_row_map_->GID(i);
    offproc_row_map_ = Teuchos::rcp(
      new Epetra_Map(-1, n_used_ - n_owned_, &offproc_gids[0], 0, ghosted_row_map_->Comm()));

    offproc_graph_ =
      Teuchos::rcp(new Epetra_CrsGraph(Copy, *offproc_row_map_, *col_map_, max_nnz_per_row, false));
    // create the exporter from offproc to onproc
    exporter_ = Teuchos::rcp(new Epetra_Export(*offproc_row_map_, *row_map_));
  }
}

GraphFE::GraphFE(const Teuchos::RCP<const Epetra_Map>& row_map,
                 const Teuchos::RCP<const Epetra_Map>& ghosted_row_map,
                 const Teuchos::RCP<const Epetra_Map>& col_map,
                 const int* max_nnz_per_row)
  : row_map_(row_map), ghosted_row_map_(ghosted_row_map), col_map_(col_map)
{
  // defaults to square matrix with elemental access patterns
  if (col_map_ == Teuchos::null) col_map_ = ghosted_row_map_;

  // create the offproc maps
  n_used_ = ghosted_row_map_->NumMyElements();
  n_owned_ = row_map_->NumMyElements();

  // offproc graph is not empty when at least one processor has a mesh
  int tmp2, tmp1(n_used_ - n_owned_);
  ghosted_row_map_->Comm().MaxAll(&tmp1, &tmp2, 1);
  includes_ghosted_ = (tmp2 > 0);

  // create the graphs
  graph_ = Teuchos::rcp(new Epetra_CrsGraph(Copy, *row_map_, *col_map_, max_nnz_per_row, false));

  if (includes_ghosted_) {
    std::vector<int> offproc_gids(n_used_ - n_owned_);
    for (int i = n_owned_; i != n_used_; ++i) offproc_gids[i - n_owned_] = ghosted_row_map_->GID(i);
    offproc_row_map_ = Teuchos::rcp(
      new Epetra_Map(-1, n_used_ - n_owned_, &offproc_gids[0], 0, ghosted_row_map_->Comm()));

    offproc_graph_ = Teuchos::rcp(
      new Epetra_CrsGraph(Copy, *offproc_row_map_, *col_map_, max_nnz_per_row + n_owned_, false));

    // create the exporter from offproc to onproc
    exporter_ = Teuchos::rcp(new Epetra_Export(*offproc_row_map_, *row_map_));
  }
}


// fill graph using local indices
int
GraphFE::InsertMyIndices(int row, int count, int* indices)
{
  int ierr(0);

  if (row < n_owned_) {
    std::vector<int> global_indices(count);
    for (int n = 0; n != count; ++n) global_indices[n] = col_map_->GID(indices[n]);
    ierr = graph_->InsertGlobalIndices(row_map_->GID(row), count, &global_indices[0]);
    EPETRA_CHK_ERR(ierr);
  } else {
    ierr = offproc_graph_->InsertMyIndices(row - n_owned_, count, indices);
    EPETRA_CHK_ERR(ierr);
  }

  return ierr;
}

// fill graph using global indices
int
GraphFE::InsertGlobalIndices(int row, int count, int* indices)
{
  int ierr(0);

  int local_row = ghosted_row_map_->LID(row);
  AMANZI_ASSERT(local_row >= 0);
  if (local_row < n_owned_) {
    ierr = graph_->InsertGlobalIndices(row, count, indices);
    EPETRA_CHK_ERR(ierr);
  } else {
    // std::vector<int> local_indices(count);
    // for (int n=0; n!=count; ++n) local_indices[n] = col_map_->LID(indices[n]);
    // ierr = offproc_graph_->InsertMyIndices(local_row-n_owned_, count, &local_indices[0]);
    ierr = offproc_graph_->InsertGlobalIndices(row, count, indices);
    EPETRA_CHK_ERR(ierr);
  }

  return ierr;
}

// fill graph using local indices
int
GraphFE::InsertMyIndices(int row_count, int* row_inds, int count, int* indices)
{
  int ierr(0);

  for (int n = 0; n != row_count; ++n) ierr |= InsertMyIndices(row_inds[n], count, indices);
  return ierr;
}

// fill graph using global indices
int
GraphFE::InsertGlobalIndices(int row_count, int* row_inds, int count, int* indices)
{
  int ierr(0);

  for (int n = 0; n != row_count; ++n) ierr |= InsertGlobalIndices(row_inds[n], count, indices);
  return ierr;
}

// finish fill
int
GraphFE::FillComplete(const Teuchos::RCP<const Epetra_Map>& domain_map,
                      const Teuchos::RCP<const Epetra_Map>& range_map)
{
  int ierr = 0;
  domain_map_ = domain_map;
  range_map_ = range_map;

  if (includes_ghosted_) {
    // fill complete the offproc graph
    ierr = offproc_graph_->FillComplete(*offproc_row_map_, *range_map_);
    EPETRA_CHK_ERR(ierr);
    AMANZI_ASSERT(!ierr);

    // scatter offproc into onproc
    ierr |= graph_->Export(*offproc_graph_, *exporter_, Insert);
    EPETRA_CHK_ERR(ierr);
    AMANZI_ASSERT(!ierr);
  }

  // fillcomplete the final graph
  ierr |= graph_->FillComplete(*domain_map_, *range_map_);
  EPETRA_CHK_ERR(ierr);
  AMANZI_ASSERT(!ierr);

  return ierr;
}

} // namespace Operators
} // namespace Amanzi
