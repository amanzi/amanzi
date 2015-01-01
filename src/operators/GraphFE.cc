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
#include "Epetra_Map.h"
#include "Epetra_CrsGraph.h"

#include "GraphFE.hh"

namespace Amanzi {
namespace Operators {

// Constructor
GraphFE::GraphFE(const Teuchos::RCP<const Epetra_Map>& row_map,
		 const Teuchos::RCP<const Epetra_Map>& ghosted_row_map,
		 const Teuchos::RCP<const Epetra_Map>& col_map,
		 int max_nnz_per_row) :
    row_map_(row_map),
    ghosted_row_map_(ghosted_row_map),
    col_map_(col_map) {

  // defaults to square matrix with elemental access patterns
  if (col_map_ == Teuchos::null) col_map_ = ghosted_row_map_;

  // create the offproc maps
  n_used = ghosted_row_map_->MyLength();
  n_owned_ = row_map_->MyLength();
  
  std::vector<int> offproc_gids(n_used - n_owned_);
  for (int i=n_owned_; i!=n_used; ++i)
    offproc_gids[i-n_owned_] = ghosted_row_map_->GID(i);
  offproc_row_map_ = Teuchos::rcp(new Epetra_Map(-1, n_used-n_owned_,
			&offproc_gids[0], 0, ghosted_row_map_->Comm()));
  
  // create the graphs
  graph_ = Teuchos::rcp(new Epetra_CrsGraph(Copy, *row_map_, *col_map_,
					    max_nnz_per_row, false));
  offproc_graph_ = Teuchos::rcp(new Epetra_CrsGraph(Copy, *offproc_row_map_,
			*col_map_, max_nnz_per_row, false));

  // create the exporter from offproc to onproc
  exporter_ = Teuchos::rcp(new Epetra_Export(*offproc_row_map_, *row_map_));
}

GraphFE::GraphFE(const Teuchos::RCP<const Epetra_Map>& row_map,
		 const Teuchos::RCP<const Epetra_Map>& ghosted_row_map,
		 const Teuchos::RCP<const Epetra_Map>& col_map,
		 const int* max_nnz_per_row) :
    row_map_(row_map),
    ghosted_row_map_(ghosted_row_map),
    col_map_(col_map) {

  // defaults to square matrix with elemental access patterns
  if (col_map_ == Teuchos::null) col_map_ = ghosted_row_map_;

  // create the offproc maps
  n_used = ghosted_row_map_->MyLength();
  n_owned_ = row_map_->MyLength();
  
  std::vector<int> offproc_gids(n_used - n_owned_);
  for (int i=n_owned_; i!=n_used; ++i)
    offproc_gids[i-n_owned_] = ghosted_row_map_->GID(i);
  offproc_row_map_ = Teuchos::rcp(new Epetra_Map(-1, n_used-n_owned_,
			&offproc_gids[0], 0, ghosted_row_map_->Comm()));
  
  // create the graphs
  graph_ = Teuchos::rcp(new Epetra_CrsGraph(Copy, *row_map_, *col_map_,
					    max_nnz_per_row, false));
  offproc_graph_ = Teuchos::rcp(new Epetra_CrsGraph(Copy, *offproc_row_map_,
			*col_map_, max_nnz_per_row+n_owned_, false));

  // create the exporter from offproc to onproc
  exporter_ = Teuchos::rcp(new Epetra_Export(*offproc_row_map_, *row_map_));

}


// fill graph
int
GraphFE::InsertMyIndices(int row, int count, int *indices) {
  int ierr(0);

  if (row < n_owned_) {
    ierr = graph_->InsertMyIndices(row, count, indices);
  } else {
    ierr = offproc_graph_->InsertMyIndices(row-n_owned_, count, indices);
  }

  return ierr;
}

// finish fill
int
GraphFE::FillComplete(const Teuchos::RCP<const Epetra_Map>& domain_map,
		      const Teuchos::RCP<const Epetra_Map>& range_map) {
  domain_map_ = domain_map;
  range_map_ = range_map;

  // fill complete the offproc graph
  int ierr = offproc_graph_->FillComplete(*domain_map_, *range_map_);
  ASSERT(!ierr);

  // scatter offproc into onproc
  ierr |= graph_->Export(*offproc_graph_, *exporter_, Insert);
  ASSERT(!ierr);

  // fillcomplete the final graph
  ierr |= graph_->FillComplete(*domain_map_, *range_map_);
  ASSERT(!ierr);

  return ierr;  
}
  

  
} // namespace Operators
} // namespace Amanzi
