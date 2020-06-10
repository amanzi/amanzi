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

#include "dbc.hh"
#include "GraphFE.hh"
#include "AmanziMatrix.hh"

namespace Amanzi {

// Constructor
GraphFE::GraphFE(const Map_ptr_type& row_map,
		 const Map_ptr_type& ghosted_row_map,
		 const Map_ptr_type& col_map,
		 std::size_t max_nnz_per_row) :
    row_map_(row_map),
    ghosted_row_map_(ghosted_row_map),
    col_map_(col_map) {

  // defaults to square matrix with elemental access patterns
  if (col_map_ == Teuchos::null) col_map_ = ghosted_row_map_;

  // create the offproc maps
  n_used_ = ghosted_row_map_->getNodeNumElements();
  n_owned_ = row_map_->getNodeNumElements();

  // offproc graph is not empty when at least one processor has a ghost cell
  int tmp2, tmp1(n_used_ - n_owned_);
  Teuchos::reduceAll(*ghosted_row_map_->getComm(), Teuchos::REDUCE_MAX, 1, &tmp1, &tmp2);
  includes_ghosted_ = (tmp2 > 0);

  // create the graphs
  graph_ = Teuchos::rcp(new Graph_type(row_map_, col_map_, max_nnz_per_row));

  // potentially create the offproc graphs
  if (includes_ghosted_) {
    std::vector<GO> offproc_gids(n_used_ - n_owned_);
    for (int i=n_owned_; i!=n_used_; ++i)
      offproc_gids[i-n_owned_] = ghosted_row_map_->getGlobalElement(i);
    offproc_row_map_ = Teuchos::rcp(new Map_type(Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid(), 
            offproc_gids.data(), n_used_-n_owned_, 0, ghosted_row_map_->getComm()));
  
    offproc_graph_ = Teuchos::rcp(new Graph_type(offproc_row_map_, col_map_, max_nnz_per_row));
    // create the exporter from offproc to onproc
    exporter_ = Teuchos::rcp(new Export_type(offproc_row_map_, row_map_));
  }
}

GraphFE::GraphFE(const Map_ptr_type& row_map,
		 const Map_ptr_type& ghosted_row_map,
		 const Map_ptr_type& col_map,
		 const Teuchos::ArrayRCP<const std::size_t>& max_nnz_per_row) :
    row_map_(row_map),
    ghosted_row_map_(ghosted_row_map),
    col_map_(col_map) {

  // defaults to square matrix with elemental access patterns
  if (col_map_ == Teuchos::null) col_map_ = ghosted_row_map_;

  // create the offproc maps
  n_used_ = ghosted_row_map_->getNodeNumElements();
  n_owned_ = row_map_->getNodeNumElements();

  // offproc graph is not empty when at least one processor has a ghost cell
  int tmp2, tmp1(n_used_ - n_owned_);
  Teuchos::reduceAll(*ghosted_row_map_->getComm(), Teuchos::REDUCE_MAX, 1, &tmp1, &tmp2);
  includes_ghosted_ = (tmp2 > 0);

  // create the graphs
  graph_ = Teuchos::rcp(new Graph_type(row_map_, col_map_, max_nnz_per_row.view(0, n_owned_)));

  if (includes_ghosted_) {
    std::vector<GO> offproc_gids(n_used_ - n_owned_);
    for (int i=n_owned_; i!=n_used_; ++i)
      offproc_gids[i-n_owned_] = ghosted_row_map_->getGlobalElement(i);

    offproc_row_map_ = Teuchos::rcp(new Map_type(Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid(), 
            offproc_gids.data(), n_used_-n_owned_, 0, ghosted_row_map_->getComm()));
    offproc_graph_ = Teuchos::rcp(new Graph_type(offproc_row_map_, col_map_, max_nnz_per_row.view(n_owned_, n_used_-n_owned_)));

    // create the exporter from offproc to onproc
    exporter_ = Teuchos::rcp(new Export_type(offproc_row_map_, row_map_));
  }
}


// fill graph using local indices
void
GraphFE::insertLocalIndices(LO row, std::size_t count, LO *indices) {
  if (row < n_owned_) {
    graph_->insertLocalIndices(row, count, indices);
  } else {
    offproc_graph_->insertLocalIndices(row-n_owned_, count, indices);
  }
}

// fill graph using global indices
void
GraphFE::insertGlobalIndices(GO row, std::size_t count, GO *indices) {
  LO local_row = ghosted_row_map_->getLocalElement(row);
  AMANZI_ASSERT(local_row >= 0);
  if (local_row < n_owned_) {
    graph_->insertGlobalIndices(row, count, indices);
  } else {
    offproc_graph_->insertGlobalIndices(row, count, indices);
  }
}

// fill graph using local indices
void
GraphFE::insertLocalIndices(std::size_t row_count, LO *row_inds, std::size_t count, LO *indices) {
  for (std::size_t n=0; n!=row_count; ++n)
    insertLocalIndices(row_inds[n], count, indices);
}

// fill graph using global indices
void
GraphFE::insertGlobalIndices(std::size_t row_count, GO *row_inds, std::size_t count, GO *indices) {
  for (std::size_t n=0; n!=row_count; ++n)
    insertGlobalIndices(row_inds[n], count, indices);
}

// start fill
void
GraphFE::resumeFill() {
  offproc_graph_->resumeFill();
  graph_->resumeFill();
}


// finish fill
void
GraphFE::fillComplete(const Map_ptr_type& domain_map,
		      const Map_ptr_type& range_map) {
  domain_map_ = domain_map;
  range_map_ = range_map;

  if (includes_ghosted_) {
    // fill complete the offproc graph
    offproc_graph_->fillComplete(offproc_row_map_, range_map_);

    // scatter offproc into onproc
    graph_->doExport(*offproc_graph_, *exporter_, Tpetra::INSERT);
  }

  // fillcomplete the final graph
  graph_->fillComplete(domain_map_, range_map_);
}
  
} // namespace Amanzi
