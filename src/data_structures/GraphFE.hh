/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:

*/

//!

/*
Author: Ethan Coon (coonet@ornl.gov)

A plausibly scalable graph for use in FE-like systems, where assembly
must be done into rows of ghost entities as well as owned entities.
This graph is taken in the construction of Amanzi's MatrixFE.

This graph uses the "construct, insert, complete fill" paradigm of all
Epetra Graphs.  The only real difference is the use of
InserMyIndices(), which may now take local indices from the GHOSTED
map, not the true row map.

*/

#ifndef AMANZI_GRAPH_FE_HH_
#define AMANZI_GRAPH_FE_HH_

#include "Teuchos_RCP.hpp"

#include "AmanziTypes.hh"

namespace Amanzi {
namespace Operators {

class GraphFE {
 public:
  // Constructor
  GraphFE(const Map_ptr_type& row_map,
	  const Map_ptr_type& ghosted_row_map,
	  const Map_ptr_type& col_map,
	  std::size_t max_nnz_per_row);

  // Constructor with nnz -- note this should include ghosted rows.
  GraphFE(const Map_ptr_type& row_map,
	  const Map_ptr_type& ghosted_row_map,
	  const Map_ptr_type& col_map,
	  const Teuchos::ArrayRCP<const std::size_t>& max_nnz_per_row);

  // does this graph include off-process entries?
  bool includes_offproc() const { return includes_ghosted_; }

  // accessors to maps
  Map_ptr_type DomainMap() const { return domain_map_; }
  Map_ptr_type RangeMap() const { return range_map_; }

  Map_ptr_type RowMap() const { return row_map_; }
  Map_ptr_type ColMap() const { return col_map_; }

  Map_ptr_type GhostedRowMap() const { return ghosted_row_map_; }

  // accessor to the importer
  const Export_type& Exporter() const { return *exporter_; }

  // accessor to graphs
  Graph_ptr_type Graph() const { return graph_; }
  Graph_ptr_type OffProcGraph() const { return offproc_graph_; }

  // fill graph
  int InsertMyIndices(LO row, std::size_t count, LO *indices);
  int InsertGlobalIndices(GO row, std::size_t count, GO *indices);
  int InsertMyIndices(std::size_t row_count, LO *row_indices, std::size_t col_count, LO *col_indices);
  int InsertGlobalIndices(std::size_t row_count, GO *row_indices, std::size_t col_count, GO *col_indices);

  // finish fill
  int ResumeFill();
  int FillComplete(const Map_ptr_type& domain_map,
                   const Map_ptr_type& range_map);

 protected:
  Map_ptr_type row_map_;
  Map_ptr_type ghosted_row_map_;
  Map_ptr_type col_map_;
  Map_ptr_type domain_map_;
  Map_ptr_type range_map_;

  Map_ptr_type offproc_row_map_;
  Graph_ptr_type graph_;
  Graph_ptr_type offproc_graph_;  

  Export_ptr_type exporter_;

  LO n_owned_;
  LO n_used_;
  bool includes_ghosted_;
};


} // namespace Operators
} // namespace Amanzi

#endif
