/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (coonet@ornl.gov)
*/

//! A graph for matrix structure in Tpetra

/*

A plausibly scalable graph for use in FE-like systems, where assembly must be
done into rows of ghost entities as well as owned entities.  This graph is
taken in the construction of Amanzi's MatrixFE.

This graph uses the "construct, insert, complete fill" paradigm of all *petra
Graphs.  The only real difference is the use of InserMyIndices(), which may now
take local indices from the GHOSTED map, not the true row map.

This greatly improves performance relative to providing global IDs, which is
necessary in vanilla *petra because no GHOSTED map is provided to the matrix
class.  When global IDs are provided, a search/lookup/export must be done each
time -- but thanks to the halo-exchange, we already have computed this export,
so we can just use it.

*/

#pragma once

#include "Teuchos_RCP.hpp"

#include "AmanziTypes.hh"

namespace Amanzi {

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
  Map_ptr_type getDomainMap() const { return domain_map_; }
  Map_ptr_type getRangeMap() const { return range_map_; }
  Map_ptr_type getRowMap() const { return row_map_; }
  Map_ptr_type getColMap() const { return col_map_; }
  Map_ptr_type getGhostedRowMap() const { return ghosted_row_map_; }

  // accessor to the importer
  Export_ptr_type getExporter() const { return exporter_; }

  // accessor to graphs
  Graph_ptr_type getGraph() const { return graph_; }
  Graph_ptr_type getOffProcGraph() const { return offproc_graph_; }

  // fill graph
  void insertLocalIndices(LO row, std::size_t count, LO *indices);
  void insertGlobalIndices(GO row, std::size_t count, GO *indices);
  void insertLocalIndices(std::size_t row_count, LO *row_indices, std::size_t col_count, LO *col_indices);
  void insertGlobalIndices(std::size_t row_count, GO *row_indices, std::size_t col_count, GO *col_indices);

  // finish fill
  void resumeFill();
  void fillComplete(const Map_ptr_type& domain_map,
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


} // namespace Amanzi
