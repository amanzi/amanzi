/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
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

// forward declarations
class Epetra_Map;
class Epetra_CrsGraph;
class Epetra_Export;

namespace Amanzi {
namespace Operators {

class GraphFE {
 public:
  // Constructor
  GraphFE(const Teuchos::RCP<const Epetra_Map>& row_map,
          const Teuchos::RCP<const Epetra_Map>& ghosted_row_map,
          const Teuchos::RCP<const Epetra_Map>& col_map,
          int max_nnz_per_row);

  // Constructor with nnz -- note this should include ghosted rows.
  GraphFE(const Teuchos::RCP<const Epetra_Map>& row_map,
          const Teuchos::RCP<const Epetra_Map>& ghosted_row_map,
          const Teuchos::RCP<const Epetra_Map>& col_map,
          const int* max_nnz_per_row);

  // does this graph include off-process entries?
  bool includes_offproc() const { return includes_ghosted_; }

  // accessors to maps
  const Epetra_Map& DomainMap() const { return *domain_map_; }
  const Epetra_Map& RangeMap() const { return *range_map_; }

  const Epetra_Map& RowMap() const { return *row_map_; }
  const Epetra_Map& ColMap() const { return *col_map_; }

  const Epetra_Map& GhostedRowMap() const { return *ghosted_row_map_; }

  // accessor to the importer
  const Epetra_Export& Exporter() const { return *exporter_; }

  // accessor to graphs
  const Epetra_CrsGraph& Graph() const { return *graph_; }
  const Epetra_CrsGraph& OffProcGraph() const { return *offproc_graph_; }

  // fill graph
  int InsertMyIndices(int row, int count, int* indices);
  int InsertGlobalIndices(int row, int count, int* indices);
  int InsertMyIndices(int row_count, int* row_indices, int col_count, int* col_indices);
  int InsertGlobalIndices(int row_count, int* row_indices, int col_count, int* col_indices);

  // finish fill
  int FillComplete(const Teuchos::RCP<const Epetra_Map>& domain_map,
                   const Teuchos::RCP<const Epetra_Map>& range_map);

 protected:
  Teuchos::RCP<const Epetra_Map> row_map_;
  Teuchos::RCP<const Epetra_Map> ghosted_row_map_;
  Teuchos::RCP<const Epetra_Map> col_map_;
  Teuchos::RCP<const Epetra_Map> domain_map_;
  Teuchos::RCP<const Epetra_Map> range_map_;

  Teuchos::RCP<Epetra_Map> offproc_row_map_;
  Teuchos::RCP<Epetra_CrsGraph> graph_;
  Teuchos::RCP<Epetra_CrsGraph> offproc_graph_;

  Teuchos::RCP<Epetra_Export> exporter_;

  int n_owned_;
  int n_used_;
  bool includes_ghosted_;
};


} // namespace Operators
} // namespace Amanzi

#endif
