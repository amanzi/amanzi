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

#ifndef AMANZI_GRAPH_FE_HH_
#define AMANZI_GRAPH_FE_HH_

#include "Teuchos_RCP.hpp"

// forward declarations
class Epetra_Map;
class Epetra_CrsGraph;

namespace Amanzi {
namespace Operators {

class GraphFE {

 public:

  // Constructor
  GraphFE(const Teuchos::RCP<const Epetra_Map>& row_map,
	  const Teuchos::RCP<const Epetra_Map>& ghosted_row_map,
	  const Teuchos::RCP<const Epetra_Map>& col_map=Teuchos::null,
	  int max_nnz_per_row);

  // Constructor with nnz -- note this should include ghosted rows.
  GraphFE(const Teuchos::RCP<const Epetra_Map>& row_map,
	  const Teuchos::RCP<const Epetra_Map>& ghosted_row_map,
	  const Teuchos::RCP<const Epetra_Map>& col_map=Teuchos::null,
	  const int* max_nnz_per_row);

  // accessors to maps
  const Epetra_Map& DomainMap() const {
    return *domain_map_; }

  const Epetra_Map& RangeMap() const {
    return *range_map_; }

  const Epetra_Map& RowMap() const {
    return *row_map_; }

  const Epetra_Map& GhostedRowMap() const {
    return *ghosted_row_map_; }

  const Epetra_Map& ColMap() const {
    return *col_map_; }

  // accessor to the importer
  const Epetra_Import& Importer() const {
    return *importer_;
  }

  // accessor to graphs
  const Epetra_CrsGraph& Graph() const {
    return *graph_; }

  const Epetra_CrsGraph& OffProcGraph() const {
    return *offproc_graph_; }

  // fill graph
  int InsertMyIndices(int row, int count, int *indices);

  // finish fill
  int FillComplete(const Epetra_Map& domain_map,
		   const Epetra_Map& range_map);

 protected:

  Teuchos::RCP<const Epetra_Map> row_map_;
  Teuchos::RCP<const Epetra_Map> ghosted_row_map_;
  Teuchos::RCP<const Epetra_Map> col_map_;
  Teuchos::RCP<const Epetra_Map> domain_map_;
  Teuchos::RCP<const Epetra_Map> range_map_;

  Teuchos::RCP<Epetra_Map> offproc_row_map_;
  Teuchos::RCP<Epetra_CrsGraph> graph_;
  Teuchos::RCP<Epetra_CrsGraph> offproc_graph_;  

  Teuchos::RCP<Epetra_Import> importer_;

  int n_owned_;
  int n_used_;
  
};


}  // namespace Operators
}  // namespace Amanzi

#endif
