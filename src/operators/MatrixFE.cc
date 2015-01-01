/*
Author: Ethan Coon (ecoon@lanl.gov)

A plausibly scalable matrix for use in FE-like systems, where assembly
must be done into rows of ghost entities as well as owned entities.

This matrix uses the "construct, insert, complete fill" paradigm of all
Epetra matrices.  The only real difference is the use of
InserMyIndices(), which may now take local indices from the GHOSTED
map, not the true row map.

*/

#include <vector>
#include "Epetra_Map.h"
#include "Epetra_CrsGraph.h"

#include "MatrixFE.hh"

namespace Amanzi {
namespace Operators {

// Constructor
MatrixFE::MatrixFE(const Teuchos::RCP<const GraphFE>& graph) :
    graph_(graph) {

  // create the crs matrices
  n_used = graph->GhostedRowMap().MyLength();
  n_owned_ = graph->RowMap().MyLength();

  matrix_ = Teuchos::rcp(new Epetra_CrsMatrix(Copy, graph->Graph()));
  offproc_matrix_ = Teuchos::rcp(new Epetra_CrsMatrix(Copy, graph->OffProcGraph()));
}

// fill graph
int
MatrixFE::SumIntoMyValues(int row, int count, double *values, int *indices) {
  int ierr(0);

  if (row < n_owned_) {
    ierr = graph_->SumIntoMyValues(row, count, values, indices);
  } else {
    ierr = offproc_graph_->SumIntoMyValues(row-n_owned_, count, values, indices);
  }

  return ierr;
}

// finish fill
int
MatrixFE::FillComplete() {
  // fill complete the offproc matrix
  int ierr = offproc_matrix_->FillComplete(graph->DomainMap(), graph->RangeMap());
  ASSERT(!ierr);

  // scatter offproc into onproc
  ierr |= matrix_->Export(*offproc_matrix_, graph->Exporter(), Add);
  ASSERT(!ierr);

  // fillcomplete the final graph
  ierr |= matrix_->FillComplete(graph->DomainMap(), graph->RangeMap());
  ASSERT(!ierr);

  return ierr;  
}
  

  
} // namespace Operators
} // namespace Amanzi
