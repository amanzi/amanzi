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
#include "Epetra_CrsMatrix.h"
#include "Epetra_SerialDenseMatrix.h"

#include "dbc.hh"
#include "MatrixFE.hh"

namespace Amanzi {
namespace Operators {

// Constructor
MatrixFE::MatrixFE(const Teuchos::RCP<const GraphFE>& graph) :
    graph_(graph) {

  // create the crs matrices
  n_used_ = graph_->GhostedRowMap().NumMyElements();
  n_owned_ = graph_->RowMap().NumMyElements();

  matrix_ = Teuchos::rcp(new Epetra_CrsMatrix(Copy, graph_->Graph()));
  offproc_matrix_ = Teuchos::rcp(new Epetra_CrsMatrix(Copy, graph_->OffProcGraph()));
}

// fill matrix
int
MatrixFE::SumIntoMyValues(int row, int count, double *values, int *indices) {
  int ierr(0);

  if (row < n_owned_) {
    ierr = matrix_->SumIntoMyValues(row, count, values, indices);
  } else {
    ierr = offproc_matrix_->SumIntoMyValues(row-n_owned_, count, values, indices);
  }

  return ierr;
}


int
MatrixFE::SumIntoMyValues_Transposed(int *indices, Epetra_SerialDenseMatrix& vals) {
  int ierr(0);
  for (int i=0; i!=vals.N(); ++i)
    ierr |= SumIntoMyValues(indices[i], vals.M(), vals[i], indices);
  return ierr;
}

// Epetra_SerialDenseMatrices are column-major!
int
MatrixFE::SumIntoMyValues(int *indices, Epetra_SerialDenseMatrix& vals) {
  int ierr(0);
  std::vector<double> row_vals(vals.N());
  for (int i=0; i!=vals.M(); ++i) {
    for (int j=0; j!=vals.N(); ++j) row_vals[j] = vals(i,j);
    ierr |= SumIntoMyValues(indices[i], vals.N(), &row_vals[0], indices);
  }
  return ierr;
}


// finish fill
int
MatrixFE::FillComplete() {
  // fill complete the offproc matrix
  int ierr = offproc_matrix_->FillComplete(graph_->DomainMap(), graph_->RangeMap());
  ASSERT(!ierr);

  // scatter offproc into onproc
  ierr |= matrix_->Export(*offproc_matrix_, graph_->Exporter(), Add);
  ASSERT(!ierr);

  // fillcomplete the final graph
  ierr |= matrix_->FillComplete(graph_->DomainMap(), graph_->RangeMap());
  ASSERT(!ierr);

  return ierr;  
}
  

  
} // namespace Operators
} // namespace Amanzi
