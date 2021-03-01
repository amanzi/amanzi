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
#include "Epetra_Vector.h"
#include "Epetra_SerialDenseMatrix.h"

#include "dbc.hh"
#include "errors.hh"
#include "DenseMatrix.hh"
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
  if (graph_->includes_offproc())
    offproc_matrix_ = Teuchos::rcp(new Epetra_CrsMatrix(Copy, graph_->OffProcGraph()));
}

// zero for summation
int
MatrixFE::Zero() {
  int ierr(0);
  ierr = matrix_->PutScalar(0.);
  if (graph_->includes_offproc())
    ierr |= offproc_matrix_->PutScalar(0.);
  return ierr;
}

// fill matrix
int
MatrixFE::SumIntoMyValues(int row, int count, const double *values, const int *indices) {
  int ierr(0);

  if (row < n_owned_) {
    ierr = matrix_->SumIntoMyValues(row, count, values, indices);
  } else {
    ierr = offproc_matrix_->SumIntoMyValues(row-n_owned_, count, values, indices);
  }

  return ierr;
}


// Epetra_SerialDenseMatrices are column-major.
int
MatrixFE::SumIntoMyValues_Transposed(const int *row_indices, const int *col_indices,
        const Epetra_SerialDenseMatrix& vals) {
  int ierr(0);
  for (int i=0; i!=vals.N(); ++i)
    ierr |= SumIntoMyValues(row_indices[i], vals.M(), vals[i], col_indices);
  return ierr;
}

// Epetra_SerialDenseMatrices are column-major.
int
MatrixFE::SumIntoMyValues(const int *row_indices, const int *col_indices,
                          const Epetra_SerialDenseMatrix& vals) {
  int ierr(0);
  std::vector<double> row_vals(vals.N());
  for (int i=0; i!=vals.M(); ++i) {
    for (int j=0; j!=vals.N(); ++j) row_vals[j] = vals(i,j);
    ierr |= SumIntoMyValues(row_indices[i], vals.N(), &row_vals[0], col_indices);
  }
  return ierr;
}

// Teuchos::SerialDenseMatrices are column-major.
int
MatrixFE::SumIntoMyValues_Transposed(const int *row_indices, const int *col_indices,
        const Teuchos::SerialDenseMatrix<int,double>& vals) {
  int ierr(0);
  for (int i=0; i!=vals.numCols(); ++i)
    ierr |= SumIntoMyValues(row_indices[i], vals.numRows(), vals[i], col_indices);
  return ierr;
}

// Teuchhos::SerialDenseMatrices are column-major.
int
MatrixFE::SumIntoMyValues(const int *row_indices, const int *col_indices,
                          const Teuchos::SerialDenseMatrix<int,double>& vals) {
  int ierr(0);
  std::vector<double> row_vals(vals.numCols());
  for (int i=0; i!=vals.numRows(); ++i) {
    for (int j=0; j!=vals.numCols(); ++j) row_vals[j] = vals(i,j);
    ierr |= SumIntoMyValues(row_indices[i], vals.numRows(), &row_vals[0], col_indices);
  }
  return ierr;
}


// WhetStone::DenseMatrices are column-major.
int
MatrixFE::SumIntoMyValues_Transposed(const int *row_indices, const int *col_indices,
        const WhetStone::DenseMatrix& vals) {
  int ierr(0);
  for (int i=0; i!=vals.NumCols(); ++i)
    ierr |= SumIntoMyValues(row_indices[i], vals.NumRows(), vals.Value(0,i), col_indices);
  return ierr;
}

// WhetStone::DenseMatrix are column-major.
int
MatrixFE::SumIntoMyValues(const int *row_indices, const int *col_indices,
                          const WhetStone::DenseMatrix& vals) {
  int ierr(0);
  std::vector<double> row_vals(vals.NumCols());
  for (int i=0; i!=vals.NumRows(); ++i) {
    for (int j=0; j!=vals.NumCols(); ++j) row_vals[j] = vals(i,j);
    ierr |= SumIntoMyValues(row_indices[i], vals.NumCols(), &row_vals[0], col_indices);
  }
  return ierr;
}

// diagonal shift for (near) singular matrices where the constant vector is the null space
int
MatrixFE::DiagonalShift(double shift) {
  int ierr(0);
  Epetra_Vector diag(RowMap());
  ierr = matrix_->ExtractDiagonalCopy(diag);
  for (int i=0; i!=diag.MyLength(); ++i) {
    shift = diag[i] > 0 ? 1.0e-6 : shift;
    diag[i] += shift;
  }
  ierr |= matrix_->ReplaceDiagonalValues(diag);  
  return ierr;
}


// Passthroughs
int
MatrixFE::InsertMyValues(int row, int count, const double *values, const int *indices) {
  int ierr(0);

  if (row < n_owned_) {
    ierr = matrix_->InsertMyValues(row, count, values, indices);
  } else {
    ierr = -1;
    Errors::Message message("MatrixFE does not support offproc Insert semantics.");
    Exceptions::amanzi_throw(message);
  }
  return ierr;
}


int
MatrixFE::ExtractMyRowCopy(int row, int size, int& count,
                           double *values, int *indices) const {
  int ierr(0);
  if (row < n_owned_) {
    ierr = matrix_->ExtractMyRowCopy(row, size, count, values, indices);
  } else {
    ierr = -1;
    Errors::Message message("MatrixFE does not support offproc Extract semantics.");
    Exceptions::amanzi_throw(message);
  }
  return ierr;
}


// finish fill
int
MatrixFE::FillComplete() {
  int ierr = 0;

  if (graph_->includes_offproc()) {
    // fill complete the offproc matrix
    ierr |= offproc_matrix_->FillComplete(graph_->DomainMap(), graph_->RangeMap());
    AMANZI_ASSERT(!ierr);

    // scatter offproc into onproc
    ierr |= matrix_->Export(*offproc_matrix_, graph_->Exporter(), Add);
    AMANZI_ASSERT(!ierr);

    // zero the offproc in case of multiple stage assembly and multiple calls to FillComplete()
    ierr |= offproc_matrix_->PutScalar(0.);
    AMANZI_ASSERT(!ierr);
  }

  // fillcomplete the final graph
  ierr |= matrix_->FillComplete(graph_->DomainMap(), graph_->RangeMap());
  AMANZI_ASSERT(!ierr);

  return ierr;  
}
  

  
} // namespace Operators
} // namespace Amanzi
