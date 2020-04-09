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
#include "Teuchos_SerialDenseMatrix.hpp"

#include "dbc.hh"
#include "errors.hh"
#include "DenseMatrix.hh"
#include "MatrixFE.hh"
#include "AmanziMatrix.hh"

namespace Amanzi {
namespace Operators {

// Constructor
MatrixFE::MatrixFE(const Teuchos::RCP<const GraphFE>& graph) :
    graph_(graph) {

  // create the crs matrices
  n_used_ = graph_->GhostedRowMap()->getNodeNumElements();
  n_owned_ = graph_->RowMap()->getNodeNumElements();

  matrix_ = Teuchos::rcp(new Matrix_type(graph_->Graph()));
  if (graph_->includes_offproc())
    offproc_matrix_ = Teuchos::rcp(new Matrix_type(graph_->OffProcGraph()));
}

// zero for summation
int
MatrixFE::Zero() {
  int ierr(0);
  matrix_->setAllToScalar(0.);
  if (graph_->includes_offproc())
    offproc_matrix_->setAllToScalar(0.);
  return ierr;
}

// fill matrix
int
MatrixFE::SumIntoMyValues(int row, int count, const double *values, const int *indices) {
  LO ierr(0);
  if (row < n_owned_) {
    LO nchanged = matrix_->sumIntoLocalValues(row, count, values, indices);
    AMANZI_ASSERT(nchanged == count);
  } else {
    LO nchanged = offproc_matrix_->sumIntoLocalValues(row-n_owned_, count, values, indices);
    AMANZI_ASSERT(nchanged == count);
  }
  return (int) ierr;
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
    ierr |= SumIntoMyValues(row_indices[i], vals.numRows(), row_vals.data(), col_indices);
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
    ierr |= SumIntoMyValues(row_indices[i], vals.NumCols(), row_vals.data(), col_indices);
  }
  return ierr;
}

// diagonal shift for (near) singular matrices where the constant vector is the null space
int
MatrixFE::DiagonalShift(double shift) {
  int ierr(0);

  matrix_->resumeFill();
  LO local_len = graph_->RowMap()->getNodeNumElements();
  for (LO i=0; i!=local_len; ++i) {
    SumIntoMyValues(i, 1, &shift, &i);
  }
  matrix_->fillComplete(graph_->DomainMap(), graph_->RangeMap());
  return ierr;
}


// Passthroughs
int
MatrixFE::InsertMyValues(int row, int count, const double *values, const int *indices) {
  int ierr(0);

  if (row < n_owned_) {
    matrix_->insertLocalValues(row, count, values, indices);
  } else {
    ierr = -1;
    Errors::Message message("MatrixFE does not support offproc Insert semantics.");
    Exceptions::amanzi_throw(message);
  }
  return ierr;
}


int
MatrixFE::ExtractMyRowCopy(int row, int& count, const double* &values, const int* &indices) const {
  int ierr(0);
  if (row < n_owned_) {
    matrix_->getLocalRowView(row, count, values, indices);
  } else {
    ierr = -1;
    Errors::Message message("MatrixFE does not support offproc Extract semantics.");
    Exceptions::amanzi_throw(message);
  }
  return ierr;
}

// start fill
int
MatrixFE::ResumeFill() {
  offproc_matrix_->resumeFill();
  matrix_->resumeFill();
  return 0;
}


// finish fill
int
MatrixFE::FillComplete() {
  int ierr = 0;

  if (graph_->includes_offproc()) {
    // fill complete the offproc matrix
    offproc_matrix_->fillComplete(graph_->DomainMap(), graph_->RangeMap());

    // scatter offproc into onproc
    matrix_->doExport(*offproc_matrix_, graph_->Exporter(), Tpetra::ADD);

    // zero the offproc in case of multiple stage assembly and multiple calls to FillComplete()
    offproc_matrix_->resumeFill();
    offproc_matrix_->setAllToScalar(0.);
  }

  // fillcomplete the final matrix
  matrix_->fillComplete(graph_->DomainMap(), graph_->RangeMap());

  return ierr;  
}
  

  
} // namespace Operators
} // namespace Amanzi
