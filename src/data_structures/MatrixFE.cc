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

// Constructor
MatrixFE::MatrixFE(const Teuchos::RCP<const GraphFE>& graph) :
    graph_(graph) {

  // create the crs matrices
  n_used_ = graph_->getGhostedRowMap()->getNodeNumElements();
  n_owned_ = graph_->getRowMap()->getNodeNumElements();

  matrix_ = Teuchos::rcp(new Matrix_type(graph_->getGraph()));
  if (graph_->includes_offproc())
    offproc_matrix_ = Teuchos::rcp(new Matrix_type(graph_->getOffProcGraph()));
}

// zero for summation
void
MatrixFE::zero() {
  matrix_->setAllToScalar(0.);
  if (graph_->includes_offproc())
    offproc_matrix_->setAllToScalar(0.);
}

// fill matrix
void
MatrixFE::sumIntoLocalValues(int row, int count, const double *values, const int *indices) {
  if (row < n_owned_) {
    LO nchanged = matrix_->sumIntoLocalValues(row, count, values, indices);
    AMANZI_ASSERT(nchanged == count);
  } else {
    LO nchanged = offproc_matrix_->sumIntoLocalValues(row-n_owned_, count, values, indices);
    AMANZI_ASSERT(nchanged == count);
  }
}



// Teuchos::SerialDenseMatrices are column-major.
void
MatrixFE::sumIntoLocalValuesTransposed(const int *row_indices, const int *col_indices,
        const Teuchos::SerialDenseMatrix<int,double>& vals) {
  for (int i=0; i!=vals.numCols(); ++i)
    sumIntoLocalValues(row_indices[i], vals.numRows(), vals[i], col_indices);
}

// Teuchhos::SerialDenseMatrices are column-major.
void
MatrixFE::sumIntoLocalValues(const int *row_indices, const int *col_indices,
                          const Teuchos::SerialDenseMatrix<int,double>& vals) {
  std::vector<double> row_vals(vals.numCols());
  for (int i=0; i!=vals.numRows(); ++i) {
    for (int j=0; j!=vals.numCols(); ++j) row_vals[j] = vals(i,j);
    sumIntoLocalValues(row_indices[i], vals.numRows(), row_vals.data(), col_indices);
  }
}


// diagonal shift for (near) singular matrices where the constant vector is the null space
void
MatrixFE::diagonalShift(double shift) {
  matrix_->resumeFill();
  LO local_len = graph_->getRowMap()->getNodeNumElements();
  for (LO i=0; i!=local_len; ++i) {
    sumIntoLocalValues(i, 1, &shift, &i);
  }
  matrix_->fillComplete(graph_->getDomainMap(), graph_->getRangeMap());
}


// Passthroughs
void
MatrixFE::insertLocalValues(int row, int count, const double *values, const int *indices) {
  if (row < n_owned_) {
    matrix_->insertLocalValues(row, count, values, indices);
  } else {
    Errors::Message message("MatrixFE does not support offproc Insert semantics.");
    Exceptions::amanzi_throw(message);
  }
}


void
MatrixFE::getLocalRowView(int row, int& count, const double* &values, const int* &indices) const {
  if (row < n_owned_) {
    matrix_->getLocalRowView(row, count, values, indices);
  } else {
    Errors::Message message("MatrixFE does not support offproc Extract semantics.");
    Exceptions::amanzi_throw(message);
  }
}

// start fill
void
MatrixFE::resumeFill() {
  if (matrix_->isFillComplete()) {
    matrix_->resumeFill();
    if (graph_->includes_offproc()) offproc_matrix_->resumeFill();
  }
}


// finish fill
void
MatrixFE::fillComplete() {
  if (graph_->includes_offproc()) {
    // fill complete the offproc matrix
    offproc_matrix_->fillComplete(graph_->getOffProcGraph()->getDomainMap(), graph_->getOffProcGraph()->getRangeMap());

    // scatter offproc into onproc
    matrix_->doExport(*offproc_matrix_, *graph_->getExporter(), Tpetra::ADD);

    // zero the offproc in case of multiple stage assembly and multiple calls to FillComplete()
    offproc_matrix_->resumeFill();
    offproc_matrix_->setAllToScalar(0.);
  }

  // fillcomplete the final matrix
  matrix_->fillComplete(graph_->getDomainMap(), graph_->getRangeMap());
}

  
} // namespace Amanzi
