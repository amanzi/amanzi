/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
A plausibly scalable matrix for use in FE-like systems, where assembly
must be done into rows of ghost entities as well as owned entities.

This matrix uses the "construct, insert, complete fill" paradigm of all
Epetra matrices.  The only real difference is the use of
InserMyIndices(), which may now take local indices from the GHOSTED
map, not the true row map.

*/

#ifndef AMANZI_MATRIX_FE_HH_
#define AMANZI_MATRIX_FE_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"

class Epetra_CrsMatrix;
class Epetra_SerialDenseMatrix;

#include "GraphFE.hh"

namespace Amanzi {

namespace WhetStone {
class DenseMatrix;
}

namespace Operators {

class MatrixFE {
 public:
  // Constructor
  MatrixFE(const Teuchos::RCP<const GraphFE>& graph);

  // accessors to maps
  const Epetra_Map& DomainMap() const { return graph_->DomainMap(); }
  const Epetra_Map& RangeMap() const { return graph_->RangeMap(); }

  const Epetra_Map& RowMap() const { return graph_->RowMap(); }
  const Epetra_Map& ColMap() const { return graph_->ColMap(); }

  const Epetra_Map& GhostedRowMap() const { return graph_->GhostedRowMap(); }

  // accessor to graphs
  const GraphFE& Graph() const { return *graph_; }

  // accessor to matrices
  Teuchos::RCP<const Epetra_CrsMatrix> Matrix() const { return matrix_; }
  Teuchos::RCP<Epetra_CrsMatrix> Matrix() { return matrix_; }

  Teuchos::RCP<const Epetra_CrsMatrix> OffProcMatrix() const { return offproc_matrix_; }
  Teuchos::RCP<Epetra_CrsMatrix> OffProcMatrix() { return offproc_matrix_; }

  // zero to allow mation
  int Zero();

  // fill graph
  int SumIntoMyValues(int row, int count, const double* values, const int* indices);

  int SumIntoMyValues(const int* indices, const Epetra_SerialDenseMatrix& vals)
  {
    return SumIntoMyValues(indices, indices, vals);
  }
  int SumIntoMyValues(const int* row_indices,
                      const int* col_indices,
                      const Epetra_SerialDenseMatrix& vals);
  int SumIntoMyValues_Transposed(const int* indices, const Epetra_SerialDenseMatrix& vals)
  {
    return SumIntoMyValues_Transposed(indices, indices, vals);
  }
  int SumIntoMyValues_Transposed(const int* row_indices,
                                 const int* col_indices,
                                 const Epetra_SerialDenseMatrix& vals);

  int SumIntoMyValues(const int* indices, const Teuchos::SerialDenseMatrix<int, double>& vals)
  {
    return SumIntoMyValues(indices, indices, vals);
  }
  int SumIntoMyValues(const int* row_inds,
                      const int* col_inds,
                      const Teuchos::SerialDenseMatrix<int, double>& vals);
  int SumIntoMyValues_Transposed(const int* indices,
                                 const Teuchos::SerialDenseMatrix<int, double>& vals)
  {
    return SumIntoMyValues_Transposed(indices, indices, vals);
  }
  int SumIntoMyValues_Transposed(const int* row_inds,
                                 const int* col_inds,
                                 const Teuchos::SerialDenseMatrix<int, double>& vals);

  int SumIntoMyValues(const int* indices, const WhetStone::DenseMatrix& vals)
  {
    return SumIntoMyValues(indices, indices, vals);
  }
  int SumIntoMyValues(const int* row_inds, const int* col_inds, const WhetStone::DenseMatrix& vals);
  int SumIntoMyValues_Transposed(const int* indices, const WhetStone::DenseMatrix& vals)
  {
    return SumIntoMyValues_Transposed(indices, indices, vals);
  }
  int SumIntoMyValues_Transposed(const int* row_inds,
                                 const int* col_inds,
                                 const WhetStone::DenseMatrix& vals);

  // hack the diagonal
  int DiagonalShift(double shift);

  int DiagonalShiftMin(double shift_min);
  // Passthroughs.
  // --
  // NOTE that currently many of these cannot work on an offproc -- the Export is
  // done with Sum, not Insert. Some could be changed by adding a sum/insert
  // flag, and only having issues when you try to both Sum and Insert without
  // FillComplete called between, but I don't see a need for this
  // functionality.  If you do, ask. --etc
  int InsertMyValues(int row, int count, const double* values, const int* indices);
  int ExtractMyRowCopy(int row, int size, int& count, double* values, int* indices) const;

  // finish fill
  int FillComplete();

 protected:
  Teuchos::RCP<const GraphFE> graph_;
  Teuchos::RCP<Epetra_CrsMatrix> matrix_;
  Teuchos::RCP<Epetra_CrsMatrix> offproc_matrix_;

  int n_owned_;
  int n_used_;
};

} // namespace Operators
} // namespace Amanzi

#endif
