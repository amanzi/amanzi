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
  Map_ptr_type DomainMap() const { return graph_->DomainMap(); }
  Map_ptr_type RangeMap() const { return graph_->RangeMap(); }

  Map_ptr_type RowMap() const { return graph_->RowMap(); }
  Map_ptr_type ColMap() const { return graph_->ColMap(); }

  Map_ptr_type GhostedRowMap() const { return graph_->GhostedRowMap(); }

  // accessor to graphs
  const GraphFE& Graph() const { return *graph_; }

  // accessor to matrices
  cMatrix_ptr_type Matrix() const { return matrix_; }
  Matrix_ptr_type Matrix() { return matrix_; }

  cMatrix_ptr_type OffProcMatrix() const { return offproc_matrix_; }
  Matrix_ptr_type OffProcMatrix() { return offproc_matrix_; }

  // zero to allow mation
  int Zero();

  // fill graph
  int
  SumIntoMyValues(int row, int count, const double* values, const int* indices);

  // local matrix sum
  int SumIntoMyValues(const int *row_inds, const int *col_inds,
                      const Teuchos::SerialDenseMatrix<int,double>& vals);
  int SumIntoMyValues_Transposed(const int *row_inds, const int *col_inds,
          const Teuchos::SerialDenseMatrix<int,double>& vals);

  // local matrix sum
  int SumIntoMyValues(const int *row_inds, const int *col_inds,
                      const WhetStone::DenseMatrix& vals);
  int SumIntoMyValues_Transposed(const int *row_inds, const int *col_inds,
          const WhetStone::DenseMatrix& vals);
  
  // hack the diagonal
  int DiagonalShift(double shift);

  // Passthroughs.
  // --
  // NOTE that currently many of these cannot work on an offproc -- the Export
  // is done with Sum, not Insert. Some could be changed by adding a sum/insert
  // flag, and only having issues when you try to both Sum and Insert without
  // FillComplete called between, but I don't see a need for this
  // functionality.  If you do, ask. --etc
  int InsertMyValues(int row, int count, const double *values, const int *indices);
  int ExtractMyRowCopy(int row, int& count, const double* &values, const int* &indices) const;
  
  // finish fill
  int ResumeFill();
  int FillComplete();

 protected:
  Teuchos::RCP<const GraphFE> graph_;
  Matrix_ptr_type matrix_;
  Matrix_ptr_type offproc_matrix_;  

  int n_owned_;
  int n_used_;
};

} // namespace Operators
} // namespace Amanzi

#endif
