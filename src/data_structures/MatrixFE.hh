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

#pragma once

#include "Teuchos_RCP.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"

#include "DenseMatrix.hh"
#include "GraphFE.hh"


namespace Amanzi {

class MatrixFE {
 public:
  // Constructor
  MatrixFE(const Teuchos::RCP<const GraphFE>& graph);

  // accessors to maps
  Map_ptr_type getDomainMap() const { return graph_->getDomainMap(); }
  Map_ptr_type getRangeMap() const { return graph_->getRangeMap(); }

  Map_ptr_type getRowMap() const { return graph_->getRowMap(); }
  Map_ptr_type getColMap() const { return graph_->getColMap(); }

  Map_ptr_type getGhostedRowMap() const { return graph_->getGhostedRowMap(); }

  // accessor to graphs
  Teuchos::RCP<const GraphFE> getGraph() const { return graph_; }

  // accessor to matrices
  cMatrix_ptr_type getMatrix() const { return matrix_; }
  Matrix_ptr_type getMatrix() { return matrix_; }

  cMatrix_ptr_type getOffProcMatrix() const { return offproc_matrix_; }
  Matrix_ptr_type getOffProcMatrix() { return offproc_matrix_; }

  // zero to allow mation
  void zero();

  // fill graph
  void
  sumIntoLocalValues(int row, int count, const double* values, const int* indices);

  // local matrix sum
  void sumIntoLocalValues(const int *row_inds, const int *col_inds,
                      const Teuchos::SerialDenseMatrix<int,double>& vals);
  void sumIntoLocalValuesTransposed(const int *row_inds, const int *col_inds,
          const Teuchos::SerialDenseMatrix<int,double>& vals);

  void sumIntoLocalValues(const int *row_inds, const int *col_inds,
                      const WhetStone::DenseMatrix<>& vals);
  void sumIntoLocalValuesTransposed(const int *row_inds, const int *col_inds,
          const WhetStone::DenseMatrix<>& vals);
  
  // hack the diagonal
  void diagonalShift(double shift);

  // Passthroughs.
  // --
  // NOTE that currently many of these cannot work on an offproc -- the Export
  // is done with Sum, not Insert. Some could be changed by adding a sum/insert
  // flag, and only having issues when you try to both Sum and Insert without
  // FillComplete called between, but I don't see a need for this
  // functionality.  If you do, ask. --etc
  void insertLocalValues(int row, int count, const double *values, const int *indices);
  void getLocalRowView(int row, int& count, const double* &values, const int* &indices) const;
  
  // finish fill
  void resumeFill();
  void fillComplete();

 protected:
  Teuchos::RCP<const GraphFE> graph_;
  Matrix_ptr_type matrix_;
  Matrix_ptr_type offproc_matrix_;  

  int n_owned_;
  int n_used_;
};

} // namespace Amanzi
