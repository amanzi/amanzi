/*
  Chemistry 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ben Andre
*/

#include <cmath>
#include <sstream>
#include <iomanip>

#include "MatrixBlock.hh"

namespace Amanzi {
namespace AmanziChemistry {

MatrixBlock::MatrixBlock() 
  : size_(0),
    cols_(0),
    A_(NULL) {
}


MatrixBlock::MatrixBlock(int size)
  : size_(size),
    cols_(size) {
  A_ = new double[size_ * size_];
}


MatrixBlock::MatrixBlock(int size, int cols)
  : size_(size),
    cols_(cols) {
  A_ = new double[size_ * cols_];
}


MatrixBlock::~MatrixBlock() {
  if (A_) delete [] A_;
  A_ = NULL;
}


/* ******************************************************************
* Change matrix size
****************************************************************** */
void MatrixBlock::Resize(int new_size) {
  if (A_) delete [] A_;
  size_ = new_size;
  cols_ = new_size;
  A_ = new double[size_ * cols_];
}


void MatrixBlock::Zero() {
  for (int i = 0; i < size_ * cols_; i++) A_[i] = 0.0;
}


void MatrixBlock::SetDiagonal(double d) {
  for (int i = 0; i < size_; i++) {
    A_[i * size_ + i] = d;
  }
}


double MatrixBlock::GetRowAbsMax(int irow) {
  double max = 0.0;
  for (int i = 0; i < cols_; i++) {
    double value = std::fabs(A_[i * size_ + irow]);
    if (value > max) {
      max = value;
    }
  }
  return max;
}


void MatrixBlock::ScaleRow(int irow, double scale) {
  for (int i = 0; i < cols_; i++) {
    A_[i * size_ + irow] *= scale;
  }
}


void MatrixBlock::ScaleColumn(int icol, double scale) {
  for (int i = 0; i < size_; i++) {
    A_[icol * size_ + i] *= scale;
  }
}


void MatrixBlock::Scale(double scale) {
  for (int i = 0; i < size_ * cols_; i++) A_[i] *= scale;
}


void MatrixBlock::AddValue(int i, int j, double value) {
  A_[j * size_ + i] += value;
}


void MatrixBlock::AddValues(MatrixBlock* b) {
  double* B = b->GetValues();
  for (int i = 0; i < size_ * cols_; i++) A_[i] += B[i];
}


void MatrixBlock::AddValues(MatrixBlock* b, double scale) {
  double* B = b->GetValues();
  for (int i = 0; i < size_ * cols_; i++) A_[i] += scale * B[i];
}


/* ******************************************************************
* We assume that matrix is square
****************************************************************** */
int MatrixBlock::Inverse() 
{
  int ierr;
  int iwork[size_];
  DGETRF_F77(&size_, &size_, A_, &size_, iwork, &ierr);
  if (ierr) return ierr;
 
  int lwork = size_ * size_;
  double dwork[lwork];
  DGETRI_F77(&size_, A_, &size_, iwork, dwork, &lwork, &ierr);
  return ierr;
}


/* ******************************************************************
* Compute product A * B.
****************************************************************** */
void Multiply(const MatrixBlock& A, const MatrixBlock& B, MatrixBlock& AB)
{
  const double* dataA = A.GetValues();
  const double* dataB = B.GetValues();

  int nrowsA = A.size(), ncolsA = A.cols();
  int nrowsB = B.size(), ncolsB = B.cols();

  for (int i = 0; i < nrowsA; ++i) {
    const double* tmpB = dataB;
    for (int j = 0; j < ncolsB; ++j) {
      const double* tmpA = dataA + i;
      double s(0.0);
      for (int k = 0; k < nrowsB; ++k) {
        s += (*tmpA) * (*tmpB);
        tmpA += ncolsA;
        tmpB++;
      }
      AB(i, j) = s;
    }
  }
}

}  // namespace AmanziChemistry
}  // namespace Amanzi
