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
    A_(NULL) {
}


MatrixBlock::MatrixBlock(const int size) : size_(size) {
  A_ = new double[size_ * size_];
}


MatrixBlock::~MatrixBlock() {
  if (A_) delete [] A_;
  A_ = NULL;
}


void MatrixBlock::Resize(const int new_size) {
  if (A_) delete [] A_;
  size_ = new_size;
  A_ = new double[size_ * size_];
}


void MatrixBlock::Zero() {
  for (int i = 0; i < size_ * size_; i++) A_[i] = 0.0;
}


void MatrixBlock::SetDiagonal(double d) {
  for (int i = 0; i < size(); i++) {
    A_[i * size_ + i] = d;
  }
}


double MatrixBlock::GetRowAbsMax(int irow) {
  double max = 0.;
  for (int i = 0; i < size(); i++) {
    double value = std::fabs(A_[i * size_ + irow]);
    if (value > max) {
      max = value;
    }
  }
  return max;
}


void MatrixBlock::ScaleRow(int irow, double scale) {
  for (int i = 0; i < size(); i++) {
    A_[i * size_ + irow] *= scale;
  }
}


void MatrixBlock::ScaleColumn(int icol, double scale) {
  for (int i = 0; i < size(); i++) {
    A_[icol * size_ + i] *= scale;
  }
}


void MatrixBlock::Scale(double scale) {
  for (int i = 0; i < size_ * size_; i++) A_[i] *= scale;
}


void MatrixBlock::AddValue(int i, int j, double value) {
  A_[j * size_ + i] += value;
}


void MatrixBlock::AddValues(MatrixBlock* b) {
  double* B = b->GetValues();
  for (int i = 0; i < size_ * size_; i++) A_[i] += B[i];
}


void MatrixBlock::AddValues(MatrixBlock* b, double scale) {
  double* B = b->GetValues();
  for (int i = 0; i < size_ * size_; i++) A_[i] += scale * B[i];
}

}  // namespace AmanziChemistry
}  // namespace Amanzi
