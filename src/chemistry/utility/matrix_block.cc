/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

#include "matrix_block.hh"

#include <cmath>

#include <iostream>
#include <iomanip>

namespace amanzi {
namespace chemistry {

MatrixBlock::MatrixBlock() : size_(0),
                             A_(NULL) {
}  // end MatrixBlock()

MatrixBlock::MatrixBlock(const int size) 
    : size_(size),
      A_(NULL) {
  AllocateMemory();
}  // end MatrixBlock(size)

MatrixBlock::~MatrixBlock() {
  FreeMemory();
}  // end ~MatrixBlock

void MatrixBlock::Resize(const int new_size) {
  FreeMemory();
  set_size(new_size);
  AllocateMemory();
}  // end Resize(new_size)

void MatrixBlock::AllocateMemory(void) {
  A_ = new double*[size()];
  for (int i = 0; i < size(); i++) {
    A_[i] = new double[size()];
  }
}  // end AllocateMemory()

void MatrixBlock::FreeMemory(void) {
  if (A_) {
    for (int i = 0; i < size(); i++) {
      delete [] A_[i];
    }
    delete [] A_;
  }
  A_ = NULL;
}  // end FreeMemory()

void MatrixBlock::Zero() {
  for (int i = 0; i < size(); i++) {
    for (int j = 0; j < size(); j++) {
      A_[i][j] = 0.;
    }
  }
}

void MatrixBlock::SetDiagonal(double d) {
  for (int i = 0; i < size(); i++) {
    A_[i][i] = d;
  }
}

double MatrixBlock::GetRowAbsMax(int irow) {
  double max = 0.;
  for (int i = 0; i < size(); i++) {
    double value = std::fabs(A_[irow][i]);
    if (value > max) {
      max = value;
    }
  }
  return max;
}

void MatrixBlock::ScaleRow(int irow, double scale) {
  for (int i = 0; i < size(); i++) {
    A_[irow][i] *= scale;
  }
}

void MatrixBlock::ScaleColumn(int icol, double scale) {
  for (int i = 0; i < size(); i++) {
    A_[i][icol] *= scale;
  }
}

void MatrixBlock::Scale(double scale) {
  for (int i = 0; i < size(); i++) {
    for (int j = 0; j < size(); j++) {
      A_[i][j] *= scale;
    }
  }
}

void MatrixBlock::SetValue(int i, int j, double value) {
  A_[i][j] = value;
}

void MatrixBlock::SetValues(double** values) {
  for (int i = 0; i < size(); i++) {
    for (int j = 0; j < size(); j++) {
      A_[i][j] = values[i][j];
    }
  }
}

void MatrixBlock::SetValues(MatrixBlock* b) {
  double** B = b->GetValues();
  for (int i = 0; i < size(); i++) {
    for (int j = 0; j < size(); j++) {
      A_[i][j] = B[i][j];
    }
  }
}

void MatrixBlock::SetValues(int ioffset, int joffset, MatrixBlock* b) {
  double** B = b->GetValues();
  for (int i = 0; i < size(); i++) {
    for (int j = 0; j < size(); j++) {
      A_[i + ioffset][j + joffset] = B[i][j];
    }
  }
}

void MatrixBlock::SetValues(MatrixBlock* b, double scale) {
  double** B = b->GetValues();
  for (int i = 0; i < size(); i++) {
    for (int j = 0; j < size(); j++) {
      A_[i][j] = scale * B[i][j];
    }
  }
}

void MatrixBlock::SetValues(double** values, double scale) {
  for (int i = 0; i < size(); i++) {
    for (int j = 0; j < size(); j++) {
      A_[i][j] = scale * values[i][j];
    }
  }
}

void MatrixBlock::SetValues(int ioffset, int joffset, MatrixBlock* b, double scale) {
  double** B = b->GetValues();
  for (int i = 0; i < size(); i++) {
    for (int j = 0; j < size(); j++) {
      A_[i + ioffset][j + joffset] = scale * B[i][j];
    }
  }
}

void MatrixBlock::AddValue(int i, int j, double value) {
  A_[i][j] += value;
}

void MatrixBlock::AddValues(double** values) {
  for (int i = 0; i < size(); i++) {
    for (int j = 0; j < size(); j++) {
      A_[i][j] += values[i][j];
    }
  }
}

void MatrixBlock::AddValues(MatrixBlock* b) {
  double** B = b->GetValues();
  for (int i = 0; i < size(); i++) {
    for (int j = 0; j < size(); j++) {
      A_[i][j] += B[i][j];
    }
  }
}

void MatrixBlock::AddValues(int ioffset, int joffset, MatrixBlock* b) {
  double** B_ = b->GetValues();
  for (int i = 0; i < size(); i++) {
    for (int j = 0; j < size(); j++) {
      A_[i + ioffset][j + joffset] += B_[i][j];
    }
  }
}

void MatrixBlock::AddValues(double** values, double scale) {
  for (int i = 0; i < size(); i++) {
    for (int j = 0; j < size(); j++) {
      A_[i][j] += scale * values[i][j];
    }
  }
}

void MatrixBlock::AddValues(MatrixBlock* b, double scale) {
  double** B = b->GetValues();
  for (int i = 0; i < size(); i++) {
    for (int j = 0; j < size(); j++) {
      A_[i][j] += scale * B[i][j];
    }
  }
}

void MatrixBlock::AddValues(int ioffset, int joffset, MatrixBlock* b, double scale) {
  double** B = b->GetValues();
  for (int i = 0; i < size(); i++) {
    for (int j = 0; j < size(); j++) {
      A_[i + ioffset][j + joffset] += scale * B[i][j];
    }
  }
}

void MatrixBlock::Print() {
  for (int i = 0; i < size(); i++) {
    for (int j = 0; j < size(); j++) {
      if (std::fabs(A_[j][i]) > 0.) {
        // TODO(bandre): is the [j][i] indexing here intentional...?
        std::cout << i << " " << j << " : "
                  << std::scientific << A_[j][i] << std::endl;
      }
    }
  }
}

void MatrixBlock::Print_ij() {
  for (int i = 0; i < size(); i++) {
    for (int j = 0; j < size(); j++) {
      if (std::fabs(A_[i][j]) > 0.) {
        std::cout << i << " " << j << " : "
                  << std::scientific << A_[i][j] << std::endl;
      }
    }
  }
}

}  // namespace chemistry
}  // namespace amanzi
