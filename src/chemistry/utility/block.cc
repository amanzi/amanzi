#include <cmath>

#include <iostream>
#include <iomanip>

#include "block.hh"

Block::Block() {
  size = 0;
  A = NULL;
}

Block::Block(int n) {
  size = n;
  A = new double*[size];
  for (int i = 0; i < size; i++)
    A[i] = new double[size];
}

void Block::zero() {
  for (int i = 0; i < size; i++) {
    for (int j = 0; j < size; j++)
      A[i][j] = 0.;
  }
}

void Block::setDiagonal(double d) {
  for (int i = 0; i < size; i++)
    A[i][i] = d;
}

double Block::getRowAbsMax(int irow) {
  double max = 0.;
  for (int i = 0; i < size; i++) {
    double value = std::fabs(A[irow][i]);
    if (value > max) max = value;
  }
  return max;
}

void Block::scaleRow(int irow, double scale) {
  for (int i = 0; i < size; i++)
    A[irow][i] *= scale;
}

void Block::scaleColumn(int icol, double scale) {
  for (int i = 0; i < size; i++)
    A[i][icol] *= scale;
}

void Block::scale(double scale) {
  for (int i = 0; i < size; i++) {
    for (int j = 0; j < size; j++)
      A[i][j] *= scale;
  }
}

void Block::setValue(int i, int j, double value) {
  A[i][j] = value;
}

void Block::setValues(double **values) {
  for (int i = 0; i < size; i++) {
    for (int j = 0; j < size; j++)
      A[i][j] = values[i][j];
  }
}

void Block::setValues(Block *b) {
  double **B = b->getValues();
  for (int i = 0; i < size; i++) {
    for (int j = 0; j < size; j++)
      A[i][j] = B[i][j];
  }
}

void Block::setValues(int ioffset, int joffset, Block *b) {
  double **A_ = b->getValues();
  for (int i = 0; i < size; i++) {
    for (int j = 0; j < size; j++)
      A[i+ioffset][j+joffset] = A_[i][j];
  }
}

void Block::setValues(Block *b, double scale) {
  double **B = b->getValues();
  for (int i = 0; i < size; i++) {
    for (int j = 0; j < size; j++)
      A[i][j] = scale*B[i][j];
  }
}

void Block::setValues(double **values, double scale) {
  for (int i = 0; i < size; i++) {
    for (int j = 0; j < size; j++)
      A[i][j] = scale*values[i][j];
  }
}

void Block::setValues(int ioffset, int joffset, Block *b, double scale) {
  double **B = b->getValues();
  for (int i = 0; i < size; i++) {
    for (int j = 0; j < size; j++)
      A[i+ioffset][j+joffset] = scale*B[i][j];
  }
}

void Block::addValue(int i, int j, double value) {
  A[i][j] += value;
}

void Block::addValues(double **values) {
  for (int i = 0; i < size; i++) {
    for (int j = 0; j < size; j++)
      A[i][j] += values[i][j];
  }
}

void Block::addValues(Block *b) {
  double **B = b->getValues();
  for (int i = 0; i < size; i++) {
    for (int j = 0; j < size; j++)
      A[i][j] += B[i][j];
  }
}

void Block::addValues(int ioffset, int joffset, Block *b) {
  double **A_ = b->getValues();
  for (int i = 0; i < size; i++) {
    for (int j = 0; j < size; j++)
      A[i+ioffset][j+joffset] += A_[i][j];
  }
}

void Block::addValues(double **values, double scale) {
  for (int i = 0; i < size; i++) {
    for (int j = 0; j < size; j++)
      A[i][j] += scale*values[i][j];
  }
}

void Block::addValues(Block *b, double scale) {
  double **B = b->getValues();
  for (int i = 0; i < size; i++) {
    for (int j = 0; j < size; j++)
      A[i][j] += scale*B[i][j];
  }
}

void Block::addValues(int ioffset, int joffset, Block *b, double scale) {
  double **B = b->getValues();
  for (int i = 0; i < size; i++) {
    for (int j = 0; j < size; j++)
      A[i+ioffset][j+joffset] += scale*B[i][j];
  }
}

void Block::print() {
  for (int i = 0; i < size; i++) {
    for (int j = 0; j < size; j++) {
      if (std::fabs(A[j][i]) > 0.) {
        std::cout << i << " " << j << " : "
                  << std::scientific << A[j][i] << std::endl;
      }
    }
  }
}

Block::~Block() {
  if (A) {
    for (int i=0; i < size; i++)
      delete [] A[i];
    delete [] A;
  }
  A = NULL;
}
