/*
 This is the discretization component of the Amanzi code. 

 Copyright 2010-20XX held jointly by LANS/LANL, LBNL, and PNNL. 
 Amanzi is released under the three-clause BSD License. 
 The terms of use and "as is" disclaimer for this license are 
 provided in the top-level COPYRIGHT file.

Version: 2.0
Release name: naka-to.
Author: Konstantin Lipnikov (lipnikov@lanl.gov)

Usage: 
*/

#include <vector>

#include "lapack.hh"
#include "DenseMatrix.hh"

namespace Amanzi {
namespace WhetStone {

/* ******************************************************************
* Constructor.
****************************************************************** */
DenseMatrix::DenseMatrix(int mrow, int ncol) 
{ 
  m_ = mrow;
  n_ = ncol;
  data_ = new double[m_ * n_]; 
  access_ = WHETSTONE_DATA_ACCESS_COPY;
}


/* ******************************************************************
* No memory check is performed: invalid read is possible.
****************************************************************** */
DenseMatrix::DenseMatrix(int mrow, int ncol, double* data, int data_access)
{
  m_ = mrow;
  n_ = ncol;
  access_ = data_access;

  if (access_ == WHETSTONE_DATA_ACCESS_COPY) {
    data_ = new double[m_ * n_]; 
    for (int i = 0; i < mrow * ncol; i++) data_[i] = data[i];
  } else {
    data_ = data;
  }
}


/* ******************************************************************
* Copy constructor creates an new matrix.
****************************************************************** */
DenseMatrix::DenseMatrix(const DenseMatrix& B)
{
  m_ = B.NumRows();
  n_ = B.NumCols();
  access_ = WHETSTONE_DATA_ACCESS_COPY;

  data_ = new double[m_ * n_];
  const double* dataB = B.Values();
  for (int i = 0; i < m_ * n_; i++) data_[i] = dataB[i];
}


/* ******************************************************************
* Matrix-matrix product. The matrices are ordered by columns.
****************************************************************** */
int DenseMatrix::Multiply(const DenseMatrix& A, 
                          const DenseMatrix& B, bool transposeA)
{
  const double* dataA = A.Values();
  const double* dataB = B.Values();

  int mrowsA = A.NumRows(), ncolsA = A.NumCols();
  int mrowsB = B.NumRows(), ncolsB = B.NumCols();

  if (! transposeA) {
    if (ncolsA != mrowsB || m_ != mrowsA || n_ != ncolsB) return 1;

    for (int i = 0; i < m_; i++) {
      const double* tmpB = dataB;
      for (int j = 0; j < n_; j++) {
        const double* tmpA = dataA + i;
        double s(0.0);
        for (int k = 0; k < mrowsB; k++) {
          s += (*tmpA) * (*tmpB);
          tmpA += m_;
          tmpB++;
        }
        (*this)(i, j) = s;
      }
    } 
  } else {
    if (mrowsA != mrowsB || m_ != ncolsA || n_ != ncolsB) return 1;

    for (int i = 0; i < m_; i++) {
      const double* tmpB = dataB;
      for (int j = 0; j < n_; j++) {
        const double* tmpA = dataA + i * mrowsA;
        double s(0.0);
        for (int k = 0; k < mrowsB; k++) {
          s += (*tmpA) * (*tmpB);
          tmpA++;
          tmpB++;
        }
        (*this)(i, j) = s;
      }
    }
  }

  return 0;
}


/* ******************************************************************
* Matrix-vector product. The matrix is ordered by columns.
****************************************************************** */
int DenseMatrix::Multiply(const DenseVector& A, DenseVector& B, bool transpose)
{
  const double* dataA = A.Values();
  double* dataB = B.Values();

  int mrowsA = A.NumRows();
  int mrowsB = B.NumRows();

  if (! transpose) {
    if (n_ != mrowsA || m_ != mrowsB) return 1;

    for (int i = 0; i < m_; i++) {
      const double* tmpM = data_ + i;
      double s(0.0);
      for (int j = 0; j < n_; j++) {
        s += (*tmpM) * dataA[j];
        tmpM += m_;
      }
      dataB[i] = s;
    } 
  } else {
    if (m_ != mrowsA || n_ != mrowsB) return 1;

    const double* tmpM = data_;
    for (int i = 0; i < n_; i++) {
      double s(0.0);
      for (int j = 0; j < m_; j++) {
        s += (*tmpM) * dataA[j];
        tmpM++;
      }
      dataB[i] = s;
    }
  }

  return 0;
}


/* ******************************************************************
* Second level routine: inversion
****************************************************************** */
int DenseMatrix::Inverse() 
{
  if (n_ != m_) return 911;

  int ierr;
  int iwork[n_];
  DGETRF_F77(&n_, &n_, data_, &n_, iwork, &ierr);
  if (ierr) return ierr;
 
  int lwork = n_ * n_;
  double dwork[lwork];
  DGETRI_F77(&n_, data_, &n_, iwork, dwork, &lwork, &ierr);
  return ierr;
}

}  // namespace WhetStone
}  // namespace Amanzi
