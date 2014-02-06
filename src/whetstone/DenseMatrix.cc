/*
  This is the discretization component of the Amanzi code. 

  Copyright 2010-20XX held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Version: 2.0
  Release name: naka-to.
  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <vector>
#include <cmath>

#include "lapack.hh"
#include "WhetStoneDefs.hh"
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
* Constructor.
****************************************************************** */
DenseMatrix::DenseMatrix() 
{ 
  m_ = 0;
  n_ = 0;
  data_ = NULL; 
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
* Second level routine: max values
****************************************************************** */
void DenseMatrix::MaxRowValue(int irow, int jmin, int jmax, int* j, double* value)
{
  double* data = data_ + jmin * m_ + irow;
  *j = jmin;
  *value = *data; 

  for (int k = jmin + 1; k < jmax + 1; k++) {
    data += m_;
    if (*data > *value) {
      *j = k;
      *value = *data; 
    }
  }
}


/* ******************************************************************
* Second level routine: max absolute value
****************************************************************** */
void DenseMatrix::MaxRowMagnitude(int irow, int jmin, int jmax, int* j, double* value)
{
  double* data = data_ + jmin * m_ + irow;
  *j = jmin;
  *value = fabs(*data); 

  for (int k = jmin + 1; k < jmax + 1; k++) {
    data += m_;
    if (fabs(*data) > *value) {
      *j = k;
      *value = fabs(*data); 
    }
  }
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


/* ******************************************************************
* Second level routine: calculates matrix D such that (*this)^T D = 0.
* The matrix (*this) must have a full rank and have more rows than
* columns.
****************************************************************** */
int DenseMatrix::NullSpace(DenseMatrix& D)
{
  // We can treat only one type of rectangular matrix.
  if (m_ <= n_) return WHETSTONE_ELEMENTAL_MATRIX_FAILED;

  // D must have proper size.
  int m = D.NumRows(), n = D.NumCols();
  if (m != m_ || n != m_ - n_) return WHETSTONE_ELEMENTAL_MATRIX_FAILED;

  // Allocate memory for Lapack routine.
  int ldv(1), lwork, info; 
  lwork = std::max(m_ + 3 * n_, 5*n_);
  double U[m_ * m_], V, S[n_], work[lwork];

  DGESVD_F77("A", "N", &m_,  &n_, data_, &m_, 
             S, U, &m, &V, &ldv, work, &lwork, &info);

  if (info != 0) return WHETSTONE_ELEMENTAL_MATRIX_FAILED;
  
  double* data = D.Values();
  int offset = m_ * n_;
  for (int i = 0; i < m * n; i++) data[i] = U[offset + i];

  return 0;
}

}  // namespace WhetStone
}  // namespace Amanzi
