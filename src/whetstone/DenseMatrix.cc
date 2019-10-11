/*
  WhetStone, Version 2.2
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Dense matrices and operations with them.
*/

#include <vector>
#include <cmath>

#include "dbc.hh"

#include "lapack.hh"
#include "WhetStoneDefs.hh"
#include "DenseMatrix.hh"

namespace Amanzi {
namespace WhetStone {

/* ******************************************************************
 * Constructor: memory is allocated.
 ****************************************************************** */
DenseMatrix::DenseMatrix(int mrow, int ncol)
{
  m_ = mrow;
  n_ = ncol;
  mem_ = m_ * n_;
  data_ = new double[mem_];
  access_ = WHETSTONE_DATA_ACCESS_COPY;
}


/* ******************************************************************
 * Constructor: empty matrix
 ****************************************************************** */
DenseMatrix::DenseMatrix()
{
  m_ = 0;
  n_ = 0;
  mem_ = 0;
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
    mem_ = m_ * n_;
    data_ = new double[mem_];
    for (int i = 0; i < mem_; i++) data_[i] = data[i];
  } else {
    mem_ = 0;
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

  mem_ = m_ * n_;
  data_ = new double[mem_];
  const double* dataB = B.Values();
  for (int i = 0; i < mem_; i++) data_[i] = dataB[i];
}


/* ******************************************************************
 * Copy constructor creates a new matrix from the submatrix of B in
 * rows m1 to m2-1 and columns n1 to n2-1.
 ****************************************************************** */
DenseMatrix::DenseMatrix(const DenseMatrix& B, int m1, int m2, int n1, int n2)
{
  m_ = m2 - m1;
  n_ = n2 - n1;
  access_ = WHETSTONE_DATA_ACCESS_COPY;

  mem_ = m_ * n_;
  data_ = new double[mem_];

  int mB = B.NumRows();
  int nB = B.NumCols();
  const double* dataB = B.Values();

  for (int j = n1; j < n2; ++j) {
    const double* tmpB = dataB + j * mB + m1;
    double* tmpA = data_ + (j - n1) * m_;
    for (int i = 0; i < m_; ++i) {
      *tmpA = *tmpB;
      tmpA++;
      tmpB++;
    }
  }
}


/* ******************************************************************
 * Smart memory management. Data destroyed in general.
 ****************************************************************** */
void
DenseMatrix::Reshape(int mrow, int ncol)
{
  AMANZI_ASSERT(access_ == WHETSTONE_DATA_ACCESS_COPY);

  m_ = mrow;
  n_ = ncol;

  if (mem_ < m_ * n_) {
    if (data_ != NULL) { delete[] data_; }
    mem_ = m_ * n_;
    data_ = new double[mem_];
  }
}


/* ******************************************************************
 * Trace operator is extended to rectangular matrices.
 ****************************************************************** */
double
DenseMatrix::Trace()
{
  double s(0.0);
  for (int i = 0; i < m_ * n_; i += m_ + 1) s += data_[i];
  return s;
}


/* ******************************************************************
 * Matrix-matrix product. The matrices are ordered by columns.
 ****************************************************************** */
int
DenseMatrix::Multiply(const DenseMatrix& A, const DenseMatrix& B,
                      bool transposeA)
{
  const double* dataA = A.Values();
  const double* dataB = B.Values();

  int mrowsA = A.NumRows(), ncolsA = A.NumCols();
  int mrowsB = B.NumRows(), ncolsB = B.NumCols();

  if (!transposeA) {
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
int
DenseMatrix::Multiply(const DenseVector& A, DenseVector& B,
                      bool transpose) const
{
  const double* dataA = A.Values();
  double* dataB = B.Values();

  int mrowsA = A.NumRows();
  int mrowsB = B.NumRows();

  if (!transpose) {
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
void
DenseMatrix::MaxRowValue(int irow, int jmin, int jmax, int* j, double* value)
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
void
DenseMatrix::MaxRowMagnitude(int irow, int jmin, int jmax, int* j,
                             double* value)
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
 * Second level routine: submatrix in rows [ib,ie) and columns [jb,je)
 ****************************************************************** */
DenseMatrix
DenseMatrix::SubMatrix(int ib, int ie, int jb, int je)
{
  int mrows(ie - ib), ncols(je - jb);
  DenseMatrix tmp(mrows, ncols);
  double* dataB = tmp.Values();

  for (int j = jb; j < je; ++j) {
    const double* dataA = data_ + j * m_ + ib;
    for (int i = ib; i < ie; ++i) {
      *dataB = *dataA;
      dataA++;
      dataB++;
    }
  }

  return tmp;
}


/* ******************************************************************
 * Second level routine: transpose
 ****************************************************************** */
void
DenseMatrix::Transpose(const DenseMatrix& A)
{
  const double* dataA = A.Values();
  int mrowsA = A.NumRows(), ncolsA = A.NumCols();

  Reshape(ncolsA, mrowsA);

  for (int j = 0; j < ncolsA; ++j) {
    for (int i = 0; i < mrowsA; ++i) {
      *(data_ + i * ncolsA + j) = *dataA;
      dataA++;
    }
  }
}


int
DenseMatrix::Transpose()
{
  if (m_ != n_) return 1;

  for (int i = 0; i < m_; ++i) {
    for (int j = i + 1; j < n_; ++j) {
      double* aij = data_ + i * n_ + j;
      double* aji = data_ + j * n_ + i;

      double tmp(*aij);
      *aij = *aji;
      *aji = tmp;
    }
  }
  return 0;
}


/* ******************************************************************
 * Second level routine: inversion
 ****************************************************************** */
int
DenseMatrix::Inverse()
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
 * Second level routine: inversion of symmetric poitive definite
 ****************************************************************** */
int
DenseMatrix::InverseSPD()
{
  if (n_ != m_) return 911;

  int ierr;
  DPOTRF_F77("U", &n_, data_, &n_, &ierr);
  if (ierr) return ierr;

  DPOTRI_F77("U", &n_, data_, &n_, &ierr);
  for (int i = 0; i < n_; i++)
    for (int j = i + 1; j < n_; j++) data_[i * n_ + j] = data_[j * n_ + i];
  return ierr;
}


/* ******************************************************************
 * Second level routine: calculates matrix D such that (*this)^T D = 0.
 * The matrix (*this) must have a full rank and have more rows than
 * columns.
 ****************************************************************** */
int
DenseMatrix::NullSpace(DenseMatrix& D)
{
  // We can treat only one type of rectangular matrix.
  if (m_ <= n_) return WHETSTONE_ELEMENTAL_MATRIX_FAILED;

  // D must have proper size.
  int m = D.NumRows(), n = D.NumCols();
  if (m != m_ || n != m_ - n_) return WHETSTONE_ELEMENTAL_MATRIX_FAILED;

  // Allocate memory for Lapack routine.
  int ldv(1), lwork, info;
  lwork = std::max(m_ + 3 * n_, 5 * n_);
  double U[m_ * m_], V, S[n_], work[lwork];

  DGESVD_F77(
    "A", "N", &m_, &n_, data_, &m_, S, U, &m, &V, &ldv, work, &lwork, &info);

  if (info != 0) return WHETSTONE_ELEMENTAL_MATRIX_FAILED;

  double* data = D.Values();
  int offset = m_ * n_;
  for (int i = 0; i < m * n; i++) data[i] = U[offset + i];

  return 0;
}


/* ******************************************************************
 * Optimized linear algebra: determinant for matrices of size<=3.
 ****************************************************************** */
double
DenseMatrix::Det()
{
  double a = 0.0;
  if (m_ == 2) {
    a = data_[0] * data_[3] - data_[1] * data_[2];
  } else if (m_ == 3) {
    a = data_[0] * data_[4] * data_[8] + data_[3] * data_[7] * data_[2] +
        data_[1] * data_[5] * data_[6] - data_[6] * data_[4] * data_[2] -
        data_[3] * data_[1] * data_[8] - data_[0] * data_[7] * data_[5];
  } else {
    a = data_[0];
  }
  return a;
}


/* ******************************************************************
 * Orthonormalize selected matrix columns.
 ****************************************************************** */
int
DenseMatrix::OrthonormalizeColumns(int n1, int n2)
{
  double sum;

  for (int i = n1; i < n2; ++i) {
    double* tmp1 = data_ + i * m_;

    for (int j = n1; j < i; ++j) {
      double* tmp2 = data_ + j * m_;

      sum = 0.0;
      for (int k = 0; k < m_; ++k) sum += tmp1[k] * tmp2[k];
      for (int k = 0; k < m_; ++k) tmp1[k] -= sum * tmp2[k];
    }

    sum = 0.0;
    for (int k = 0; k < m_; ++k) sum += tmp1[k] * tmp1[k];
    if (sum == 0.0) return -1;

    sum = std::pow(1.0 / sum, 0.5);
    for (int k = 0; k < m_; ++k) tmp1[k] *= sum;
  }

  return 0;
}


/* ******************************************************************
 * Permutation of matrix columns.
 ****************************************************************** */
void
DenseMatrix::SwapColumns(int n1, int n2)
{
  if (n1 != n2) {
    double tmp;
    double* row1 = data_ + n1 * m_;
    double* row2 = data_ + n2 * m_;

    for (int i = 0; i < m_; i++) {
      tmp = row1[i];
      row1[i] = row2[i];
      row2[i] = tmp;
    }
  }
}

} // namespace WhetStone
} // namespace Amanzi
