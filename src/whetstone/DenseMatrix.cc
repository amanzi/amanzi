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
DenseMatrix::DenseMatrix(const DenseMatrix& other)
{
  m_ = other.NumRows();
  n_ = other.NumCols();
  access_ = WHETSTONE_DATA_ACCESS_COPY;

  mem_ = m_ * n_;
  data_ = new double[mem_];
  const double* dataB = other.Values();
  for (int i = 0; i < mem_; i++) data_[i] = dataB[i];
}


/* ******************************************************************
* Move constructor re-assigns memory
****************************************************************** */
DenseMatrix::DenseMatrix(DenseMatrix&& other) noexcept
{
  m_ = other.m_;
  n_ = other.n_;
  mem_ = other.mem_;
  access_ = other.access_;
  data_ = other.data_;
  other.data_ = NULL;
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
* Assignment operators
****************************************************************** */
DenseMatrix& DenseMatrix::operator=(const DenseMatrix& other)
{
  if (this != &other) {
    if (mem_ < other.m_ * other.n_) {
      if (data_ != NULL) {
        delete [] data_;
      }
      mem_ = other.m_ * other.n_;
      data_ = new double[mem_];
    }
    n_ = other.n_;
    m_ = other.m_;
    const double *b = other.Values();
    for (int i = 0; i < m_ * n_; i++) data_[i] = b[i];
  }
  return *this;
}


DenseMatrix& DenseMatrix::operator=(DenseMatrix&& other) noexcept
{
  if (this != &other) {
    n_ = other.n_;
    m_ = other.m_;
    mem_ = other.mem_;
    access_ = other.access_;
    data_ = other.data_;
    other.data_ = NULL;
  }
  return *this;
}


/* ******************************************************************
* Smart memory management. Data destroyed in general.
****************************************************************** */
void DenseMatrix::Reshape(int mrow, int ncol)
{
  AMANZI_ASSERT(access_ == WHETSTONE_DATA_ACCESS_COPY);

  m_ = mrow;
  n_ = ncol;

  if (mem_ < m_ * n_) {
    if (data_ != NULL) {
      delete [] data_;
    }
    mem_ = m_ * n_;
    data_ = new double[mem_];
  }
}


/* ******************************************************************
* Trace operator is extended to rectangular matrices.
****************************************************************** */
double DenseMatrix::Trace()
{
  double s(0.0);
  for (int i = 0; i < m_ * n_; i += m_ + 1) s += data_[i];
  return s;
}


/* ******************************************************************
* Ring algebra
****************************************************************** */
DenseMatrix operator+(const DenseMatrix& A, const DenseMatrix& B)
{
  int m = A.NumRows();
  int n = A.NumCols();
  AMANZI_ASSERT(m == B.NumRows() && n == B.NumCols()); 

  DenseMatrix AB(m, n);
  
  const double* dataA = A.Values();
  const double* dataB = B.Values();
  double* dataAB = AB.Values();

  for (int i = 0; i < m * n; ++i) {
    *dataAB = *dataA + *dataB;
    dataA++;
    dataB++;
    dataAB++;
  }
  return AB;
}


DenseMatrix operator-(const DenseMatrix& A, const DenseMatrix& B)
{
  int m = A.NumRows();
  int n = A.NumCols();
  AMANZI_ASSERT(m == B.NumRows() && n == B.NumCols()); 

  DenseMatrix AB(m, n);
  
  const double* dataA = A.Values();
  const double* dataB = B.Values();
  double* dataAB = AB.Values();

  for (int i = 0; i < m * n; ++i) {
    *dataAB = *dataA - *dataB;
    dataA++;
    dataB++;
    dataAB++;
  }
  return AB;
}


DenseMatrix operator*(const DenseMatrix& A, const DenseMatrix& B)
{
  const double* dataA = A.Values();
  const double* dataB = B.Values();

  int mrowsA = A.NumRows();
  int mrowsB = B.NumRows(), ncolsB = B.NumCols();

  DenseMatrix AB(mrowsA, ncolsB);

  for (int i = 0; i < mrowsA; i++) {
    const double* tmpB = dataB;
    for (int j = 0; j < ncolsB; j++) {
      const double* tmpA = dataA + i;
      double s(0.0);
      for (int k = 0; k < mrowsB; k++) {
        s += (*tmpA) * (*tmpB);
        tmpA += mrowsA;
        tmpB++;
      }
      AB(i, j) = s;
    }
  } 
  return AB;
}


/* ******************************************************************
* Kroneker product (A ^ K)
****************************************************************** */
DenseMatrix operator^(const DenseMatrix& A, const Tensor& K)
{
  int d = K.dimension();
  int nrows = A.NumRows();
  int ncols = A.NumCols();

  DenseMatrix AK(nrows * d, ncols * d);

  for (int m = 0; m < nrows; ++m) {
    for (int n = 0; n < ncols; ++n) {
      for (int i = 0; i < d; ++i) {
        for (int j = 0; j < d; ++j) {
          AK(d * m + i, d * n + j) = A(m, n) * K(i, j);
        }
      }
    }
  }

  return AK;
}


/* ******************************************************************
* Kroneker product (K ^ A)
****************************************************************** */
DenseMatrix operator^(const Tensor& K, const DenseMatrix& A)
{
  int d = K.dimension();
  int nrows = A.NumRows();
  int ncols = A.NumCols();

  DenseMatrix AK(nrows * d, ncols * d);

  for (int m = 0; m < nrows; ++m) {
    for (int n = 0; n < ncols; ++n) {
      for (int i = 0; i < d; ++i) {
        for (int j = 0; j < d; ++j) {
          AK(nrows * i + m, ncols * j + n) = A(m, n) * K(i, j);
        }
      }
    }
  }

  return AK;
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
int DenseMatrix::Multiply(const DenseVector& A, DenseVector& B, bool transpose) const
{
  const double* dataA = A.Values();
  double* dataB = B.Values();

  int mrowsA = A.NumRows();
  int mrowsB = B.NumRows();

  return Multiply_(dataA, mrowsA, dataB, mrowsB, transpose);
}


/* ******************************************************************
* Block matrix-vector product: single matrix, multiple subvectors.
****************************************************************** */
int DenseMatrix::BlockMultiply(const DenseVector& A, DenseVector& B, bool transpose) const
{
  const double* dataA = A.Values();
  double* dataB = B.Values();

  int mrowsA = A.NumRows();
  int mrowsB = B.NumRows();

  int k1, k2;
  if (! transpose) {
    if (mrowsA % n_ != 0) return 1;
    if (mrowsB % m_ != 0) return 1;

    k1 = mrowsA / n_;
    k2 = mrowsB / m_;
  } else {
    if (mrowsA % m_ != 0) return 1;
    if (mrowsB % n_ != 0) return 1;

    k1 = mrowsA / m_;
    k2 = mrowsB / n_;
  }
  if (k1 != k2) return 1;

  mrowsA /= k1;
  mrowsB /= k2;

  int i(0);
  for (int k = 0; k < k1; ++k) {
    i |= Multiply_(dataA, mrowsA, dataB, mrowsB, transpose);
    dataA += mrowsA;
    dataB += mrowsB;
  }
  return i;
}


/* ******************************************************************
* Matrix-vector product. The matrix is ordered by columns.
****************************************************************** */
int DenseMatrix::Multiply_(
    const double* A, int mrowsA, double* B, int mrowsB, bool transpose) const
{
  if (! transpose) {
    if (n_ != mrowsA || m_ != mrowsB) return 1;

    for (int i = 0; i < m_; i++) {
      const double* tmpM = data_ + i;
      double s(0.0);
      for (int j = 0; j < n_; j++) {
        s += (*tmpM) * A[j];
        tmpM += m_;
      }
      B[i] = s;
    } 
  } else {
    if (m_ != mrowsA || n_ != mrowsB) return 1;

    const double* tmpM = data_;
    for (int i = 0; i < n_; i++) {
      double s(0.0);
      for (int j = 0; j < m_; j++) {
        s += (*tmpM) * A[j];
        tmpM++;
      }
      B[i] = s;
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
* Second level routine: submatrix in rows [ib,ie) and columns [jb,je)
****************************************************************** */
DenseMatrix DenseMatrix::SubMatrix(int ib, int ie, int jb, int je) 
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
* Second level routine: insert submatrix of B at position (is, js)
****************************************************************** */
void DenseMatrix::InsertSubMatrix(
    const DenseMatrix& B, int ib, int ie, int jb, int je, int is, int js)
{
  int ks(js), mB(B.NumRows());

  for (int j = jb; j < je; ++j) {
    double* dataA = data_ + ks * m_ + is;
    const double* dataB = B.Values() + j * mB + ib;
    for (int i = ib; i < ie; ++i) {
      *dataA = *dataB;
      dataA++; 
      dataB++; 
    }
    ks++;
  } 
}


/* ******************************************************************
* Second level routine: transpose
****************************************************************** */
void DenseMatrix::Transpose(const DenseMatrix& A) 
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


int DenseMatrix::Transpose() 
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
* Second level routine: inversion of symmetric positive definite
****************************************************************** */
int DenseMatrix::InverseSPD() 
{
  if (n_ != m_) return 911;

  int ierr;
  DPOTRF_F77("U", &n_, data_, &n_, &ierr); 
  if (ierr) return ierr;
 
  DPOTRI_F77("U", &n_, data_, &n_, &ierr);
  for (int i = 0; i < n_; i++)
    for (int j = i + 1; j < n_; j++)
      data_[i * n_ + j] = data_[j * n_ + i];
  return ierr;
}


/* ******************************************************************
* Second level routine: calculates matrix D such that (*this)^T D = 0.
* The matrix (*this) must have a full rank and have more rows than
* columns.
****************************************************************** */
DenseMatrix DenseMatrix::NullSpace()
{
  // We can treat only one type of rectangular matrix.
  if (m_ <= n_) AMANZI_ASSERT(false);

  // D must have proper size.
  DenseMatrix D(m_, m_ - n_);

  // Allocate memory for Lapack routine.
  int ldv(1), lwork, info; 
  lwork = m_ + 5 * n_;
  double U[m_ * m_], V[ldv], S[n_], work[lwork];

  DGESVD_F77("A", "N", &m_,  &n_, data_, &m_, 
             S, U, &m_, V, &ldv, work, &lwork, &info);

  if (info != 0) AMANZI_ASSERT(false);
  
  double* data = D.Values();
  int offset = m_ * n_;
  for (int i = 0; i < m_ * (m_ - n_); i++) data[i] = U[offset + i];

  return D;
}


/* ******************************************************************
* Second level routine: calculates pseudo-inverse using Moore-Penrose
* algorithm.
****************************************************************** */
int DenseMatrix::InverseMoorePenrose()
{
  // We can treat only square matrices.
  if (m_ != n_) return WHETSTONE_ELEMENTAL_MATRIX_FAILED;

  // Allocate memory for Lapack routine.
  int mn, lwork, info; 
  mn = std::min(n_, m_);
  lwork = std::max(m_ + 3 * n_, 5*n_);
  double U[m_ * m_], V[n_ * n_], S[mn], work[lwork];

  DGESVD_F77("A", "A", &m_, &n_, data_, &m_, 
             S, U, &m_, V, &n_, work, &lwork, &info);

  if (info != 0) return WHETSTONE_ELEMENTAL_MATRIX_FAILED;

  // inverse of S
  for (int i = 0; i < mn; ++i)
    if (fabs(S[i]) > 1e-12) S[i] = 1.0 / S[i];

  // assemble preuso-inverse
  DenseMatrix UU(m_, m_, U, WHETSTONE_DATA_ACCESS_VIEW);
  DenseMatrix VV(n_, n_, V, WHETSTONE_DATA_ACCESS_VIEW);

  for (int i = 0; i < n_; ++i) {
    for (int j = 0; j < m_; ++j) {
      double tmp(0.0);
      for (int k = 0; k < m_; ++k)
        tmp += VV(k, i) * S[k] * UU(j, k);

      data_[j * m_ + i] = tmp;
    }
  }

  return 0;
}


/* ******************************************************************
* Optimized linear algebra: determinant for matrices of size<=3.
****************************************************************** */
double DenseMatrix::Det()
{
  double a = 0.0;
  if (m_ == 2) {
    a = data_[0] * data_[3] - data_[1] * data_[2];
  } else if (m_ == 3) {
    a = data_[0] * data_[4] * data_[8]
      + data_[3] * data_[7] * data_[2]
      + data_[1] * data_[5] * data_[6]
      - data_[6] * data_[4] * data_[2]
      - data_[3] * data_[1] * data_[8]
      - data_[0] * data_[7] * data_[5];
  } else {
    a = data_[0];
  }
  return a;
}


/* ******************************************************************
* Orthonormalize selected matrix columns.
****************************************************************** */
int DenseMatrix::OrthonormalizeColumns(int n1, int n2)
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
void DenseMatrix::SwapColumns(int n1, int n2) 
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

}  // namespace WhetStone
}  // namespace Amanzi
