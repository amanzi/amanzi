/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  WhetStone, Version 2.2
  Release name: naka-to.

  Dense matrices and operations with them.
*/

#ifndef AMANZI_DENSE_MATRIX_HH_
#define AMANZI_DENSE_MATRIX_HH_

#include <cstdlib>
#include <iostream>
#include <cstdio>
#include <iomanip>
#include <Kokkos_Core.hpp>

#include "lapack.hh"
#include "DenseVector.hh"
#include "WhetStoneDefs.hh"

namespace Amanzi {
namespace WhetStone {


template <class MEMSPACE = DefaultHostMemorySpace>
class DenseMatrix {
 public:
  KOKKOS_INLINE_FUNCTION DenseMatrix() : m_(0), n_(0), mem_(0) {}

  DenseMatrix(const int& mrow, const int& ncol) : m_(mrow), n_(ncol), mem_(m_ * n_)
  {
    Kokkos::resize(data_, mem_);
  } // memory is not initialized

  template <class MEMSPACE_OTHER>
  DenseMatrix(const DenseMatrix<MEMSPACE_OTHER>& other)
    : m_(other.m_), n_(other.n_), mem_(other.mem_)
  {
    if (mem_ != 0) {
      Kokkos::resize(data_, mem_);
      Kokkos::deep_copy(data_, other.Values());
    }
  }

  DenseMatrix(const DenseMatrix<MEMSPACE>& other) : m_(other.m_), n_(other.n_), mem_(other.mem_)
  {
    if (mem_ != 0) {
      Kokkos::resize(data_, mem_);
      Kokkos::deep_copy(data_, other.Values());
    }
  }

#if WHETSTONE_VIEW_SEMANTICS
  KOKKOS_INLINE_FUNCTION
  DenseMatrix(Kokkos::View<double*, MEMSPACE> data, const int& mrow, const int& ncol)
    : m_(mrow), n_(ncol), mem_(m_ * n_), data_(data)
  {}

  KOKKOS_INLINE_FUNCTION DenseMatrix(DenseMatrix<MEMSPACE>&& other)
    : m_(other.m_), n_(other.n_), mem_(other.mem_), data_(other.data_)
  {}

  // Default assigment implies view semantics
  KOKKOS_INLINE_FUNCTION DenseMatrix<MEMSPACE>& operator=(const DenseMatrix<MEMSPACE>&& other)
  {
    m_ = other.m_;
    n_ = other.n_;
    data_ = other.data_;
    mem_ = other.mem_;
    return *this;
  }

  // Default assigment implies view semantics
  DenseMatrix<MEMSPACE>& operator=(const DenseMatrix<MEMSPACE>& B) = default;

#else
  DenseMatrix<MEMSPACE>& operator=(const DenseMatrix<MEMSPACE>& B) = delete;

#endif

  /* ******************************************************************
   * No memory check is performed: invalid read is possible.
   ****************************************************************** */
  DenseMatrix(int mrow, int ncol, double* data)
  {
    m_ = mrow;
    n_ = ncol;
    mem_ = m_ * n_;
    Kokkos::resize(data_, mem_);
    for (int i = 0; i < mem_; i++) data_[i] = data[i];
  }

  DenseMatrix(int mrow, int ncol, const Kokkos::View<double*, MEMSPACE>& data)
  {
    m_ = mrow;
    n_ = ncol;
    mem_ = m_ * n_;
    Kokkos::resize(data_, mem_);
    Kokkos::deep_copy(data_, data);
  }

  /* ******************************************************************
   * Copy constructor creates a new matrix from the submatrix of B in
   * rows m1 to m2-1 and columns n1 to n2-1.
   ****************************************************************** */
  DenseMatrix(const DenseMatrix<MEMSPACE>& B, int m1, int m2, int n1, int n2)
  {
    m_ = m2 - m1;
    n_ = n2 - n1;
    // access_ = WHETSTONE_DATA_ACCESS_COPY;

    mem_ = m_ * n_;
    Kokkos::resize(data_, mem_);
    // data_ = new double[mem_];

    int mB = B.NumRows();
    int nB = B.NumCols();
    // const double* dataB = B.Values();

    for (int j = n1; j < n2; ++j) {
      // const double* tmpB = dataB + j * mB + m1;
      // double* tmpA = data_ + (j - n1) * m_;
      for (int i = m1; i < m2; ++i) {
        operator()(i, j) = B(i, j);
        //*tmpA = *tmpB;
        // tmpA++;
        // tmpB++;
      }
    }
  }
  /* ******************************************************************
   * Smart memory management. Data destroyed in general.
   ****************************************************************** */
  void reshape(int mrow, int ncol)
  {
    m_ = mrow;
    n_ = ncol;
    mem_ = m_ * n_;
    Kokkos::resize(data_, mem_);
  }

  KOKKOS_INLINE_FUNCTION
  void assign(const DenseMatrix<MEMSPACE>& other)
  {
    if (this != &other) {
      assert(n_ == other.n_ && m_ == other.m_);
      //if (n_ != other.n_ || m_ != other.m_) reshape(other.m_, other.n_);
      for (int i = 0; i < n_ * m_; ++i) { data_(i) = other.data_(i); }
    }
  }

  KOKKOS_INLINE_FUNCTION double& operator()(int i, int j) { return data_(j * m_ + i); }
  KOKKOS_INLINE_FUNCTION const double& operator()(int i, int j) const { return data_(j * m_ + i); }


  KOKKOS_INLINE_FUNCTION DenseMatrix<MEMSPACE>& operator=(double val)
  {
    for (int i = 0; i < m_ * n_; i++) data_(i) = val;
    return *this;
  }

  KOKKOS_INLINE_FUNCTION DenseMatrix<MEMSPACE>& operator*=(double val)
  {
    for (int i = 0; i < m_ * n_; i++) data_(i) *= val;
    return *this;
  }

  KOKKOS_INLINE_FUNCTION DenseMatrix<MEMSPACE>& operator/=(double val)
  {
    for (int i = 0; i < m_ * n_; i++) data_(i) /= val;
    return *this;
  }

  KOKKOS_INLINE_FUNCTION DenseMatrix<MEMSPACE>& operator+=(const DenseMatrix<MEMSPACE>& A)
  {
    assert(m_ == A.m_ && n_ == A.n_);
    for (int j = 0; j < n_; j++) {
      for (int i = 0; i < m_; i++) { operator()(i, j) += A(i, j); }
    }
    return *this;
  }

  KOKKOS_INLINE_FUNCTION DenseMatrix<MEMSPACE>& operator-=(const DenseMatrix<MEMSPACE>& A)
  {
    assert(m_ == A.m_ && n_ == A.n_);
    for (int j = 0; j < n_; j++) {
      for (int i = 0; i < m_; i++) { operator()(i, j) -= A(i, j); }
    }
    return *this;
  }

  /* ******************************************************************
   * Matrix-matrix product. The matrices are ordered by columns.
   ****************************************************************** */

  // calculates either A * B to A^T * B
  KOKKOS_INLINE_FUNCTION int
  Multiply(const DenseMatrix<MEMSPACE>& A, const DenseMatrix<MEMSPACE>& B, bool transposeA)
  {
    if (!transposeA)
      assert(m_ == A.m_ && n_ == B.n_ && A.n_ == B.m_);
    else
      assert(m_ == A.n_ && n_ == B.n_ && A.m_ == B.m_);
    int mrowsA = A.NumRows(), ncolsA = A.NumCols();
    int mrowsB = B.NumRows(), ncolsB = B.NumCols();

    if (!transposeA) {
      if (ncolsA != mrowsB || m_ != mrowsA || n_ != ncolsB) return 1;

      for (int i = 0; i < m_; i++) {
        // const double* tmpB = dataB;
        for (int j = 0; j < n_; j++) {
          // const double* tmpA = dataA + i;
          double s(0.0);
          for (int k = 0; k < mrowsB; k++) {
            s += A(i, k) * B(k, j);
            // s += (*tmpA) * (*tmpB);
            // tmpA += m_;
            // tmpB++;
          }
          (*this)(i, j) = s;
        }
      }
    } else {
      if (mrowsA != mrowsB || m_ != ncolsA || n_ != ncolsB) return 1;

      for (int i = 0; i < m_; i++) {
        // const double* tmpB = dataB;
        for (int j = 0; j < n_; j++) {
          // const double* tmpA = dataA + i * mrowsA;
          double s(0.0);
          for (int k = 0; k < mrowsB; k++) {
            s += A(k, i) * B(k, j);
            // s += (*tmpA) * (*tmpB);
            // tmpA++;
            // tmpB++;
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
  // calculates B = *this * A
  KOKKOS_INLINE_FUNCTION
  int Multiply(const DenseVector<MEMSPACE>& A, DenseVector<MEMSPACE>& B, bool transpose) const
  {
    if (transpose) {
      assert(B.NumRows() == NumCols());
      assert(A.NumRows() == NumRows());
    } else {
      assert(A.NumRows() == NumCols());
      assert(B.NumRows() == NumRows());
    }

    int mrowsA = A.NumRows();
    int mrowsB = B.NumRows();

    if (!transpose) {
      if (n_ != mrowsA || m_ != mrowsB) return 1;

      for (int i = 0; i < m_; i++) {
        // const double* tmpM = data_ + i;
        double s(0.0);
        for (int j = 0; j < n_; j++) {
          s += operator()(i, j) * A(j);
          // s += (*tmpM) * dataA[j];
          // tmpM += m_;
        }
        B(i) = s;
      }
    } else {
      if (m_ != mrowsA || n_ != mrowsB) return 1;

      // const double* tmpM = data_;
      for (int i = 0; i < n_; i++) {
        double s(0.0);
        for (int j = 0; j < m_; j++) {
          s += operator()(j, i) * A(j);
          // s += (*tmpM) * dataA[j];
          // tmpM++;
        }
        B(i) = s;
      }
    }

    return 0;
  }

  KOKKOS_INLINE_FUNCTION void putScalar(double val)
  {
    for (int i = 0; i < m_ * n_; i++) data_(i) = val;
  }

  // access: the data are ordered by columns
  KOKKOS_INLINE_FUNCTION int NumRows() const { return m_; }
  KOKKOS_INLINE_FUNCTION int NumCols() const { return n_; }

  KOKKOS_INLINE_FUNCTION Kokkos::View<double*, MEMSPACE> Values() { return data_; }
  KOKKOS_INLINE_FUNCTION double* Values_ptr() { return &data_[0]; }
  KOKKOS_INLINE_FUNCTION const double* Values_ptr() const { return &data_[0]; }

  KOKKOS_INLINE_FUNCTION double& Value(int i, int j)
  {
#if WHETSTONE_ABC
    AMANZI_ASSERT(0 <= i && i < NumRows());
    AMANZI_ASSERT(0 <= j && j < NumCols());
#endif
    return data_(j * m_ + i);
  }
  KOKKOS_INLINE_FUNCTION const Kokkos::View<double*, MEMSPACE> Values() const { return data_; }

  KOKKOS_INLINE_FUNCTION const double& Value(int i, int j) const
  {
#if WHETSTONE_ABC
    AMANZI_ASSERT(0 <= i && i < NumRows());
    AMANZI_ASSERT(0 <= j && j < NumCols());
#endif
    return data_(j * m_ + i);
  }
  // KOKKOS_INLINE_FUNCTION double* Values_ptr() { return &data_(0); }
  // KOKKOS_INLINE_FUNCTION const double* Values_ptr() const { return &data_(0); }


  /* ******************************************************************
   * Trace operator is extended to rectangular matrices.
   ****************************************************************** */
  KOKKOS_INLINE_FUNCTION double Trace()
  {
    double s(0.0);
    for (int i = 0; (i < m_) && (i < n_); ++i) s += operator()(i, i);
    return s;
  }

  // -- extrema in rows and columns
  KOKKOS_INLINE_FUNCTION void MaxRowValue(int irow, int* j, double* value)
  {
    MaxRowValue(irow, 0, n_, j, value);
  }

  /*******************************************************************
   * Second level routine: max values
   ****************************************************************** */
  KOKKOS_INLINE_FUNCTION void MaxRowValue(int irow, int jmin, int jmax, int* j, double* value)
  {
    // double* data = data_ + jmin * m_ + irow;
    //*j = jmin;
    //*value = *data;
    *j = jmin;
    *value = operator()(irow, jmin);

    for (int k = jmin + 1; k < jmax + 1; k++) {
      // data += m_;
      double v = operator()(irow, k);
      if (v > *value) {
        // if (*data > *value) {
        *j = k;
        *value = v;
        //*value = *data;
      }
    }
  }

  KOKKOS_INLINE_FUNCTION void MaxRowMagnitude(int irow, int* j, double* value)
  {
    MaxRowMagnitude(irow, 0, n_, j, value);
  }

  /* ******************************************************************
   * Second level routine: max absolute value
   ****************************************************************** */
  KOKKOS_INLINE_FUNCTION void MaxRowMagnitude(int irow, int jmin, int jmax, int* j, double* value)
  {
    // double* data = data_ + jmin * m_ + irow;
    *j = jmin;
    *value = fabs(operator()(irow, jmin));

    for (int k = jmin + 1; k < jmax + 1; k++) {
      // data += m_;
      double v = fabs(operator()(irow, k));
      if (v > *value) {
        *j = k;
        *value = v;
      }
    }
  }

  KOKKOS_INLINE_FUNCTION double normInf() const
  {
    double a = 0.0;
    for (int i = 0; i < m_ * n_; i++)
      if (a < data_[i]) a = data_[i];
    return a;
  }

  KOKKOS_INLINE_FUNCTION double Norm2() const
  {
    double a = 0.0;
    for (int i = 0; i < m_ * n_; i++) a += data_[i] * data_[i];
    return sqrt(a);
  }

  KOKKOS_INLINE_FUNCTION void Scale(double value)
  {
    for (int i = 0; i < m_ * n_; i++) data_[i] *= value;
  }

  /* ******************************************************************
   * Optimized linear algebra: determinant for matrices of size<=3.
   ****************************************************************** */
  KOKKOS_INLINE_FUNCTION double Det()
  {
    assert(n_ == m_);
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
  KOKKOS_INLINE_FUNCTION int OrthonormalizeColumns(int n1, int n2)
  {
    double sum;

    for (int i = n1; i < n2; ++i) {
      // double* tmp1 = data_ + i * m_;

      for (int j = n1; j < i; ++j) {
        // double* tmp2 = data_ + j * m_;

        sum = 0.0;
        for (int k = 0; k < m_; ++k)
          sum += operator()(k, i) * operator()(k, j); // tmp1[k] * tmp2[k];
        for (int k = 0; k < m_; ++k) operator()(k, i) -= sum * operator()(k, j);
      }

      sum = 0.0;
      for (int k = 0; k < m_; ++k) sum += operator()(k, i) * operator()(k, i);
      if (sum == 0.0) return -1;

      sum = std::pow(1.0 / sum, 0.5);
      for (int k = 0; k < m_; ++k) operator()(k, i) *= sum; // tmp1[k] *= sum;
    }

    return 0;
  }

  /* ******************************************************************
   * Permutation of matrix columns.
   ****************************************************************** */
  KOKKOS_INLINE_FUNCTION void SwapColumns(int n1, int n2)
  {
    if (n1 != n2) {
      double tmp;
      // double* row1 = data_ + n1 * m_;
      // double* row2 = data_ + n2 * m_;

      for (int i = 0; i < m_; i++) {
        tmp = operator()(i, n1);
        operator()(i, n1) = operator()(i, n2);
        operator()(i, n2) = tmp;
        // row1[i] = row2[i];
        // row2[i] = tmp;
      }
    }
  }

  /* ******************************************************************
   * Second level routine: submatrix in rows [ib,ie) and columns [jb,je)
   ****************************************************************** */
  DenseMatrix<MEMSPACE> SubMatrix(int ib, int ie, int jb, int je)
  {
    int mrows(ie - ib), ncols(je - jb);
    DenseMatrix<MEMSPACE> tmp(mrows, ncols);
    for (int j = jb; j < je; ++j) {
      for (int i = ib; i < ie; ++i) { tmp(i - ib, j - jb) = operator()(i, j); }
    }
    return tmp;
  }

  /* ******************************************************************
   * Second level routine: transpose
   ****************************************************************** */
  void Transpose(const DenseMatrix<MEMSPACE>& A)
  {
    int mrowsA = A.NumRows(), ncolsA = A.NumCols();
    AMANZI_ASSERT(n_ == A.m_ && m_ == A.n_);
    //reshape(ncolsA, mrowsA);

    for (int j = 0; j < ncolsA; ++j) {
      for (int i = 0; i < mrowsA; ++i) { operator()(j, i) = A(i, j); }
    }
  }

  int Transpose()
  {
    assert(m_ == n_);
    if (m_ != n_) return 1;

    for (int i = 0; i < m_; ++i) {
      for (int j = i + 1; j < n_; ++j) {
        double tmp = operator()(i, j);
        operator()(i, j) = operator()(j, i);
        operator()(j, i) = tmp;
      }
    }
    return 0;
  }

  /* ******************************************************************
   * Second level routine: inversion
   ****************************************************************** */
  int Inverse()
  {
    assert(m_ == n_);
    if (n_ != m_) return 911;

    int ierr;
    int iwork[n_];
    DGETRF_F77(&n_, &n_, &data_(0), &n_, iwork, &ierr);
    //assert(!ierr);
    if (ierr) return ierr;

    int lwork = n_ * n_;
    double dwork[lwork];
    DGETRI_F77(&n_, &data_(0), &n_, iwork, dwork, &lwork, &ierr);
    //assert(!ierr);
    return ierr;
  }


  /* ******************************************************************
   * Second level routine: inversion of symmetric poitive definite
   ****************************************************************** */
  int InverseSPD()
  {
    assert(n_ == m_);
    if (n_ != m_) return 911;

    int ierr;
    DPOTRF_F77("U", &n_, &data_(0), &n_, &ierr);
    if (ierr) return ierr;

    DPOTRI_F77("U", &n_, &data_(0), &n_, &ierr);
    for (int i = 0; i < n_; i++)
      for (int j = i + 1; j < n_; j++) data_[i * n_ + j] = data_[j * n_ + i];
    return ierr;
  }


  /* ******************************************************************
   * Second level routine: calculates matrix D such that (*this)^T D = 0.
   * The matrix (*this) must have a full rank and have more rows than
   * columns.
   ****************************************************************** */
  DenseMatrix<MEMSPACE> NullSpace()
  {
    // We can treat only one type of rectangular matrix.
    AMANZI_ASSERT(m_ > n_);

    int m = m_, n = m_ - n_;

    // Allocate memory for Lapack routine.
    int ldv(1), lwork, info;
    lwork = std::max(m_ + 3 * n_, 5 * n_);
    double U[m_ * m_], V, S[n_], work[lwork];

    DGESVD_F77("A", "N", &m_, &n_, &data_(0), &m_, S, U, &m, &V, &ldv, work, &lwork, &info);

    AMANZI_ASSERT(!info);

    DenseMatrix<MEMSPACE> D(m, n);
    {
      Kokkos::View<double*, MEMSPACE> data = D.Values();
      int offset = m_ * n_;
      for (int i = 0; i < m * n; i++) data[i] = U[offset + i];
    }
    return D;
  }

 private:
  int m_, n_, mem_;
  Kokkos::View<double*, MEMSPACE> data_;
};

// non-member functions
template <typename MEMSPACE>
KOKKOS_INLINE_FUNCTION bool
operator==(const DenseMatrix<MEMSPACE>& A, const DenseMatrix<MEMSPACE>& B)
{
  if (A.NumRows() != B.NumRows()) return false;
  if (A.NumCols() != B.NumCols()) return false;
  for (int i = 0; i != A.NumRows() * A.NumCols(); ++i)
    if (A.Values()[i] != B.Values()[i]) return false;
  return true;
}


template <typename MEMSPACE>
KOKKOS_INLINE_FUNCTION bool
operator!=(const DenseMatrix<MEMSPACE>& A, const DenseMatrix<MEMSPACE>& B)
{
  return !(A == B);
}


template <typename MEMSPACE>
KOKKOS_INLINE_FUNCTION void
PrintMatrix(const DenseMatrix<MEMSPACE>& A, const char* format = "%12.5f", int mmax = 0)
{
  int m = A.NumRows();
  int n = A.NumCols();

  if (mmax > 0) {
    m = Kokkos::min(mmax, m);
    n = Kokkos::min(mmax, n);
  }

  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) printf(format, A(i, j));
    printf("\n");
  }
  printf("\n");
}

// output
inline std::ostream&
operator<<(std::ostream& os, const DenseMatrix<Kokkos::HostSpace>& A)
{
  os << "NumRows: " << A.NumRows() << " NumCols: " << A.NumCols() << "\n";
  for (int i = 0; i < A.NumRows(); i++) {
    for (int j = 0; j < A.NumCols(); j++) {
      os << std::setw(12) << std::setprecision(12) << A(i, j) << " ";
    }
    os << std::endl;
  }
  return os;
}


} // namespace WhetStone
} // namespace Amanzi

#endif
