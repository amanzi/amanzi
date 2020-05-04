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


template<class MEMSPACE = DefaultHostMemorySpace>
class DenseMatrix {
 public:
  KOKKOS_INLINE_FUNCTION DenseMatrix()
  {
    m_ = 0;
    n_ = 0;
    mem_ = 0;
  }
  DenseMatrix(const int& mrow, const int& ncol)
  {
    m_ = mrow;
    n_ = ncol;
    mem_ = m_ * n_;
    Kokkos::resize(data_,mem_); 
  } // memory is not initialized

  KOKKOS_INLINE_FUNCTION
  DenseMatrix(Kokkos::View<double*,MEMSPACE> data, const int& mrow, const int& ncol){
    m_ = mrow;
    n_ = ncol;
    mem_ = m_ * n_;
    data_ = data; 
  }

  DenseMatrix(const DenseMatrix& other)
  {
    m_ = other.m_;
    n_ = other.n_;
    if (mem_ != other.mem_) {
      mem_ = other.mem_;
      Kokkos::resize(data_, mem_);
      Kokkos::deep_copy(data_, other.data_);
    }
  }

  /* ******************************************************************
  * No memory check is performed: invalid read is possible.
  ****************************************************************** */
  DenseMatrix(int mrow, int ncol, double* data)
  {
    m_ = mrow;
    n_ = ncol;
    mem_ = m_ * n_;
    Kokkos::resize(data_,mem_); 
    for (int i = 0; i < mem_; i++) data_[i] = data[i];
  }

  DenseMatrix(int mrow, int ncol,
    const Kokkos::View<double*, MEMSPACE>& data)
  {
    m_ = mrow;
    n_ = ncol;
    mem_ = m_ * n_;
    Kokkos::resize(data_,mem_); 
    Kokkos::deep_copy(data_,data); 
  }

  /* ******************************************************************
  * Copy constructor creates a new matrix from the submatrix of B in
  * rows m1 to m2-1 and columns n1 to n2-1.
  ****************************************************************** */
  DenseMatrix(const DenseMatrix& B, int m1, int m2, int n1, int n2)
  {
    m_ = m2 - m1;
    n_ = n2 - n1;
    //access_ = WHETSTONE_DATA_ACCESS_COPY;

    mem_ = m_ * n_;
    Kokkos::resize(data_,mem_); 
    //data_ = new double[mem_];

    int mB = B.NumRows();
    int nB = B.NumCols();
    //const double* dataB = B.Values();

    for (int j = n1; j < n2; ++j) {
      //const double* tmpB = dataB + j * mB + m1;
      //double* tmpA = data_ + (j - n1) * m_;
      for (int i = 0; i < m_; ++i) {
        data_((j-n1)*m_+i) = j*mB+m1+i;
        //*tmpA = *tmpB;
        //tmpA++;
        //tmpB++;
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
    Kokkos::resize(data_,mem_); 
  }
  
  KOKKOS_INLINE_FUNCTION 
  void assign(const DenseMatrix& other) {
    if (this != &other) {
      assert(n_ == other.n_ && m_ == other.m_);
      for(int i = 0 ; i < n_*m_; ++i){
        data_(i) = other.data_(i); 
      }
    }
  }

  KOKKOS_INLINE_FUNCTION double& operator()(int i, int j) { return data_(j * m_ + i); }
  KOKKOS_INLINE_FUNCTION const double& operator()(int i, int j) const { return data_(j * m_ + i); }


  KOKKOS_INLINE_FUNCTION DenseMatrix& operator=(double val)
  {
    for (int i = 0; i < m_ * n_; i++) data_(i) = val;
    return *this;
  }

  KOKKOS_INLINE_FUNCTION DenseMatrix& operator*=(double val)
  {
    for (int i = 0; i < m_ * n_; i++) data_(i) *= val;
    return *this;
  }

  KOKKOS_INLINE_FUNCTION DenseMatrix& operator/=(double val)
  {
    for (int i = 0; i < m_ * n_; i++) data_(i) /= val;
    return *this;
  }

  KOKKOS_INLINE_FUNCTION DenseMatrix& operator+=(const DenseMatrix& A)
  {
    for (int i = 0; i < m_ * n_; i++) data_(i) += A.data_[i];
    return *this;
  }

  KOKKOS_INLINE_FUNCTION DenseMatrix& operator-=(const DenseMatrix& A)
  {
    for (int i = 0; i < m_ * n_; i++) data_(i) -= A.data_[i];
    return *this;
  }

  /* ******************************************************************
  * Matrix-matrix product. The matrices are ordered by columns.
  ****************************************************************** */

  // calculates either A * B to A^T * B
  KOKKOS_INLINE_FUNCTION int Multiply(
    const DenseMatrix& A, const DenseMatrix& B, bool transposeA)
  {
    Kokkos::View<double*,MEMSPACE> dataA = A.Values(); 
    Kokkos::View<double*,MEMSPACE> dataB = B.Values(); 
    //const double* dataA = A.Values();
    //const double* dataB = B.Values();

    int mrowsA = A.NumRows(), ncolsA = A.NumCols();
    int mrowsB = B.NumRows(), ncolsB = B.NumCols();

    if (!transposeA) {
      if (ncolsA != mrowsB || m_ != mrowsA || n_ != ncolsB) return 1;

      for (int i = 0; i < m_; i++) {
        //const double* tmpB = dataB;
        for (int j = 0; j < n_; j++) {
          //const double* tmpA = dataA + i;
          double s(0.0);
          for (int k = 0; k < mrowsB; k++) {
            s += dataA(i+k*m_) * dataB(j*mrowsB+k);
            //s += (*tmpA) * (*tmpB);
            //tmpA += m_;
            //tmpB++;
          }
          (*this)(i, j) = s;
        }
      }
    } else {
      if (mrowsA != mrowsB || m_ != ncolsA || n_ != ncolsB) return 1;

      for (int i = 0; i < m_; i++) {
        //const double* tmpB = dataB;
        for (int j = 0; j < n_; j++) {
          //const double* tmpA = dataA + i * mrowsA;
          double s(0.0);
          for (int k = 0; k < mrowsB; k++) {
            s += dataA(i*mrowsA+k) * dataB(j*mrowsB+k); 
            //s += (*tmpA) * (*tmpB);
            //tmpA++;
            //tmpB++;
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
  template<class DV_MEMSPACE>
  KOKKOS_INLINE_FUNCTION 
  int Multiply(const DenseVector<DV_MEMSPACE>& A, DenseVector<DV_MEMSPACE>& B, bool transpose) const
  {
    auto dataA = A.Values();
    auto dataB = B.Values();

    int mrowsA = A.NumRows();
    int mrowsB = B.NumRows();

    if (!transpose) {
      if (n_ != mrowsA || m_ != mrowsB) return 1;

      for (int i = 0; i < m_; i++) {
        //const double* tmpM = data_ + i;
        double s(0.0);
        for (int j = 0; j < n_; j++) {
          s += data_(i*n_+j) * dataA[j]; 
          //s += (*tmpM) * dataA[j];
          //tmpM += m_;
        }
        dataB[i] = s;
      }
    } else {
      if (m_ != mrowsA || n_ != mrowsB) return 1;

      //const double* tmpM = data_;
      for (int i = 0; i < n_; i++) {
        double s(0.0);
        for (int j = 0; j < m_; j++) {
          s += data_(i*m_+j)*dataA[j]; 
          //s += (*tmpM) * dataA[j];
          //tmpM++;
        }
        dataB[i] = s;
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

  KOKKOS_INLINE_FUNCTION Kokkos::View<double*,MEMSPACE> Values() { return data_; }
  KOKKOS_INLINE_FUNCTION double& Value(int i, int j) { return data_(j * m_ + i); }
  KOKKOS_INLINE_FUNCTION const Kokkos::View<double*,MEMSPACE> Values() const { return data_; }
  KOKKOS_INLINE_FUNCTION const double& Value(int i, int j) const { return data_(j * m_ + i); }
  KOKKOS_INLINE_FUNCTION double* Values_ptr() { return &data_(0);}
  KOKKOS_INLINE_FUNCTION const double* Values_ptr() const {return &data_(0); }

  // output
  friend std::ostream& operator<<(std::ostream& os, const DenseMatrix& A)
  {
    os<<"NumRows: "<<A.NumRows()<<" NumCols: "<<A.NumCols()<<"\n";
    for (int i = 0; i < A.NumRows(); i++) {
      for (int j = 0; j < A.NumCols(); j++) {
        os << std::setw(12) << std::setprecision(12) << A(i, j) << " ";
      }
      os << "\n";
    }
    return os;
  }

  /* ******************************************************************
  * Trace operator is extended to rectangular matrices.
  ****************************************************************** */
  KOKKOS_INLINE_FUNCTION double Trace()
  {
    double s(0.0);
    for (int i = 0; i < m_ * n_; i += m_ + 1) s += data_[i];
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
  KOKKOS_INLINE_FUNCTION void MaxRowValue(
    int irow, int jmin, int jmax, int* j, double* value)
  {
    //double* data = data_ + jmin * m_ + irow;
    //*j = jmin;
    //*value = *data;
    *j = jmin; 
    *value = data_(jmin*m_+irow); 
    
    for (int k = jmin + 1; k < jmax + 1; k++) {
      //data += m_; 
      double v = data_(jmin*m_+irow+k*m_);
      if(v > *value){
      //if (*data > *value) {
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
  KOKKOS_INLINE_FUNCTION void MaxRowMagnitude(
    int irow, int jmin, int jmax, int* j, double* value)
  {
    //double* data = data_ + jmin * m_ + irow;
    *j = jmin;
    *value = fabs(data_(jmin*m_+irow));

    for (int k = jmin + 1; k < jmax + 1; k++) {
      //data += m_;
      double v = fabs(data_(jmin*m_+irow+k*m_));
      if (v > *value) {
        *j = k;
        *value = v;
      }
    }
  }

  KOKKOS_INLINE_FUNCTION double NormInf() const
  {
    double a = 0.0;
    for (int i = 0; i < m_ * n_; i++) 
      if(a < data_[i])
        a = data_[i];
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
  KOKKOS_INLINE_FUNCTION double Det(){
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
      //double* tmp1 = data_ + i * m_;

      for (int j = n1; j < i; ++j) {
        //double* tmp2 = data_ + j * m_;

        sum = 0.0;
        for (int k = 0; k < m_; ++k) sum += data_(i*m_+k)*data_(j*m_+k);//tmp1[k] * tmp2[k];
        for (int k = 0; k < m_; ++k) data_(i*m_+k) -= sum * data_(j*m_+k);//tmp1[k] -= sum * tmp2[k];
      }

      sum = 0.0;
      for (int k = 0; k < m_; ++k) sum += data_(i*m_+k)*data_(i*m_+k);//tmp1[k] * tmp1[k];
      if (sum == 0.0) return -1;

      sum = std::pow(1.0 / sum, 0.5);
      for (int k = 0; k < m_; ++k) data_(i*m_+k) *= sum;//tmp1[k] *= sum;
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
      //double* row1 = data_ + n1 * m_;
      //double* row2 = data_ + n2 * m_;

      for (int i = 0; i < m_; i++) {
        tmp = data_(n1*m_+i);//row1[i];
        data_(n1*m_+i) = data_(n2*m_+i);
        data_(n2*m_+i) = tmp;  
        //row1[i] = row2[i];
        //row2[i] = tmp;
      }
    }
  }

  /* ******************************************************************
  * Second level routine: submatrix in rows [ib,ie) and columns [jb,je)
  ****************************************************************** */
  DenseMatrix SubMatrix(int ib, int ie, int jb, int je)
  {
    int mrows(ie - ib), ncols(je - jb);
    DenseMatrix tmp(mrows, ncols);
    Kokkos::View<double*,MEMSPACE> dataB = tmp.Values(); 
    for (int j = jb; j < je; ++j) {
      for (int i = ib; i < ie; ++i) {
        dataB(j*ie+i) = data_(j*m_+ib+i); 
      }
    }
    return tmp;
  }

  /* ******************************************************************
  * Second level routine: transpose
  ****************************************************************** */
  void Transpose(const DenseMatrix& A)
  {
    const Kokkos::View<double*,MEMSPACE> dataA = A.Values(); 
    int mrowsA = A.NumRows(), ncolsA = A.NumCols();

    reshape(ncolsA, mrowsA);

    for (int j = 0; j < ncolsA; ++j) {
      for (int i = 0; i < mrowsA; ++i) {
        data_(i*ncolsA+j) = dataA(j*mrowsA+i); 
      }
    }
  }

  int Transpose()
  {
    if (m_ != n_) return 1;

    for (int i = 0; i < m_; ++i) {
      for (int j = i + 1; j < n_; ++j) {
        double& aij = data_(i * n_ + j);
        double& aji = data_(j * n_ + i);

        double tmp = aij;
        aij = aji;
        aji = tmp;
      }
    }
    return 0;
  }

  /* ******************************************************************
  * Second level routine: inversion
  ****************************************************************** */
  int Inverse()
  {
    if (n_ != m_) return 911;

    int ierr;
    int iwork[n_];
    DGETRF_F77(&n_, &n_, &data_(0), &n_, iwork, &ierr);
    if (ierr) return ierr;

    int lwork = n_ * n_;
    double dwork[lwork];
    DGETRI_F77(&n_, &data_(0), &n_, iwork, dwork, &lwork, &ierr);
    return ierr;
  }


  /* ******************************************************************
  * Second level routine: inversion of symmetric poitive definite
  ****************************************************************** */
  int InverseSPD()
  {
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
  int NullSpace(DenseMatrix& D)
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
      "A", "N", &m_, &n_, &data_(0), &m_, S, U, &m, &V, &ldv, work, &lwork, &info);

    if (info != 0) return WHETSTONE_ELEMENTAL_MATRIX_FAILED;

    Kokkos::View<double*,MEMSPACE> data = D.Values();
    int offset = m_ * n_;
    for (int i = 0; i < m * n; i++) data[i] = U[offset + i];

    return 0;
  }



  // Default assigment implies view semantics
  DenseMatrix& operator=(const DenseMatrix& B) = default;

 private:

  int m_, n_, mem_; 
  //access_;
  Kokkos::View<double*,MEMSPACE> data_; 
  //double* data_;
};

// non-member functions
template<typename MEMSPACE> 
KOKKOS_INLINE_FUNCTION bool
operator==(const DenseMatrix<MEMSPACE>& A, const DenseMatrix<MEMSPACE>& B)
{
  if (A.NumRows() != B.NumRows()) return false;
  if (A.NumCols() != B.NumCols()) return false;
  for (int i = 0; i != A.NumRows() * A.NumCols(); ++i)
    if (A.Values()[i] != B.Values()[i]) return false;
  return true;
}


template<typename MEMSPACE> 
KOKKOS_INLINE_FUNCTION bool
operator!=(const DenseMatrix<MEMSPACE>& A, const DenseMatrix<MEMSPACE>& B)
{
  return !(A == B);
}

template<typename MEMSPACE> 
KOKKOS_INLINE_FUNCTION void
PrintMatrix(const DenseMatrix<MEMSPACE>& A, const char* format = "%12.5f")
{
  int m = A.NumRows();
  int n = A.NumCols();

  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) printf(format, A(i, j));
    printf("\n");
  }
  printf("\n");
}
} // namespace WhetStone
} // namespace Amanzi

#endif
