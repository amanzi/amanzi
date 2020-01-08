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

namespace Amanzi {
namespace WhetStone {

//const int WHETSTONE_DATA_ACCESS_COPY = 1;
//const int WHETSTONE_DATA_ACCESS_VIEW = 2;

class DenseMatrix {
 public:
  KOKKOS_INLINE_FUNCTION DenseMatrix()
  {
    m_ = 0;
    n_ = 0;
    mem_ = 0;
    //data_ = NULL;
    //access_ = WHETSTONE_DATA_ACCESS_COPY;
  }
  DenseMatrix(const int& mrow, const int& ncol)
  {
    m_ = mrow;
    n_ = ncol;
    mem_ = m_ * n_;
    Kokkos::resize(data_,mem_); 
    //data_ = new double[mem_];
    //access_ = WHETSTONE_DATA_ACCESS_COPY;
  } // memory is not initialized

  DenseMatrix(int mrow, int ncol, double* data
              );//, int data_access = WHETSTONE_DATA_ACCESS_COPY);
  DenseMatrix(int mrow, int ncol, const Kokkos::View<double*>& data);
  //DenseMatrix(const DenseMatrix& B);
  DenseMatrix(const DenseMatrix& B, int m1, int m2, int n1, int n2);
  //~DenseMatrix(){}
  
  // primary members
  // -- reshape can be applied only to a matrix that owns data
  // -- data are not remapped to the new matrix shape
  void reshape(int mrow, int ncol);

  KOKKOS_INLINE_FUNCTION double& operator()(int i, int j) { return data_(j * m_ + i); }
  KOKKOS_INLINE_FUNCTION const double& operator()(int i, int j) const { return data_(j * m_ + i); }

#if 0
  KOKKOS_INLINE_FUNCTION DenseMatrix& operator=(const DenseMatrix& B)
  {
    if (this != &B) {
      //if (mem_ != B.m_ * B.n_) {
        //if (data_ != NULL) { delete[] data_; }
      //  mem_ = B.m_ * B.n_;
      //  Kokkos::resize(data_,mem_);
        //data_ = new double[mem_];
      //}
      n_ = B.n_;
      m_ = B.m_;
      mem_ = n_*m_; 
      data_ = B.Values(); 
      //Kokkos::deep_copy(data_,B.data_);  
    }
    return (*this);
  }
#endif 

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
    Kokkos::View<double*> dataA = A.Values(); 
    Kokkos::View<double*> dataB = B.Values(); 
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

#if 0 
  /* ******************************************************************
   * Matrix-vector product. The matrix is ordered by columns.
  ****************************************************************** */
  // calculates B = *this * A
  KOKKOS_INLINE_FUNCTION int Multiply(const DenseVector& A, DenseVector& B, bool transpose) const;
    //Kokkos::View<double*> dataA = A.Values(); 
    //Kokkos::View<double*> dataB = B.Values();
    const double* dataA = A.Values();
    double* dataB = B.Values();

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
#endif 

  KOKKOS_INLINE_FUNCTION void putScalar(double val)
  {
    for (int i = 0; i < m_ * n_; i++) data_(i) = val;
  }

  // access: the data are ordered by columns
  KOKKOS_INLINE_FUNCTION int NumRows() const { return m_; }
  KOKKOS_INLINE_FUNCTION int NumCols() const { return n_; }

  KOKKOS_INLINE_FUNCTION Kokkos::View<double*> Values() { return data_; }
  KOKKOS_INLINE_FUNCTION double& Value(int i, int j) { return data_(j * m_ + i); }
  KOKKOS_INLINE_FUNCTION const Kokkos::View<double*> Values() const { return data_; }
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

  // Second level routines
  // -- submatrix in rows [ib, ie) and colums [jb, je)
  DenseMatrix SubMatrix(int ib, int ie, int jb, int je);
  // -- transpose creates new matrix
  void Transpose(const DenseMatrix& A);
  // -- transpose modifies square matrix
  int Transpose();

  // -- inversion is applicable for square matrices only
  int Inverse();
  int InverseSPD();
  int NullSpace(DenseMatrix& D);
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

#if 0 
\TODO transform for Kokkos::View
  /* ******************************************************************
  * Orthonormalize selected matrix columns.
  ****************************************************************** */
  KOKKOS_INLINE_FUNCTIONint OrthonormalizeColumns(int n1, int n2)
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
#endif 

#if 0 
  \TODO transform for Kokkos::View
  /* ******************************************************************
  * Permutation of matrix columns.
  ****************************************************************** */
  KOKKOS_INLINE_FUNCTION void SwapColumns(int n1, int n2)
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
#endif 

 private:
  int m_, n_, mem_; 
  //access_;
  Kokkos::View<double*> data_; 
  //double* data_;
};

// non-member functions
KOKKOS_INLINE_FUNCTION bool
operator==(const DenseMatrix& A, const DenseMatrix& B)
{
  if (A.NumRows() != B.NumRows()) return false;
  if (A.NumCols() != B.NumCols()) return false;
  for (int i = 0; i != A.NumRows() * A.NumCols(); ++i)
    if (A.Values()[i] != B.Values()[i]) return false;
  return true;
}


KOKKOS_INLINE_FUNCTION bool
operator!=(const DenseMatrix& A, const DenseMatrix& B)
{
  return !(A == B);
}


KOKKOS_INLINE_FUNCTION void
PrintMatrix(const DenseMatrix& A, const char* format = "%12.5f")
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
