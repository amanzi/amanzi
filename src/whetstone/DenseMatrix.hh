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

  //DenseMatrix(int mrow, int ncol, double* data
  //            );//, int data_access = WHETSTONE_DATA_ACCESS_COPY);
  //DenseMatrix(int mrow, int ncol, const Kokkos::View<double*>& data);
  //DenseMatrix(const DenseMatrix& B);
  //DenseMatrix(const DenseMatrix& B, int m1, int m2, int n1, int n2);
  //~DenseMatrix(){}

  // primary members
  // -- reshape can be applied only to a matrix that owns data
  // -- data are not remapped to the new matrix shape
  //void reshape(int mrow, int ncol);

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

#if 0 
  DenseMatrix& operator+=(const DenseMatrix& A)
  {
    for (int i = 0; i < m_ * n_; i++) data_(i) += A.data_[i];
    return *this;
  }

  DenseMatrix& operator-=(const DenseMatrix& A)
  {
    for (int i = 0; i < m_ * n_; i++) data_(i) -= A.data_[i];
    return *this;
  }

  // calculates either A * B to A^T * B
  int Multiply(const DenseMatrix& A, const DenseMatrix& B, bool transposeA);
  // calculates B = *this * A
  int Multiply(const DenseVector& A, DenseVector& B, bool transpose) const;
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
    for (int i = 0; i < A.NumRows(); i++) {
      for (int j = 0; j < A.NumCols(); j++) {
        os << std::setw(12) << std::setprecision(12) << A(i, j) << " ";
      }
      os << "\n";
    }
    return os;
  }

#if 0 
  // First level routines
  // -- trace of a rectangular matrix
  double Trace();

  // -- extrema in rows and columns
  void MaxRowValue(int irow, int* j, double* value)
  {
    MaxRowValue(irow, 0, n_, j, value);
  }
  void MaxRowValue(int irow, int jmin, int jmax, int* j, double* value);

  void MaxRowMagnitude(int irow, int* j, double* value)
  {
    MaxRowMagnitude(irow, 0, n_, j, value);
  }
  void MaxRowMagnitude(int irow, int jmin, int jmax, int* j, double* value);

  double NormInf() const
  {
    double a = 0.0;
    for (int i = 0; i < m_ * n_; i++) a = std::max(a, data_[i]);
    return a;
  }

  double Norm2() const
  {
    double a = 0.0;
    for (int i = 0; i < m_ * n_; i++) a += data_[i] * data_[i];
    return std::sqrt(a);
  }

  void Scale(double value)
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
  double Det(); // limited capabilities

  // -- orthonormalize matrix columns between n1 and n2-1.
  //    Returns 0 is sucessful.
  int OrthonormalizeColumns(int n1, int n2);

  // -- permutations
  void SwapColumns(int n1, int n2);
#endif 
 private:
  int m_, n_, mem_; 
  //access_;
  Kokkos::View<double*> data_; 
  //double* data_;
};

#if 0 
// non-member functions
inline bool
operator==(const DenseMatrix& A, const DenseMatrix& B)
{
  if (A.NumRows() != B.NumRows()) return false;
  if (A.NumCols() != B.NumCols()) return false;
  for (int i = 0; i != A.NumRows() * A.NumCols(); ++i)
    if (A.Values()[i] != B.Values()[i]) return false;
  return true;
}


inline bool
operator!=(const DenseMatrix& A, const DenseMatrix& B)
{
  return !(A == B);
}


inline void
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
#endif 
} // namespace WhetStone
} // namespace Amanzi

#endif
