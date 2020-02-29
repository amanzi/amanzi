/*
  WhetStone, Version 2.2
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Replacement of linear vectors. It may go away after code upgrade.
*/

#ifndef AMANZI_DENSE_VECTOR_HH_
#define AMANZI_DENSE_VECTOR_HH_

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <vector>

#include "lapack.hh"

namespace Amanzi {
namespace WhetStone {

class DenseVector {
 public:
  DenseVector() : m_(0), mem_(0), data_(NULL) {};
  explicit DenseVector(int mrow);
  DenseVector(int mrow, double* data);
  DenseVector(const DenseVector& other);
  DenseVector(DenseVector&& other);
  DenseVector(const std::vector<double>& B);

  ~DenseVector() { if (data_ != NULL) { delete[] data_; } }

  // primary members 
  // -- smart memory management that preserves data
  void Reshape(int mrow);
  void Regroup(int stride1, int stride2);

  // -- assignment operators
  DenseVector& operator=(const DenseVector& other);
  DenseVector& operator=(DenseVector&& other) noexcept;

  // -- initialization
  void PutScalar(double val) {
    for (int i = 0; i < m_; i++) data_[i] = val;
  }

  // -- the size of the vector is not changed. For a short
  //    vector v, remaining data are initialized to val
  void PutVector(const DenseVector& v, double val);

  // -- access to components
  double& operator()(int i) { return data_[i]; }
  const double& operator()(int i) const { return data_[i]; }

  // -- dot products
  int Dot(const DenseVector& B, double* result) {
    if (m_ != B.m_) return -1;

    const double *b = B.Values();
    *result = 0.0;
    for (int i = 0; i < m_; i++) *result += data_[i] * b[i];

    return 0;
  }

  friend double operator*(const DenseVector& A, const DenseVector& B) {
    double s = 0.0; 
    for (int i = 0; i < A.NumRows(); i++ ) s += A(i) * B(i);
    return s;
  }

  // -- scalar type behaviour
  DenseVector& operator*=(double val) {
    for (int i = 0; i < m_; ++i) data_[i] *= val;
    return *this;
  }

  DenseVector& operator/=(double val) {
    for (int i = 0; i < m_; ++i) data_[i] /= val;
    return *this;
  }

  DenseVector& operator-=(double val) {
    for (int i = 0; i < m_; ++i) data_[i] -= val;
    return *this;
  }

  DenseVector& operator=(double val) {
    if (data_ == NULL) {
      m_ = 1;
      mem_ = 1;
      data_ = new double[mem_];
    }
    for (int i = 0; i < m_; ++i) data_[i] = val;
    return *this;
  }

  // ring algebra
  friend DenseVector operator*(double val, const DenseVector& v) {
    DenseVector tmp(v);
    tmp *= val;
    return tmp;
  }

  friend DenseVector operator+(const DenseVector& v1, const DenseVector& v2);
  friend DenseVector operator-(const DenseVector& v1, const DenseVector& v2);

  DenseVector& operator+=(const DenseVector& v) {
    const double* datav = v.Values();  
    for (int i = 0; i < m_; ++i) data_[i] += datav[i];
    return *this;
  }

  DenseVector& operator-=(const DenseVector& A) {
    const double* dataA = A.Values();  
    for (int i = 0; i < m_; ++i) data_[i] -= dataA[i];
    return *this;
  }

  // compatibility members
  // -- for nonlinear solvers: this = sa * A + sthis * this
  DenseVector& Update(double sa, const DenseVector& A, double sthis) {
    const double* dataA = A.Values();  
    for (int i = 0; i < m_; ++i) data_[i] = sa * dataA[i] + sthis * data_[i];
    return *this;
  }

  // -- for nonlinear solvers: this = sa * A + sb * B + sthis * this
  DenseVector& Update(double sa, const DenseVector& A,
                      double sb, const DenseVector& B, double sthis) {
    const double* dataA = A.Values();  
    const double* dataB = B.Values();  
    for (int i = 0; i < m_; ++i) {
      data_[i] = sa * dataA[i] + sb * dataB[i] + sthis * data_[i];
    }
    return *this;
  }

  // -- scale
  void Scale(double val) {
    for (int i = 0; i < m_; ++i) data_[i] *= val;
  }

  // access to private data
  int NumRows() const { return m_; }
  double* Values() { return data_; }
  const double* Values() const { return data_; }

  // output 
  friend std::ostream& operator << (std::ostream& os, const DenseVector& A) {
    for (int i = 0; i < A.NumRows(); i++)
        os << std::setw(12) << std::setprecision(12) << A(i) << " ";
    os << "\n";
    return os;
  }

  // First level routines
  void Norm2(double* result) const {
    *result = 0.0;
    for (int i = 0; i < m_; i++) *result += data_[i] * data_[i];
    *result = std::pow(*result, 0.5);
  }

  // -- we use 'inf' instead of 'max' for compatibility with solvers
  void NormInf(double* result) const {
    *result = 0.0;
    for (int i = 0; i < m_; ++i) {
      *result = std::max(*result, std::fabs(data_[i]));
    }
  }

  void SwapRows(int m1, int m2) {
    double tmp = data_[m2];
    data_[m2] = data_[m1];
    data_[m1] = tmp;
  }
 
 private:
  int m_, mem_;
  double* data_;                       
};


// non-member functions
// find rows with the maximum value
inline int VectorMaxValuePosition(const DenseVector& v) {
  int imax(0);
  const double* data = v.Values();
  for (int i = 1; i < v.NumRows(); ++i) {
    if (data[imax] < data[i]) imax = i;
  }
  return imax;
}

}  // namespace WhetStone
}  // namespace Amanzi

#endif
