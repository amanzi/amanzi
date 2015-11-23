/*
  This is the mimetic discretization component of the Amanzi code. 

  Copyright 2010-20XX held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Version: 2.0
  Release name: naka-to.
  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#ifndef AMANZI_DENSE_VECTOR_HH_
#define AMANZI_DENSE_VECTOR_HH_

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iomanip>

#include "lapack.hh"

namespace Amanzi {
namespace WhetStone {

class DenseVector {
 public:
  DenseVector() {
    m_ = 0;
    data_ = NULL;
  }

  DenseVector(int mrow) {
    m_ = mrow;
    data_ = new double[m_];
  }

  DenseVector(int mrow, double* data) {
    m_ = mrow;
    data_ = new double[m_]; 
    for (int i = 0; i < m_; i++) data_[i] = data[i];
  }

  DenseVector(const DenseVector& B) {
    m_ = B.NumRows();
    data_ = new double[m_];
    const double* dataB = B.Values();
    for (int i = 0; i < m_; i++) data_[i] = dataB[i];
  }

  ~DenseVector() { if (data_ != NULL) { delete[] data_; } }

  // primary members 
  void clear() { for (int i = 0; i < m_; i++) data_[i] = 0.0; } 

  double& operator()(int i) { return data_[i]; }
  const double& operator()(int i) const { return data_[i]; }

  DenseVector& operator=(const DenseVector& B) {
    if (this != &B) {
      if (m_ != B.m_ && B.m_ != 0) {
	if (data_ != NULL) {
	  delete [] data_;
	}
	data_ = new double[B.m_];
      }
      m_ = B.m_;
      const double *b = B.Values();
      for (int i = 0; i < m_; ++i) data_[i] = b[i];
    }
    return (*this);
  }

  void PutScalar(double val) {
    for (int i = 0; i < m_; i++) data_[i] = val;
  }

  int Dot(const DenseVector& B, double* result) {
    if (m_ != B.m_) return -1;

    const double *b = B.Values();
    *result = 0.0;
    for (int i = 0; i < m_; i++) *result += data_[i] * b[i];

    return 0;
  }

  // access
  int NumRows() const { return m_; }
  double* Values() { return data_; }
  const double* Values() const { return data_; }

  // output 
  friend std::ostream& operator << (std::ostream& os, DenseVector& A) {
    for (int i = 0; i < A.NumRows(); i++)
        os << std::setw(12) << std::setprecision(12) << A(i) << " ";
    os << "\n";
    return os;
  }

  // First level routines
  void Norm2(double* result) {
    *result = 0.0;
    for (int i = 0; i < m_; i++) *result += data_[i] * data_[i];
    *result = std::pow(*result, 0.5);
  }
 
 private:
  int m_;
  double* data_;                       
};

}  // namespace WhetStone
}  // namespace Amanzi

#endif
