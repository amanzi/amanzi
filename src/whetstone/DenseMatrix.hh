/*
 This is the mimetic discretization component of the Amanzi code. 

 Copyright 2010-20XX held jointly by LANS/LANL, LBNL, and PNNL. 
 Amanzi is released under the three-clause BSD License. 
 The terms of use and "as is" disclaimer for this license are 
 provided in the top-level COPYRIGHT file.

Version: 2.0
Release name: naka-to.
Author: Konstantin Lipnikov (lipnikov@lanl.gov)

Usage: 
*/

#ifndef __DENSE_MATRIX_HH__
#define __DENSE_MATRIX_HH__

#include <cstdlib>
#include <iostream>
#include <iomanip>

#include "lapack.hh"

namespace Amanzi {
namespace WhetStone {

class DenseMatrix {
 public:
  DenseMatrix() { 
    n_ = 0;
    m_ = 0;
    data_ = NULL;
  }
  DenseMatrix(int nrow, int ncol) { Init(nrow, ncol); }
  ~DenseMatrix() { delete[] data_; }

  // primary members 
  void clear() { for (int i = 0; i < n_ * m_; i++) data_[i] = (double)0; } 

  double& operator()(int i, int j) { return data_[j * n_ + i]; }
  const double& operator()(int i, int j) const { return data_[j * n_ + i]; }

  DenseMatrix& operator=(const DenseMatrix& m) {        
    double *a = (*this).Values();
    const double *b = m.Values();

    for (int i = 0; i < n_ * m_; i++) a[i] = b[i];
    return (*this);
  }

  DenseMatrix& operator*=(double val) {
    for (int i = 0; i < n_ * m_; i++) data_[i] *= val;
    return *this;
  }

  void PutScalar(double val) {
    for (int i = 0; i < n_ * m_; i++) data_[i] = val;
  }

  // access
  int NumCols() { return n_; }
  int NumRows() { return m_; }

  double* Values() { return data_; }
  double* Value(int i, int j)  { return data_ + j * n_ + i; } 
  const double* Values() const { return data_; }

  // output 
  friend std::ostream& operator << (std::ostream& os, DenseMatrix& m) {
    for (int i = 0; i < m.NumRows(); i++) {
      for (int j = 0; j < m.NumCols(); j++) {
        os << std::setw(12) << std::setprecision(12) << m(i, j) << " ";
      }
      os << "\n";
    }
    return os;
  }
 
  // second level routines
  int Inverse();

 private:
  void Init(int nrow, int ncol) { 
    n_ = nrow;
    m_ = ncol;
    data_ = new double[n_ * m_]; 
  }

 private:
  int n_, m_;
  double* data_;                       
};

}  // namespace WhetStone
}  // namespace Amanzi

#endif
