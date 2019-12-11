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

#include "Teuchos_RCP.hpp"
#include "lapack.hh"

namespace Amanzi {
namespace WhetStone {

class DenseVector {
 public:
  DenseVector() : m_(0), mem_(0), data_(NULL){};

  explicit DenseVector(int mrow) : m_(mrow), mem_(mrow)
  {
    data_ = new double[mem_];
    map_ = Teuchos::rcp(new int(mem_)); 
  }

  DenseVector(int mrow, double* data)
  {
    m_ = mrow;
    mem_ = mrow;
    data_ = new double[mem_];
    for (int i = 0; i < m_; i++) data_[i] = data[i];
    map_ = Teuchos::rcp(new int(mem_)); 
  }

  DenseVector(const DenseVector& B)
  {
    m_ = B.NumRows();
    mem_ = m_;
    data_ = NULL;
    if (m_ > 0) {
      data_ = new double[m_];
      const double* dataB = B.Values();
      for (int i = 0; i < m_; i++) data_[i] = dataB[i];
    }
    map_ = Teuchos::rcp(new int(mem_)); 
  }

  DenseVector(const Teuchos::RCP<const int>& map){
    map_ = map; 
    mem_ = *(map.get()); 
    m_ = mem_; 
    data_ = NULL; 
    if(m_ > 0){
      data_ = new double[mem_]; 
    }
  }

  ~DenseVector()
  {
    if (data_ != NULL) { delete[] data_; }
  }

  // primary members
  // -- smart memory management: preserves data only for vector reduction
  void reshape(int mrow);

  // -- initialization
  DenseVector& operator=(const DenseVector& B);

  void putScalar(double val)
  {
    for (int i = 0; i < m_; i++) data_[i] = val;
  }

  // -- access to components
  double& operator()(int i) { return data_[i]; }
  const double& operator()(int i) const { return data_[i]; }

  // -- dot products
  int dot(const DenseVector& B, double* result)
  {
    if (m_ != B.m_) return -1;

    const double* b = B.Values();
    *result = 0.0;
    for (int i = 0; i < m_; i++) *result += data_[i] * b[i];

    return 0;
  }

  friend double operator*(const DenseVector& A, const DenseVector& B)
  {
    double s = 0.0;
    for (int i = 0; i < A.NumRows(); i++) s += A(i) * B(i);
    return s;
  }

  // -- scalar type behaviour
  DenseVector& operator*=(double val)
  {
    for (int i = 0; i < m_; ++i) data_[i] *= val;
    return *this;
  }

  DenseVector& operator/=(double val)
  {
    for (int i = 0; i < m_; ++i) data_[i] /= val;
    return *this;
  }

  DenseVector& operator-=(double val)
  {
    for (int i = 0; i < m_; ++i) data_[i] -= val;
    return *this;
  }

  DenseVector& operator=(double val)
  {
    if (data_ == NULL) {
      m_ = 1;
      mem_ = 1;
      data_ = new double[mem_];
    }
    for (int i = 0; i < m_; ++i) data_[i] = val;
    return *this;
  }

  // ring algebra
  friend DenseVector operator*(double val, const DenseVector& v)
  {
    DenseVector tmp(v);
    tmp *= val;
    return tmp;
  }

  // -- vector type behaviour (no checks for compatiility)
  DenseVector& operator+=(const DenseVector& v)
  {
    const double* datav = v.Values();
    for (int i = 0; i < m_; ++i) data_[i] += datav[i];
    return *this;
  }

  DenseVector& operator-=(const DenseVector& A)
  {
    const double* dataA = A.Values();
    for (int i = 0; i < m_; ++i) data_[i] -= dataA[i];
    return *this;
  }

  // compatibility members
  // -- for nonlinear solvers: this = sa * A + sthis * this
  DenseVector& update(double sa, const DenseVector& A, double sthis)
  {
    const double* dataA = A.Values();
    for (int i = 0; i < m_; ++i) data_[i] = sa * dataA[i] + sthis * data_[i];
    return *this;
  }

  // -- for nonlinear solvers: this = sa * A + sb * B + sthis * this
  DenseVector& update(double sa, const DenseVector& A, double sb,
                      const DenseVector& B, double sthis)
  {
    const double* dataA = A.Values();
    const double* dataB = B.Values();
    for (int i = 0; i < m_; ++i) {
      data_[i] = sa * dataA[i] + sb * dataB[i] + sthis * data_[i];
    }
    return *this;
  }

  // -- scale
  void scale(double val)
  {
    for (int i = 0; i < m_; ++i) data_[i] *= val;
  }

  // access to private data
  int NumRows() const { return m_; }
  double* Values() { return data_; }
  const double* Values() const { return data_; }

  // output
  friend std::ostream& operator<<(std::ostream& os, const DenseVector& A)
  {
    for (int i = 0; i < A.NumRows(); i++)
      os << std::setw(12) << std::setprecision(12) << A(i) << " ";
    os << "\n";
    return os;
  }

  // First level routines
  double norm2() const
  {
    double result = 0.0;
    for (int i = 0; i < m_; i++) result += data_[i] * data_[i];
    return std::pow(result, 0.5);
  }

  // -- we use 'inf' instead of 'max' for compatibility with solvers
  double normInf() const
  {
    double result = 0.0;
    for (int i = 0; i < m_; ++i) {
      result = std::max(result, std::fabs(data_[i]));
    }
    return result; 
  }

  void SwapRows(int m1, int m2)
  {
    double tmp = data_[m2];
    data_[m2] = data_[m1];
    data_[m1] = tmp;
  }

  const Teuchos::RCP<const int>& getMap() const {
    return map_;
  }

 private:
  Teuchos::RCP<const int> map_;
  int m_, mem_;
  double* data_;
};

} // namespace WhetStone
} // namespace Amanzi

#endif
