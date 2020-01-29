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

#include <Kokkos_Core.hpp>

namespace Amanzi {
namespace WhetStone {

class DenseVector {
 public:
  KOKKOS_INLINE_FUNCTION DenseVector() : m_(0), mem_(0) {};

  explicit DenseVector(int mrow) : m_(mrow), mem_(mrow)
  {
    Kokkos::resize(data_,mem_);
  }

  KOKKOS_INLINE_FUNCTION 
  DenseVector(Kokkos::View<double*> base, int mrow): m_(mrow), mem_(mrow){
    //assert(base.extent(0)>=mrow); 
    data_ = base; 
  }

  KOKKOS_INLINE_FUNCTION 
  DenseVector(Kokkos::View<double*> base)
      : m_(base.extent(0)),
        mem_(base.extent(0)),
        data_(base) {}
  
  DenseVector(int mrow, double* data)
  {
    m_ = mrow;
    mem_ = mrow;
    Kokkos::resize(data_,mem_); 
    for (int i = 0; i < m_; i++) data_[i] = data[i];
  }

  DenseVector(const DenseVector& other)
  {
    m_ = other.NumRows();
    if (m_ != mem_) {
      mem_ = m_;
      Kokkos::resize(data_, mem_);
      Kokkos::deep_copy(data_, other.data_);
    }
  }


#if 0 
  DenseVector(const Teuchos::RCP<const int>& map){
    //map_ = map; 
    mem_ = *(map.get()); 
    m_ = mem_; 
    data_ = NULL; 
    if(m_ > 0){
      data_ = new double[mem_]; 
    }
  }
#endif 

  // primary members
  // -- smart memory management: preserves data only for vector reduction
  void reshape(int mrow);

  // -- initialization
  KOKKOS_INLINE_FUNCTION
  void assign(const DenseVector& other) {
    if (this != &other) {
      assert(mem_ == other.mem_);
      for(int i = 0 ; i < mem_; ++i){
        data_[i] = other.data_[i]; 
      }
    }
  }
  
  KOKKOS_INLINE_FUNCTION void putScalar(double val)
  {
    for (int i = 0; i < m_; i++) data_[i] = val;
  }

  // -- access to components
  KOKKOS_INLINE_FUNCTION double& operator()(int i) { return data_[i]; }
  KOKKOS_INLINE_FUNCTION const double& operator()(int i) const { return data_[i]; }

  // -- dot products
  KOKKOS_INLINE_FUNCTION int dot(const DenseVector& B, double* result)
  {
    if (m_ != B.m_) return -1;

    const Kokkos::View<double*> b = B.Values();
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
  KOKKOS_INLINE_FUNCTION DenseVector& operator*=(double val)
  {
    for (int i = 0; i < m_; ++i) data_[i] *= val;
    return *this;
  }

  KOKKOS_INLINE_FUNCTION DenseVector& operator/=(double val)
  {
    for (int i = 0; i < m_; ++i) data_[i] /= val;
    return *this;
  }

  KOKKOS_INLINE_FUNCTION DenseVector& operator-=(double val)
  {
    for (int i = 0; i < m_; ++i) data_[i] -= val;
    return *this;
  }

  KOKKOS_INLINE_FUNCTION DenseVector& operator=(double val)
  {
    assert(data_.extent(0) != 0); 
    //if (data_ == NULL) {
    //  m_ = 1;
    //  mem_ = 1;
    //  data_ = new double[mem_];
    //}
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
  KOKKOS_INLINE_FUNCTION DenseVector& operator+=(const DenseVector& v)
  {
    const Kokkos::View<double*> datav = v.Values();
    for (int i = 0; i < m_; ++i) data_[i] += datav[i];
    return *this;
  }

  KOKKOS_INLINE_FUNCTION DenseVector& operator-=(const DenseVector& A)
  {
    const Kokkos::View<double*> dataA = A.Values();
    for (int i = 0; i < m_; ++i) data_[i] -= dataA[i];
    return *this;
  }

  // compatibility members
  // -- for nonlinear solvers: this = sa * A + sthis * this
  KOKKOS_INLINE_FUNCTION DenseVector& update(double sa, const DenseVector& A, double sthis)
  {
    const Kokkos::View<double*> dataA = A.Values();
    for (int i = 0; i < m_; ++i) data_[i] = sa * dataA[i] + sthis * data_[i];
    return *this;
  }

  // -- for nonlinear solvers: this = sa * A + sb * B + sthis * this
  KOKKOS_INLINE_FUNCTION DenseVector& update(double sa, const DenseVector& A, double sb,
                      const DenseVector& B, double sthis)
  {
    const Kokkos::View<double*> dataA = A.Values();
    const Kokkos::View<double*> dataB = B.Values();
    for (int i = 0; i < m_; ++i) {
      data_[i] = sa * dataA[i] + sb * dataB[i] + sthis * data_[i];
    }
    return *this;
  }

  // -- scale
  KOKKOS_INLINE_FUNCTION void scale(double val)
  {
    for (int i = 0; i < m_; ++i) data_[i] *= val;
  }

  // access to private data
  KOKKOS_INLINE_FUNCTION int NumRows() const { return m_; }
  KOKKOS_INLINE_FUNCTION Kokkos::View<double*> Values() { return data_; }
  KOKKOS_INLINE_FUNCTION const Kokkos::View<double*> Values() const { return data_; }
  KOKKOS_INLINE_FUNCTION double* Values_ptr() {return &data_[0];}
  KOKKOS_INLINE_FUNCTION const double* Values_ptr() const { return &data_[0];}
  // output
  friend std::ostream& operator<<(std::ostream& os, const DenseVector& A)
  {
    for (int i = 0; i < A.NumRows(); i++)
      os << std::setw(12) << std::setprecision(12) << A(i) << " ";
    return os;
  }

  // First level routines
  KOKKOS_INLINE_FUNCTION double norm2() const
  {
    double result = 0.0;
    for (int i = 0; i < m_; i++) result += data_[i] * data_[i];
    return pow(result, 0.5);
  }

  // -- we use 'inf' instead of 'max' for compatibility with solvers
  KOKKOS_INLINE_FUNCTION double normInf() const
  {
    double result = 0.0;
    for (int i = 0; i < m_; ++i) {
      if(result < fabs(data_[i]))
        result = fabs(data_[i]);
    }
    return result; 
  }

  KOKKOS_INLINE_FUNCTION void SwapRows(int m1, int m2)
  {
    double tmp = data_[m2];
    data_[m2] = data_[m1];
    data_[m1] = tmp;
  }

  //const Teuchos::RCP<const int>& getMap() const {
    //return Teuchos::rcp(new int(mem_));
  //  return Teuchos::null; 
  //}

  // Default assigment implies view semantics
  DenseVector& operator=(const DenseVector& B) = default;

 private:

  int m_, mem_;
  Kokkos::View<double*> data_;
};

} // namespace WhetStone
} // namespace Amanzi

#endif
