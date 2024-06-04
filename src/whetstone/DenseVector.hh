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

  Replacement of linear vectors. It may go away after code upgrade.
*/

#ifndef AMANZI_DENSE_VECTOR_HH_
#define AMANZI_DENSE_VECTOR_HH_

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iomanip>

#include "Kokkos_Core.hpp"
#include "Teuchos_RCP.hpp"
#include "lapack.hh"

#include "dbc.hh"
#include "AmanziTypes.hh"
#include "WhetStoneDefs.hh"

namespace Amanzi {
namespace WhetStone {

template <class MEMSPACE = DefaultHostMemorySpace>
class DenseVector {
 public:
  KOKKOS_INLINE_FUNCTION DenseVector() : m_(0), mem_(0), unmanaged_(false) {};

  explicit DenseVector(int m) : m_(m), mem_(m), unmanaged_(false) { Kokkos::resize(data_, mem_); }

  DenseVector(int mrow, double* data) : m_(mrow), mem_(mrow), unmanaged_(false)
  {
    Kokkos::resize(data_, mem_);
    for (int i = 0; i < m_; i++) data_[i] = data[i];
  }

  template <class MEMSPACE_OTHER>
  DenseVector(const DenseVector<MEMSPACE_OTHER>& other)
    : m_(other.NumRows()), mem_(other.Capacity()), unmanaged_(false)
  {
    if (mem_ != 0) {
      Kokkos::resize(data_, mem_);
      Kokkos::deep_copy(data_, other.Values());
    }
  }

  DenseVector(const DenseVector<MEMSPACE>& other)
    : m_(other.m_), mem_(other.mem_), unmanaged_(other.unmanaged_)
  {
    if (mem_ != 0) {
      Kokkos::resize(data_, mem_);
      Kokkos::deep_copy(data_, other.Values());
    }
  }

#if WHETSTONE_VIEW_SEMANTICS
  // NOTE, it would be nice for this to be truly unmanaged, but would require
  // templating this class which would likely break most of whetstone.  Try
  // this later!  Then resizing would fail to compile.
  KOKKOS_INLINE_FUNCTION
  DenseVector(Kokkos::View<double*, MEMSPACE> base)
    : m_(base.extent(0)), mem_(base.extent(0)), data_(base), unmanaged_(true)
  {}

  // would have to template this as well, and what would the semantics be?
  KOKKOS_INLINE_FUNCTION DenseVector(DenseVector<MEMSPACE>&& other)
    : m_(other.m_), mem_(other.mem_), data_(other.data_), unmanaged_(other.unmanaged_)
  {}

  // Default assigment implies view semantics
  KOKKOS_INLINE_FUNCTION DenseVector<MEMSPACE>& operator=(const DenseVector<MEMSPACE>&& other)
  {
    m_ = other.m_;
    mem_ = other.mem_;
    data_ = other.data_;
    unmanaged_ = other.unmanaged_;
    return *this;
  }

#else
  DenseVector<MEMSPACE>& operator=(const DenseVector<MEMSPACE>& B) = delete;
#endif

  /* ******************************************************************
   * Smart memory management: preserves data only for vector reduction
   ****************************************************************** */
  void reshape(int mrow)
  {
    AMANZI_ASSERT(!unmanaged_);
    AMANZI_ASSERT(mrow >= 0);

    // simple implementation
    if (mem_ != mrow) {
      // note that this copies the first min(mrow, m_) entries into the new view
      Kokkos::resize(data_, mrow);
      m_ = mrow;
      mem_ = mrow;
    }

    // never-shrink implementation
    // m_ = mrow;
    // if (mrow > mem_) {
    //   Kokkos::resize(data_, mrow);
    //   mem_ = mrow;
    // }
  }

  // -- initialization
  KOKKOS_INLINE_FUNCTION
  void assign(const DenseVector<MEMSPACE>& other)
  {
    if (this != &other) {
      if (m_ != other.m_) reshape(other.m_);
      for (int i = 0; i < m_; ++i) { data_[i] = other.data_[i]; }
    }
  }

  KOKKOS_INLINE_FUNCTION DenseVector<MEMSPACE>& operator=(double val)
  {
    assert(data_.extent(0) != 0);
    for (int i = 0; i < m_; ++i) data_[i] = val;
    return *this;
  }

  KOKKOS_INLINE_FUNCTION void putScalar(double val)
  {
    for (int i = 0; i < m_; i++) data_[i] = val;
  }

  KOKKOS_INLINE_FUNCTION void putVector(const DenseVector<MEMSPACE>& v, double val)
  {
    int mmin = std::min(m_, v.NumRows());
    for (int i = 0; i < mmin; i++) data_[i] = v(i);
    for (int i = mmin; i < m_; i++) data_[i] = val;
  }

  // -- access to components
  KOKKOS_INLINE_FUNCTION double& operator()(int i)
  {
#if WHETSTONE_ABC
    AMANZI_ASSERT(i <= m_);
#endif
    return data_[i];
  }
  KOKKOS_INLINE_FUNCTION const double& operator()(int i) const
  {
#if WHETSTONE_ABC
    AMANZI_ASSERT(i <= m_);
#endif
    return data_[i];
  }

  // -- dot products
  KOKKOS_INLINE_FUNCTION int dot(const DenseVector<MEMSPACE>& B, double* result) const
  {
    if (m_ != B.m_) return -1;
    *result = 0.0;
    for (int i = 0; i < m_; i++) *result += data_[i] * B(i);
    return 0;
  }

  // -- scalar type behaviour
  KOKKOS_INLINE_FUNCTION DenseVector<MEMSPACE>& operator*=(double val)
  {
    for (int i = 0; i < m_; ++i) data_[i] *= val;
    return *this;
  }

  KOKKOS_INLINE_FUNCTION DenseVector<MEMSPACE>& operator/=(double val)
  {
    for (int i = 0; i < m_; ++i) data_[i] /= val;
    return *this;
  }

  KOKKOS_INLINE_FUNCTION DenseVector<MEMSPACE>& operator-=(double val)
  {
    for (int i = 0; i < m_; ++i) data_[i] -= val;
    return *this;
  }

  // -- vector type behaviour (no checks for compatiility)
  KOKKOS_INLINE_FUNCTION DenseVector<MEMSPACE>& operator+=(const DenseVector<MEMSPACE>& v)
  {
    for (int i = 0; i < m_; ++i) data_[i] += v(i);
    return *this;
  }

  KOKKOS_INLINE_FUNCTION DenseVector<MEMSPACE>& operator-=(const DenseVector<MEMSPACE>& A)
  {
    for (int i = 0; i < m_; ++i) data_[i] -= A(i);
    return *this;
  }

  // compatibility members
  // -- for nonlinear solvers: this = sa * A + sthis * this
  KOKKOS_INLINE_FUNCTION DenseVector<MEMSPACE>&
  update(double sa, const DenseVector<MEMSPACE>& A, double sthis)
  {
    for (int i = 0; i < m_; ++i) data_[i] = sa * A(i) + sthis * data_[i];
    return *this;
  }

  // -- for nonlinear solvers: this = sa * A + sb * B + sthis * this
  KOKKOS_INLINE_FUNCTION DenseVector<MEMSPACE>& update(double sa,
                                                       const DenseVector<MEMSPACE>& A,
                                                       double sb,
                                                       const DenseVector<MEMSPACE>& B,
                                                       double sthis)
  {
    for (int i = 0; i < m_; ++i) { data_[i] = sa * A(i) + sb * B(i) + sthis * data_[i]; }
    return *this;
  }

  // -- scale
  KOKKOS_INLINE_FUNCTION void scale(double val)
  {
    for (int i = 0; i < m_; ++i) data_[i] *= val;
  }

  // access to private data
  KOKKOS_INLINE_FUNCTION int NumRows() const { return m_; }
  KOKKOS_INLINE_FUNCTION int Capacity() const { return mem_; }
  KOKKOS_INLINE_FUNCTION Kokkos::View<double*, MEMSPACE> Values() { return data_; }
  KOKKOS_INLINE_FUNCTION const Kokkos::View<double*, MEMSPACE> Values() const { return data_; }
  KOKKOS_INLINE_FUNCTION double* Values_ptr() { return &data_[0]; }
  KOKKOS_INLINE_FUNCTION const double* Values_ptr() const { return &data_[0]; }

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
      if (result < fabs(data_[i])) result = fabs(data_[i]);
    }
    return result;
  }

  KOKKOS_INLINE_FUNCTION void SwapRows(int m1, int m2)
  {
    double tmp = data_[m2];
    data_[m2] = data_[m1];
    data_[m1] = tmp;
  }


 private:
  int m_, mem_;
  bool unmanaged_;
  Kokkos::View<double*, MEMSPACE> data_;
};


template <class MEMSPACE>
KOKKOS_INLINE_FUNCTION double
operator*(const DenseVector<MEMSPACE>& A, const DenseVector<MEMSPACE>& B)
{
  double res(0);
  A.dot(B, &res);
  return res;
}

// ring algebra
template <class MEMSPACE>
KOKKOS_INLINE_FUNCTION DenseVector<MEMSPACE>
operator*(double val, const DenseVector<MEMSPACE>& v)
{
  DenseVector<MEMSPACE> tmp(v);
  tmp *= val;
  return tmp;
}

// output
inline std::ostream&
operator<<(std::ostream& os, const DenseVector<Kokkos::HostSpace>& A)
{
  for (int i = 0; i < A.NumRows(); i++) os << std::setw(12) << std::setprecision(12) << A(i) << " ";
  return os;
}


} // namespace WhetStone
} // namespace Amanzi

#endif
