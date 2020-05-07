/*
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/
// Reductions for min-loc and max-loc based on those from
// Teuchos_TimeMonitor.cpp

#pragma once

#include "Teuchos_ReductionOp.hpp"

namespace Amanzi {
namespace Reductions {


template<typename Scalar, typename GO>
struct MinLoc {
  std::pair<Scalar,GO> val;

  KOKKOS_INLINE_FUNCTION
  MinLoc() { init(); }

  KOKKOS_INLINE_FUNCTION
  MinLoc(Scalar s, GO g) {
    val.first = s;
    val.second = g;
  }

  KOKKOS_INLINE_FUNCTION
  MinLoc(const MinLoc& other) : val(other.val) {}

  KOKKOS_INLINE_FUNCTION
  void init() {
    val.first = 1.0e99;
    val.second = -1;
  }
  
  KOKKOS_INLINE_FUNCTION
  MinLoc& operator+=(const MinLoc& other) {
    if (other.val.first < val.first) {
      val = other.val;
    }
    return *this;
  }

  KOKKOS_INLINE_FUNCTION
  void operator += (const volatile MinLoc& other) volatile {
    if (other.val.first < val.first) {
      val = other.val;
    }
  }
};


template<typename Ordinal, typename Scalar, typename GO>
class MinLocArray :
      public Teuchos::ValueTypeReductionOp<Ordinal, std::pair<Scalar, GO> > {
 public:
  void
  reduce (const Ordinal count,
          const std::pair<Scalar, GO> inBuffer[],
          std::pair<Scalar, GO> inoutBuffer[]) const
  {
    for (Ordinal ind = 0; ind < count; ++ind) {
      if (inBuffer[ind].first < inoutBuffer[ind].first) {
        inoutBuffer[ind] = inBuffer[ind];
      }
    }
  }
};


template<typename Scalar, typename GO>
struct MaxLoc {
  std::pair<Scalar,GO> val;

  KOKKOS_INLINE_FUNCTION
  MaxLoc() { init(); }

  KOKKOS_INLINE_FUNCTION
  MaxLoc(Scalar s, GO g) {
    val.first = s;
    val.second = g;
  }
  
  KOKKOS_INLINE_FUNCTION
  MaxLoc(const MaxLoc& other) : val(other.val) {}

  KOKKOS_INLINE_FUNCTION
  void init() {
    val.first = -1.0e99;
    val.second = -1;
  }
  
  KOKKOS_INLINE_FUNCTION
  MaxLoc& operator+=(const MaxLoc& other) {
    if (other.val.first > val.first) {
      val = other.val;
    }
    return *this;
  }

  KOKKOS_INLINE_FUNCTION
  void operator += (const volatile MaxLoc& other) volatile {
    if (other.val.first > val.first) {
      val = other.val;
    }
  }
};


template<typename Ordinal, typename Scalar, typename GO>
class MaxLocArray :
      public Teuchos::ValueTypeReductionOp<Ordinal, std::pair<Scalar, GO> > {
 public:
  void
  reduce (const Ordinal count,
          const std::pair<Scalar, GO> inBuffer[],
          std::pair<Scalar, GO> inoutBuffer[]) const
  {
    for (Ordinal ind = 0; ind < count; ++ind) {
      if (inBuffer[ind].first > inoutBuffer[ind].first) {
        inoutBuffer[ind] = inBuffer[ind];
      }
    }
  }
};


} // namespace Reductions
} // namespace Amanzi
