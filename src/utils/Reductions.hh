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
  Scalar val;
  GO loc;

  KOKKOS_INLINE_FUNCTION
  MinLoc() { init(); }

  KOKKOS_INLINE_FUNCTION
  MinLoc(Scalar s, GO g) : val(s), loc(g) {}

  KOKKOS_INLINE_FUNCTION
  MinLoc(const MinLoc& other) : val(other.val), loc(other.loc) {}

  KOKKOS_INLINE_FUNCTION
  void init() {
    val = 1.0e99;
    loc = -1;
  }
  
  KOKKOS_INLINE_FUNCTION
  MinLoc& operator+=(const MinLoc& other) {
    if (other.val < val) {
      val = other.val;
      loc = other.loc;
    }
    return *this;
  }

  KOKKOS_INLINE_FUNCTION
  void operator += (const volatile MinLoc& other) volatile {
    if (other.val < val) {
      val = other.val;
      loc = other.loc;
    }
  }
};


// template<typename Ordinal, typename Scalar, typename GO>
// struct MinLocArray :
//       public Teuchos::ValueTypeReductionOp<Ordinal, std::pair<Scalar, GO> > {


  
//   void
//   reduce (const Ordinal count,
//           const std::pair<Scalar, GO> inBuffer[],
//           std::pair<Scalar, GO> inoutBuffer[]) const
//   {
//     for (Ordinal ind = 0; ind < count; ++ind) {
//       if (inBuffer[ind].first < inoutBuffer[ind].first) {
//         inoutBuffer[ind] = inBuffer[ind];
//       }
//     }
//   }
// };


template<typename Scalar, typename GO>
struct MaxLoc {
  Scalar val;
  GO loc;

  KOKKOS_INLINE_FUNCTION
  MaxLoc() { init(); }

  KOKKOS_INLINE_FUNCTION
  MaxLoc(Scalar s, GO g) : val(s), loc(g) {}
  
  KOKKOS_INLINE_FUNCTION
  MaxLoc(const MaxLoc& other) : val(other.val), loc(other.loc) {}

  KOKKOS_INLINE_FUNCTION
  void init() {
    val = -1.0e99;
    loc = -1;
  }
  
  KOKKOS_INLINE_FUNCTION
  MaxLoc& operator+=(const MaxLoc& other) {
    if (other.val > val) {
      val = other.val;
      loc = other.loc;
    }
    return *this;
  }

  KOKKOS_INLINE_FUNCTION
  void operator += (const volatile MaxLoc& other) volatile {
    if (other.val > val) {
      val = other.val;
      loc = other.loc;
    }
  }
};

template<typename Scalar, typename GO>
struct ValLoc {
  Scalar val;
  GO loc;
};


template<typename Ordinal, typename Scalar, typename GO>
struct MaxLocArray :
      public Teuchos::ValueTypeReductionOp<Ordinal, std::pair<Scalar, GO> > {

  void
  reduce (const Ordinal count,
          const std::pair<Scalar, GO> inBuffer[],
          std::pair<Scalar, GO> inoutBuffer[]) const
  {
    for (Ordinal ind = 0; ind < count; ++ind) {
      if (inBuffer[ind].first > inoutBuffer[ind].first) {
        inoutBuffer[ind].first = inBuffer[ind].first;
        inoutBuffer[ind].second = inBuffer[ind].second;
      }
    }
  }
};


} // namespace Reductions
} // namespace Amanzi
