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

template<typename Ordinal, typename Scalar, typename GO>
class MinLoc :
      public Teuchos::ValueTypeReductionOp<Ordinal, std::pair<Scalar, GO> > {
 public:
  void
  reduce (const Ordinal count,
          const std::pair<Scalar, GO> inBuffer[],
          std::pair<Scalar, GO> inoutBuffer[]) const
  {
    for (Ordinal ind = 0; ind < count; ++ind) {
      const std::pair<Scalar, GO>& in = inBuffer[ind];
      std::pair<Scalar, GO>& inout = inoutBuffer[ind];
      if (in.first < inout.first) {
        inout.first = in.first;
        inout.second = in.second;
      } else if (in.first > inout.first) {
        // Don't need to do anything; inout has the values.
      } else { // equal, or at least one is NaN.
        inout.first = in.first;
        inout.second = std::min (in.second, inout.second);
      }
    }
  }
};

template<typename Ordinal, typename Scalar, typename GO>
class MaxLoc :
      public Teuchos::ValueTypeReductionOp<Ordinal, std::pair<Scalar, GO> > {
 public:
  void
  reduce (const Ordinal count,
          const std::pair<Scalar, GO> inBuffer[],
          std::pair<Scalar, GO> inoutBuffer[]) const
  {
    for (Ordinal ind = 0; ind < count; ++ind) {
      const std::pair<Scalar, GO>& in = inBuffer[ind];
      std::pair<Scalar, GO>& inout = inoutBuffer[ind];
      if (in.first > inout.first) {
        inout.first = in.first;
        inout.second = in.second;
      } else if (in.first < inout.first) {
        // Don't need to do anything; inout has the values.
      } else { // equal, or at least one is NaN.
        inout.first = in.first;
        inout.second = std::max (in.second, inout.second);
      }
    }
  }
};

} // namespace Reductions
} // namespace Amanzi
