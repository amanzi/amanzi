/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

#include "Key.hh"
#include "DataStructuresHelpers.hh"
#include "CompositeVectorFunction.hh"

namespace Amanzi {
namespace Functions {

void
CompositeVectorFunction::Compute(double time, CompositeVector& vec)
{
  for (auto [compname, ps, functor] : *this) {
    if (vec.hasComponent(compname)) {
      // compute directly into vector?
      Impl::computeFunction(*functor, time, *ps, vec);

      // compute into patch then copy to vector?
      // Patch<double> p(ps);
      // Impl::computeFunction(*functor, time, p);
      // patchToCompositeVector(p, compname, vec);
    }
  }
}

} // namespace Functions
} // namespace Amanzi
