/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

#pragma once

#include <string>
#include <utility>
#include <vector>

#include "Teuchos_RCP.hpp"
#include "MeshFunction.hh"
#include "CompositeVector.hh"

namespace Amanzi {
namespace Functions {

class CompositeVectorFunction : public MeshFunction {
 public:
  using MeshFunction::MeshFunction;
  void Compute(double time, CompositeVector& vec);

  void FlagToVector(CompositeVector_<int>& flag_vec);
};


} // namespace Functions
} // namespace Amanzi
