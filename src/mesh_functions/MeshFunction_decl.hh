/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon
*/

//! <MISSING_ONELINE_DOCSTRING>

#pragma once

#include <utility>
#include <vector>
#include <string>

#include "Teuchos_RCP.hpp"

#include "Mesh.hh"
#include "MultiFunction.hh"
#include "Patch.hh"

namespace Amanzi {
namespace Functions {

//
// Computes function on a patch.
//
template<class Device=AmanziDefaultDevice>
inline void
computeMeshFunction(const MultiFunction& f, double time, Patch<Device>& p);

//
// Compute functions on a multi-patch.
//
template<class Device=AmanziDefaultDevice>
inline void
computeMeshFunction(const std::vector<Teuchos::RCP<const MultiFunction>>& f,
                    double time, MultiPatch<Device>& mp);

//
// Compute set of functions on CompositeVector
//
template<class Device=AmanziDefaultDevice>
inline void
computeMeshFunction(const std::vector<Teuchos::RCP<const MultiFunction>>& f,
                         double time, const MultiPatchSpace& mp, CompositeVector& cv);


using Spec = std::tuple<Teuchos::Array<std::string>,
                        Teuchos::Array<std::string>,
                        Teuchos::RCP<MultiFunction>>;
             

//
// process a list for regions, components, and functions
//
Spec
processSpecWithFunction(const Teuchos::ParameterList& list,
                        std::string function_name);

//
// parse a list of region-based specs
std::pair<MultiPatchSpace,
          std::vector<Teuchos::RCP<const MultiFunction>>>
processListWithFunction(const Teuchos::ParameterList& list,
                        std::string function_name="function");

} // namespace Functions
} // namespace Amanzi


