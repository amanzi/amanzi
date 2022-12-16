/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: William Perkins, others
*/

#pragma once

#include <string>
#include <vector>

namespace Amanzi {
namespace AmanziMesh {

enum struct Framework { SIMPLE=0, MSTK, MOAB };

using Preference = std::vector<Framework>;


inline
std::string to_string(const Framework framework) {
  switch(framework) {
    case(Framework::MSTK): return "MSTK";
    case(Framework::SIMPLE): return "Simple";
    case(Framework::MOAB): return "MOAB";
    default: return "unknown";
  }
}

Preference default_preference();
bool framework_enabled(Framework f);
Preference filter_preference(const Preference& pref);


} // namespace AmanziMesh
} // namespace Amanzi

