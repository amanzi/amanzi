/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      William Perkins, others
*/

//! <MISSING_ONELINE_DOCSTRING>

#ifndef AMANZI_MESH_FRAMEWORK_HH_
#define AMANZI_MESH_FRAMEWORK_HH_

#include <string>
#include <vector>
#include <map>

namespace Amanzi {
namespace AmanziMesh {

enum struct Framework { SIMPLE = 0, MSTK, MOAB, STK };

using Preference = std::vector<Framework>;

static const std::map<Framework, std::string> framework_names = {
  { Framework::MSTK, std::string("MSTK") },
  { Framework::MOAB, std::string("MOAB") },
  { Framework::STK, std::string("stk:mesh") },
  { Framework::SIMPLE, std::string("Simple") }
};


Preference
default_preference();
bool
framework_enabled(Framework f);
Preference
filter_preference(const Preference& pref);


} // namespace AmanziMesh
} // namespace Amanzi

#endif
