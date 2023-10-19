/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: William Perkins, Ethan Coon
*/

#include "MeshFrameworkTraits.hh"

namespace Amanzi {
namespace AmanziMesh {

// -------------------------------------------------------------
// A list of frameworks
// -------------------------------------------------------------
Preference
default_preference()
{
  return Preference{ Framework::MSTK, Framework::MOAB, Framework::SIMPLE };
}


bool
framework_enabled(Framework f)
{
  if (f == Framework::SIMPLE) {
    return true;

#ifdef HAVE_MESH_MSTK
  } else if (f == Framework::MSTK) {
    return true;
#endif

#ifdef HAVE_MESH_MOAB
  } else if (f == Framework::MOAB) {
    return true;
#endif
  }
  return false;
}


Preference
filter_preference(const Preference& pref)
{
  Preference result;
  for (auto p : pref) {
    if (framework_enabled(p)) { result.push_back(p); }
  }
  return result;
}


} // namespace AmanziMesh
} // namespace Amanzi
