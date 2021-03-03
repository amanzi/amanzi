/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/*
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL.
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
Preference default_preference() {
  return Preference{Framework::MSTK, Framework::MOAB, Framework::STK, Framework::SIMPLE};
}


bool framework_enabled(Framework f) {
  if (f == Framework::SIMPLE) {
    return true;

#ifdef HAVE_MSTK_MESH
  } else if (f == Framework::MSTK) {
    return true;
#endif

#ifdef HAVE_MOAB_MESH
  } else if (f == Framework::MOAB) {
    return true;
#endif

#ifdef HAVE_STK_MESH
  } else if (f == Framework::STK) {
    return true;
#endif

  }
  return false;
}


Preference filter_preference(const Preference& pref) {
  Preference result;
  for (auto p : pref) {
    if (framework_enabled(p)) {
      result.push_back(p);
    }
  }
  return result;
}


} // namespace AmanziMesh
} // namespace Amanzi
