/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
   ATS

   License: see $ATS_DIR/COPYRIGHT
   Author: Ethan Coon

   Some basic typedefs for State and company.
   ------------------------------------------------------------------------- */

#ifndef AMANZI_STATE_DEFS_HH_
#define AMANZI_STATE_DEFS_HH_

#include <set>
#include <string>
#include <vector>

#include "Key.hh"

namespace Amanzi {

// Keys and containers
typedef std::string Units;
typedef bool NullFactory; // placeholder object for no factory required

} // namespace Amanzi

#endif
