/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon
*/

//!

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


// Tag type for derivatives in models.
template <int>
struct Deriv {};


} // namespace Amanzi

#endif
