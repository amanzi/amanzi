/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

#pragma once

#include "PK_PhysicalBDF_Default.hh"
#include "MPCWeak.hh"
#include "MPCStrong.hh"
#include "PKFactory.hh"

namespace Amanzi {

template<> REGISTER_PK(MPCStrong<PK_BDF_Default>);

REGISTER_PK(MPCWeak);

} // namespace Amanzi
