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

// this needs to be defined somewhere, even if we don't register it
template<> const std::string MPCStrong<PK_PhysicalBDF_Default>::pk_type_ = "invalid";


template<> const std::string MPCStrong<PK_BDF_Default>::pk_type_ = "strong mpc";
template<> REGISTER_PK(MPCStrong<PK_BDF_Default>);

REGISTER_PK(MPCWeak);

} // namespace Amanzi
