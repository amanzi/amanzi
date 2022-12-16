/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

// Dummy_PK registration
#include "Dummy_PK.hh"

namespace Amanzi {

RegisteredPKFactory<Dummy_PK> Dummy_PK::reg_("dummy");

} // namespace Amanzi
