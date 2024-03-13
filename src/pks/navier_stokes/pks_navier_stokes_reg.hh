/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

/*
  Navier Stokes PK

  NavierStokes PK registration.
*/

#include "NavierStokes_PK.hh"

namespace Amanzi {
namespace NavierStokes {

RegisteredPKFactory<NavierStokes_PK> NavierStokes_PK::reg_("navier stokes");

} // namespace NavierStokes
} // namespace Amanzi
