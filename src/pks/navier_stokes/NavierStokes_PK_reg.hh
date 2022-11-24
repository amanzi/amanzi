/*
  Navier Stokes PK
 
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  NavierStokes PK registration.
*/

#include "NavierStokes_PK.hh"

namespace Amanzi {
namespace NavierStokes {

RegisteredPKFactory<NavierStokes_PK> NavierStokes_PK::reg_("navier stokes");

} // namespace NavierStokes
} // namespace Amanzi
