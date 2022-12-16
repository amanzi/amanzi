/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

/*
  This is the mpc_pk component of the Amanzi code.

  Implementation for the derived PK_MPCWeak class.  Provides only the advance()
  method missing from MPC.hh.  In weak coupling, we simply loop over the
  sub-PKs, calling their advance() methods and returning failure if any fail.

  See additional documentation in the base class src/pks/mpc_pk/PK_MPC.hh
*/

#include "PK_MPCWeak.hh"

namespace Amanzi {

RegisteredPKFactory<PK_MPCWeak> PK_MPCWeak::reg_("weak MPC");

} // namespace Amanzi
