/*
  This is the mpc_pk component of the Amanzi code. 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon

  Implementation for the derived MPCWeak class.  Provides only the advance()
  method missing from MPC.hh.  In weak coupling, we simply loop over the
  sub-PKs, calling their advance() methods and returning failure if any fail.

  See additional documentation in the base class src/pks/mpc_pk/MPC_PK.hh
*/

#include "MPCWeak.hh"

namespace Amanzi {

RegisteredPKFactory<MPCWeak> MPCWeak::reg_("weak MPC");

}  // namespace Amanzi

