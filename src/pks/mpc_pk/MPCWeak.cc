/*
  This is the mpc_pk component of the Amanzi code. 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon

  Implementation for the derived MPCWeak class.  Provides only the Advance()
  method missing from MPC.hh.  In weak coupling, we simply loop over the
  sub-PKs, calling their advance() methods and returning failure if any fail.

  Simplest form of sequential coupling.

  See additional documentation in the base class src/pks/mpc_pk/MPC_PK.hh
*/

#include "MPCWeak.hh"

namespace Amanzi {

// -----------------------------------------------------------------------------
// Calculate the min of sub PKs timestep sizes.
// -----------------------------------------------------------------------------
double MPCWeak::get_dt() {
  double dt = 1.0e99;
  for (MPC_PK<PK>::SubPKList::iterator pk = sub_pks_.begin();
       pk != sub_pks_.end(); ++pk) {
    dt = std::min<double>(dt, (*pk)->get_dt());
  }
  return dt;
}


// -----------------------------------------------------------------------------
// Advance each sub-PK individually, returning a failure as soon as possible.
// -----------------------------------------------------------------------------
bool MPCWeak::AdvanceStep(double t_old, double t_new) {
  bool fail = false;
  for (MPC_PK<PK>::SubPKList::iterator pk = sub_pks_.begin();
       pk != sub_pks_.end(); ++pk) {
    fail = (*pk)->AdvanceStep(t_old, t_new);
    if (fail) {
      return fail;
    }
  }
  return fail;
}

}  // namespace Amanzi
