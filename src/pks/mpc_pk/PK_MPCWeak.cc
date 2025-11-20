/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

/*
  MPC PK

  Implementation for the derived PK_MPCWeak class.  Provides only the Advance()
  method missing from MPC.hh.  In weak coupling, we simply loop over the
  sub-PKs, calling their advance() methods and returning failure if any fail.

  Simplest form of sequential coupling.

  See additional documentation in the base class src/pks/mpc_pk/PK_MPC.hh
*/

#include "PK_MPCWeak.hh"

namespace Amanzi {

// -----------------------------------------------------------------------------
// Calculate the min of sub PKs timestep sizes.
// -----------------------------------------------------------------------------
double
PK_MPCWeak::get_dt()
{
  double dt = 1.0e99;
  for (PK_MPC<PK>::SubPKList::iterator pk = sub_pks_.begin(); pk != sub_pks_.end(); ++pk) {
    dt = std::min<double>(dt, (*pk)->get_dt());
  }
  return dt;
}


// -----------------------------------------------------------------------------
// Set time step for sub PKs.
// -----------------------------------------------------------------------------
void
PK_MPCWeak::set_dt(double dt)
{
  for (PK_MPC<PK>::SubPKList::iterator pk = sub_pks_.begin(); pk != sub_pks_.end(); ++pk) {
    (*pk)->set_dt(dt);
  }
}


// -----------------------------------------------------------------------------
// Advance each sub-PK individually, returning a failure as soon as possible.
// -----------------------------------------------------------------------------
bool
PK_MPCWeak::AdvanceStep(double t_old, double t_new, bool reinit)
{
  bool fail = false;
  for (int i = 0; i < sub_pks_.size(); ++i) {
    fail = sub_pks_[i]->AdvanceStep(t_old, t_new, reinit);
    if (fail) {
      for (int k = 0; k < i; ++k) {
        sub_pks_[k]->FailStep(t_old, t_new, Tags::DEFAULT);
      }
      return fail;
    }
  }
  return fail;
}

} // namespace Amanzi
