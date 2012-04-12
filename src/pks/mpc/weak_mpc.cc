/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Implementation for the derived WeakMPC class.  Provides only the advance()
method missing from MPC.hh.  In weak coupling, we simply loop over the
sub-PKs, calling their advance() methods and returning failure if any fail.

See additional documentation in the base class src/pks/mpc/MPC.hh
------------------------------------------------------------------------- */

#include <string>

#include "errors.hh"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"

#include "weak_mpc.hh"

namespace Amanzi {

RegisteredPKFactory<WeakMPC> WeakMPC::reg_("weak MPC");

WeakMPC::WeakMPC(Teuchos::ParameterList& mpc_plist,
                 const Teuchos::RCP<State>& S, const Teuchos::RCP<TreeVector>& soln) :
    MPC::MPC(mpc_plist, S, soln) {};

// Advance each sub-PK individually.
bool WeakMPC::advance(double dt) {
  bool fail = false;
  for (std::vector< Teuchos::RCP<PK> >::iterator pk = sub_pks_.begin();
       pk != sub_pks_.end(); ++pk) {
    fail = (*pk)->advance(dt);
    if (fail) {
      return fail;
    }
  }
  return fail;
};
} // namespace
