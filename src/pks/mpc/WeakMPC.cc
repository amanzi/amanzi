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

#include "WeakMPC.hh"

namespace Amanzi {

WeakMPC::WeakMPC(Teuchos::ParameterList& mpc_plist,
                 Teuchos::RCP<State>& S, Teuchos::RCP<TreeVector>& soln) :
    MPC::MPC(mpc_plist, S, soln) {};

// Advance each sub-PK individually.
bool WeakMPC::advance(double dt, Teuchos::RCP<TreeVector>& solution) {
  bool fail = false;
  for (std::vector< Teuchos::RCP<PK> >::iterator pk = sub_pks_.begin();
       pk != sub_pks_.end(); ++pk) {
    Teuchos::RCP<TreeVector> subvec;
    int subvec_not_found = solution->SubVector((*pk)->name(), subvec);
    if (subvec_not_found) {
      Errors::Message message("WeakMPC: vector structure does not match PK structure");
      Exceptions::amanzi_throw(message);
    } else {
      fail = (*pk)->advance(dt, subvec);
      if (fail) {
        return fail;
      }
    }
  }
  return fail;
};
} // namespace
