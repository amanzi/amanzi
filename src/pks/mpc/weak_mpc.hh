/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Interface for the derived WeakMPC class.  Provides only the advance() method
missing from MPC.hh.  In weak coupling, we simply loop over the sub-PKs,
calling their advance() methods and returning failure if any fail.

See additional documentation in the base class src/pks/mpc/MPC.hh
------------------------------------------------------------------------- */

#ifndef PKS_MPC_WEAKMPC_HH_
#define PKS_MPC_WEAKMPC_HH_

#include <vector>

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "state.hh"
#include "tree_vector.hh"
#include "mpc.hh"

namespace Amanzi {

class WeakMPC : public MPC {

public:
  WeakMPC(Teuchos::ParameterList& mpc_plist,
          const Teuchos::RCP<State>& S, const Teuchos::RCP<TreeVector>& soln);

  virtual bool advance(double dt);

  // factory registration
  static RegisteredPKFactory<WeakMPC> reg_;

};
} // close namespace Amanzi

#endif
