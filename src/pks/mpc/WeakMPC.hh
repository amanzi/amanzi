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

#include "State.hh"
#include "TreeVector.hh"

#include "MPC.hh"

namespace Amanzi {

class WeakMPC : public MPC {

public:
  WeakMPC(Teuchos::ParameterList& mpc_plist,
          Teuchos::RCP<State>& S, Teuchos::RCP<TreeVector>& soln);

  ~WeakMPC() {};

  bool advance(double dt, Teuchos::RCP<TreeVector> &solution);
};
} // close namespace Amanzi

#endif
