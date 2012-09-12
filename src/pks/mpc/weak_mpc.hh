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

#include "PK.hh"
#include "mpc.hh"

namespace Amanzi {

class WeakMPC : public MPC<PK> {

public:
  WeakMPC(Teuchos::ParameterList& plist,
          const Teuchos::RCP<TreeVector>& soln) :
      PKDefaultBase(plist,soln),
      MPC<PK>(plist, soln) {};

  // PK methods
  // -- dt is the minimum of the sub pks
  virtual double get_dt();

  // -- advance each sub pk dt.
  virtual bool advance(double dt);

private:
  // factory registration
  static RegisteredPKFactory<WeakMPC> reg_;


};
} // close namespace Amanzi

#endif
