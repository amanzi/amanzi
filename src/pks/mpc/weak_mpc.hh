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
<<<<<<< HEAD
  WeakMPC(Teuchos::Ptr<State> S,const Teuchos::RCP<Teuchos::ParameterList>& plist,
          Teuchos::ParameterList& FElist,
          const Teuchos::RCP<TreeVector>& soln) :
    PKDefaultBase(S, plist, FElist, soln),
    MPC<PK>(S,plist, FElist, soln) {};
=======
  WeakMPC(const Teuchos::RCP<Teuchos::ParameterList>& plist,
          Teuchos::ParameterList& FElist,
          const Teuchos::RCP<TreeVector>& soln) :
      PKDefaultBase(plist, FElist, soln),
      MPC<PK>(plist, FElist, soln) {};
>>>>>>> 3712d1ddeb1cfe9f074d84ba39b930e7f970357e

  // Virtual destructor
  virtual ~WeakMPC() {}

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
