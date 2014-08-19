/* -------------------------------------------------------------------------
Amanzi

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Interface for the derived MPCWeak class.  Provides only the Advance() method
missing from MPC.hh.  In weak coupling, we simply loop over the sub-PKs,
calling their Advance() methods and returning failure if any fail.

Simplest form of sequential coupling.

See additional documentation in the base class src/pks/mpc/MPC.hh
------------------------------------------------------------------------- */

#ifndef ARCOS_WEAK_MPC_HH_
#define ARCOS_WEAK_MPC_HH_

#include "PK.hh"
#include "MPC_tmp.hh"

namespace Amanzi {

class MPCWeak : public MPCTmp<PK> {

public:
  MPCWeak(const Teuchos::RCP<Teuchos::ParameterList>& plist,
          Teuchos::ParameterList& FElist,
          const Teuchos::RCP<TreeVector>& soln) :
      MPCTmp<PK>(plist, FElist, soln) {};

  // PK methods
  // -- dt is the minimum of the sub pks
  virtual double get_dt();

  // -- advance each sub pk dt.
  virtual bool Advance(double dt);

private:
  // factory registration
  static RegisteredPKFactory<MPCWeak> reg_;

};
} // close namespace Amanzi

#endif
