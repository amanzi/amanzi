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

#ifndef AMANZI_WEAK_MPC_HH_
#define AMANZI_WEAK_MPC_HH_

#include "PK.hh"
#include "MPC_tmp.hh"

namespace Amanzi {

class MPCWeak : public MPCTmp<PK> {

public:
  MPCWeak(Teuchos::ParameterList& pk_tree,
          const Teuchos::RCP<Teuchos::ParameterList>& global_list,
          const Teuchos::RCP<State>& S,
          const Teuchos::RCP<TreeVector>& soln) :
      MPCTmp<PK>(pk_tree, global_list, S, soln) {};

  // PK methods
  // -- dt is the minimum of the sub pks
  virtual double get_dt();

  // -- advance each sub pk dt.
  virtual bool AdvanceStep(double t_old, double t_new);

private:
  // factory registration
  static RegisteredPKFactory<MPCWeak> reg_;

};
} // close namespace Amanzi

#endif
