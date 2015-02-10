/* -------------------------------------------------------------------------
Amanzi

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Class for subcycling a slave step within a master step.
Assumes that intermediate_time() can be used (i.e. this is not nestable?)

See additional documentation in the base class src/pks/mpc/MPC.hh
------------------------------------------------------------------------- */

#ifndef AMANZI_SUBCYCLED_MPC_HH_
#define AMANZI_SUBCYCLED_MPC_HH_

#include "PK.hh"
#include "MPC_tmp.hh"

namespace Amanzi {

class MPCSubcycled : public MPCTmp<PK> {

public:
  MPCSubcycled(Teuchos::ParameterList& pk_tree,
               const Teuchos::RCP<Teuchos::ParameterList>& global_list,
               const Teuchos::RCP<State>& S,
               const Teuchos::RCP<TreeVector>& soln);

  // PK methods
  // -- dt is the minimum of the sub pks
  virtual double get_dt();
  virtual void set_dt(double dt){};

  // -- advance each sub pk dt.
  virtual bool AdvanceStep(double t_old, double t_new);

protected:
  int master_;
  int slave_;
  double master_dt_;
  double slave_dt_;
  double min_dt_;

private:
  // factory registration
  static RegisteredPKFactory<MPCSubcycled> reg_;

};
} // close namespace Amanzi

#endif
