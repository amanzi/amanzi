/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

/*
  This is the mpc_pk component of the Amanzi code.

  Class for subcycling a slave step within a master step.
  Assumes that intermediate_time() can be used (i.e. this is not nestable?)

  See additional documentation in the base class src/pks/mpc_pk/PK_MPC.hh
*/

#ifndef AMANZI_PK_MPC_SUBCYCLED_HH_
#define AMANZI_PK_MPC_SUBCYCLED_HH_

#include "PK_MPC.hh"
#include "PK.hh"

namespace Amanzi {

class PK_MPCSubcycled : public PK_MPC<PK> {
 public:
  PK_MPCSubcycled(Teuchos::ParameterList& pk_tree,
                  const Teuchos::RCP<Teuchos::ParameterList>& global_list,
                  const Teuchos::RCP<State>& S,
                  const Teuchos::RCP<TreeVector>& soln);

  // PK methods
  // -- dt is the minimum of the sub pks
  virtual double get_dt();
  virtual void set_dt(double dt){};

  // -- advance each sub pk dt.
  virtual bool AdvanceStep(double t_old, double t_new, bool reinit = false);

 protected:
  int master_;
  int slave_;
  double master_dt_;
  double slave_dt_;
  double min_dt_;

 private:
  // factory registration
  static RegisteredPKFactory<PK_MPCSubcycled> reg_;
};

} // namespace Amanzi

#endif
