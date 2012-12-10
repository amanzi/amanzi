/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Interface for the derived MPC for flow and energy.  This couples using a
block-diagonal coupler.
------------------------------------------------------------------------- */

#ifndef MPC_FROZEN_PREC_COUPLED_FLOW_ENERGY_HH_
#define MPC_FROZEN_PREC_COUPLED_FLOW_ENERGY_HH_

#include "mpc_prec_coupled_flow_energy.hh"

namespace Amanzi {
class MPCFrozenCoupledFlowEnergy : public MPCCoupledFlowEnergy {
 public:
  MPCFrozenCoupledFlowEnergy(Teuchos::ParameterList& plist,
                             const Teuchos::RCP<TreeVector>& soln) :
      PKDefaultBase(plist, soln),
      MPCCoupledFlowEnergy(plist, soln) {}

  // Virtual destructor
  virtual ~MPCFrozenCoupledFlowEnergy() {}

  virtual void fun(double t_old, double t_new, Teuchos::RCP<TreeVector> u_old,
                   Teuchos::RCP<TreeVector> u_new, Teuchos::RCP<TreeVector> g);

  // update the predictor to be physically consistent
  virtual bool modify_predictor(double h, Teuchos::RCP<TreeVector> up);
  virtual bool modify_predictor_temp(double h, Teuchos::RCP<TreeVector> up);

  // -- Initialize owned (dependent) variables.
  virtual void initialize(const Teuchos::Ptr<State>& S);

  virtual bool is_admissible(Teuchos::RCP<const TreeVector> up);

 private:
  // factory registration
  static RegisteredPKFactory<MPCFrozenCoupledFlowEnergy> reg_;

  double the_res_norm_;
  bool modify_thaw_to_prev_;
};
} // namespace

#endif
