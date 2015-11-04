/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Interface for the derived MPC for coupling energy and water in the subsurface,
with freezing.

------------------------------------------------------------------------- */

#ifndef MPC_SUBSURFACE_HH_
#define MPC_SUBSURFACE_HH_

#include "TreeOperator.hh"
#include "pk_physical_bdf_base.hh"
#include "strong_mpc.hh"

namespace Amanzi {

class MPCDelegateEWCSubsurface;

namespace Operators {
class OperatorDiffusion;
class OperatorDiffusionWithGravity;
class OperatorAdvection;
class OperatorAccumulation;
class Operator;
class UpwindTotalFlux;
class Upwinding;
}

namespace Flow {
class Richards;
}

class MPCSubsurface : public StrongMPC<PKPhysicalBDFBase> {

 public:

  MPCSubsurface(const Teuchos::RCP<Teuchos::ParameterList>& plist,
                Teuchos::ParameterList& FElist,
                const Teuchos::RCP<TreeVector>& soln) :
      PKDefaultBase(plist, FElist, soln),
      StrongMPC<PKPhysicalBDFBase>(plist, FElist, soln) {
    dump_ = plist->get<bool>("dump preconditioner", false);
  }

  // -- Initialize owned (dependent) variables.
  virtual void setup(const Teuchos::Ptr<State>& S);
  virtual void initialize(const Teuchos::Ptr<State>& S);

  virtual void set_states(const Teuchos::RCP<const State>& S,
                          const Teuchos::RCP<State>& S_inter,
                          const Teuchos::RCP<State>& S_next);

  virtual void commit_state(double dt, const Teuchos::RCP<State>& S);

  // update the predictor to be physically consistent
  virtual bool ModifyPredictor(double h, Teuchos::RCP<const TreeVector> up0,
          Teuchos::RCP<TreeVector> up);

  // updates the preconditioner
  virtual void UpdatePreconditioner(double t, Teuchos::RCP<const TreeVector> up, double h) {
    UpdatePreconditioner(t, up, h, true);
  }

  virtual void UpdatePreconditioner(double t, Teuchos::RCP<const TreeVector> up, double h, bool assemble);
  
  // preconditioner application
  virtual int ApplyPreconditioner(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu);
  virtual AmanziSolvers::FnBaseDefs::ModifyCorrectionResult
      ModifyCorrection(double h, Teuchos::RCP<const TreeVector> res,
                       Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> du);

  Teuchos::RCP<Operators::TreeOperator> preconditioner() { return preconditioner_; }
  
 protected:

  enum PreconditionerType {
    PRECON_NONE = 0,
    PRECON_BLOCK_DIAGONAL = 1,
    PRECON_PICARD = 2,
    PRECON_EWC = 3
  };

  Teuchos::RCP<Operators::TreeOperator> preconditioner_;
  Teuchos::RCP<Operators::TreeOperator> linsolve_preconditioner_;
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;

  // preconditioner methods
  PreconditionerType precon_type_;

  // off-diagonal precon terms
  // equations are given by:
  // dWC/dt + div q = 0
  // dE/dt + div K grad T + div hq = 0
  Teuchos::RCP<Operators::OperatorDiffusionWithGravity> ddivq_dT_;
  Teuchos::RCP<Operators::Upwinding> uw_dkrdT_;

  Teuchos::RCP<Operators::OperatorAccumulation> dWC_dT_;
  Teuchos::RCP<Operators::Operator> dWC_dT_block_;

  Teuchos::RCP<Operators::OperatorDiffusion> ddivKgT_dp_;
  Teuchos::RCP<Operators::OperatorDiffusionWithGravity> ddivhq_dp_;
  Teuchos::RCP<Operators::OperatorDiffusionWithGravity> ddivhq_dT_;
  Teuchos::RCP<Operators::OperatorAccumulation> dE_dp_;
  Teuchos::RCP<Operators::Operator> dE_dp_block_;

  Teuchos::RCP<Operators::UpwindTotalFlux> upwinding_hkr_;
  Teuchos::RCP<Operators::UpwindTotalFlux> upwinding_dhkr_dp_;
  Teuchos::RCP<Operators::UpwindTotalFlux> upwinding_dhkr_dT_;

  // friend sub-pk Richards (need K_, some flags from private data)
  Teuchos::RCP<Flow::Richards> richards_pk_;

  Key domain_name_;
  Key temp_key_;
  Key pres_key_;
  Key e_key_;
  Key wc_key_;
  Key tc_key_;
  Key uw_tc_key_;
  Key kr_key_;
  Key uw_kr_key_;
  Key hkr_key_;
  Key uw_hkr_key_;
  Key mass_flux_key_;
  Key mass_flux_dir_key_;
  Key rho_key_;
  
  // EWC delegate
  Teuchos::RCP<MPCDelegateEWCSubsurface> ewc_;

  // cruft for easier global debugging
  bool dump_;
  Teuchos::RCP<Debugger> db_;
  
private:
  // factory registration
  static RegisteredPKFactory<MPCSubsurface> reg_;

};


} // namespace

#endif


