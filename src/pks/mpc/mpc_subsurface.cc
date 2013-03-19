/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Interface for the derived MPC for coupling energy and water in the subsurface,
with freezing.

------------------------------------------------------------------------- */

#include "permafrost_model.hh"
#include "mpc_delegate_ewc.hh"
#include "mpc_subsurface.hh"

namespace Amanzi {

RegisteredPKFactory<MPCSubsurface> MPCSubsurface::reg_("subsurface permafrost");


// -- Initialize owned (dependent) variables.
void MPCSubsurface::setup(const Teuchos::Ptr<State>& S) {
  // off-diagonal terms needed by MPCCoupledCells
  plist_.set("conserved quantity A", "water_content");
  plist_.set("conserved quantity B", "energy");
  plist_.set("primary variable A", "pressure");
  plist_.set("primary variable B", "temperature");

  plist_.set("mesh key", "domain");
  MPCCoupledCells::setup(S);

  // select the method used for preconditioning
  std::string precon_string = plist_.get<std::string>("preconditioner type", "picard");
  if (precon_string == "none") {
    precon_type_ = PRECON_NONE;
  } else if (precon_string == "block diagonal") {
    precon_type_ = PRECON_BLOCK_DIAGONAL;
  } else if (precon_string == "picard") {
    precon_type_ = PRECON_PICARD;
  } else if (precon_string == "ewc") {
    precon_type_ = PRECON_EWC;
  } else if (precon_string == "smart ewc") {
    precon_type_ = PRECON_SMART_EWC;
  } else {
    Errors::Message message(std::string("Invalid preconditioner type ")+precon_string);
    Exceptions::amanzi_throw(message);
  }

  // select the method used for nonlinear prediction
  std::string predictor_string = plist_.get<std::string>("predictor type", "none");
  if (predictor_string == "none") {
    predictor_type_ = PREDICTOR_NONE;
  } else if (predictor_string == "heuristic") {
    //    predictor_type_ = PREDICTOR_HEURISTIC;
    //    modify_thaw_to_prev_ = plist_.get<bool>("modify thawing cells to previous temp",true);
  } else if (predictor_string == "ewc") {
    predictor_type_ = PREDICTOR_EWC;
  } else if (predictor_string == "smart ewc") {
    predictor_type_ = PREDICTOR_SMART_EWC;
  } else {
    Errors::Message message(std::string("Invalid predictor type ")+predictor_string);
    Exceptions::amanzi_throw(message);
  }

  // create the EWC delegate if requested.
  if (precon_type_ == PRECON_EWC || precon_type_ == PRECON_SMART_EWC ||
      predictor_type_ == PREDICTOR_EWC || predictor_type_ == PREDICTOR_SMART_EWC) {
    ewc_ = Teuchos::rcp(new MPCDelegateEWC(plist_));
    Teuchos::RCP<PermafrostModel> model = Teuchos::rcp(new PermafrostModel());
    ewc_->set_model(model);
    ewc_->setup(S);
  }


}

void MPCSubsurface::initialize(const Teuchos::Ptr<State>& S) {
  MPCCoupledCells::initialize(S);
  if (ewc_ != Teuchos::null) ewc_->initialize(S);
}


void MPCSubsurface::set_states(const Teuchos::RCP<const State>& S,
        const Teuchos::RCP<State>& S_inter,
        const Teuchos::RCP<State>& S_next) {
  MPCCoupledCells::set_states(S,S_inter,S_next);
  if (ewc_ != Teuchos::null) ewc_->set_states(S,S_inter,S_next);
}

void MPCSubsurface::commit_state(double dt, const Teuchos::RCP<State>& S) {
  MPCCoupledCells::commit_state(dt,S);
  if (ewc_ != Teuchos::null) ewc_->commit_state(dt,S);
}

// update the predictor to be physically consistent
bool MPCSubsurface::modify_predictor(double h, Teuchos::RCP<TreeVector> up) {
  bool modified(false);
  if (predictor_type_ == PREDICTOR_EWC || predictor_type_ == PREDICTOR_SMART_EWC) {
    modified = ewc_->modify_predictor(h, up);
  } else if (predictor_type_ == PREDICTOR_HEURISTIC) {
    modified = modify_predictor_heuristic_(h, up);
  }

  // potentiall update faces
  modified |= MPCCoupledCells::modify_predictor(h, up);
  return modified;
}


// updates the preconditioner
void MPCSubsurface::update_precon(double t, Teuchos::RCP<const TreeVector> up, double h) {
  if (precon_type_ == PRECON_NONE) {
    // nothing to do
  } else if (precon_type_ == PRECON_BLOCK_DIAGONAL) {
    StrongMPC::update_precon(t,up,h);
  } else if (precon_type_ == PRECON_PICARD) {
    MPCCoupledCells::update_precon(t,up,h);
  } else if (precon_type_ == PRECON_EWC) {
    MPCCoupledCells::update_precon(t,up,h);
    ewc_->update_precon(t,up,h);
  } else if (precon_type_ == PRECON_SMART_EWC) {
    MPCCoupledCells::update_precon(t,up,h);
    ewc_->update_precon(t,up,h);
  }
}


// -----------------------------------------------------------------------------
// Wrapper to call the requested preconditioner.
// -----------------------------------------------------------------------------
void MPCSubsurface::precon(Teuchos::RCP<const TreeVector> u,
        Teuchos::RCP<TreeVector> Pu) {
  if (precon_type_ == PRECON_NONE) {
    *Pu = *u;
  } else if (precon_type_ == PRECON_BLOCK_DIAGONAL) {
    StrongMPC::precon(u,Pu);
  } else if (precon_type_ == PRECON_PICARD) {
    MPCCoupledCells::precon(u,Pu);
  } else if (precon_type_ == PRECON_EWC) {
    MPCCoupledCells::precon(u,Pu);
    ewc_->precon(u, Pu);
  } else if (precon_type_ == PRECON_SMART_EWC) {
    MPCCoupledCells::precon(u,Pu);
    ewc_->precon(u, Pu);
  }
}


bool MPCSubsurface::modify_predictor_heuristic_(double h, Teuchos::RCP<TreeVector> up) {
  ASSERT(0);
}


} // namespace
