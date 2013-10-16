/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Interface for the derived MPC for coupling energy and water in the subsurface,
with freezing.

------------------------------------------------------------------------- */
#include "Epetra_FEVbrMatrix.h"
#include "EpetraExt_RowMatrixOut.h"

#include "permafrost_model.hh"
#include "mpc_delegate_ewc.hh"
#include "mpc_subsurface.hh"

#define DEBUG_FLAG 1

namespace Amanzi {

RegisteredPKFactory<MPCSubsurface> MPCSubsurface::reg_("subsurface permafrost");


// -- Initialize owned (dependent) variables.
void MPCSubsurface::setup(const Teuchos::Ptr<State>& S) {
  dumped_ = false;

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
    if (modified) changed_solution();
  } else if (predictor_type_ == PREDICTOR_HEURISTIC) {
    modified = modify_predictor_heuristic_(h, up);
    if (modified) changed_solution();
  }

  // potentially update faces
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

  if (S_next_->cycle() == 437 && !dumped_) {
    Teuchos::RCP<const Epetra_FEVbrMatrix> sc = mfd_preconditioner_->Schur();
    std::stringstream filename_s;
    filename_s << "schur_" << S_next_->cycle() << ".txt";
    EpetraExt::RowMatrixToMatlabFile(filename_s.str().c_str(), *sc);
    dumped_ = true;
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
    mfd_preconditioner_->UpdateConsistentFaceCorrection(*u, Pu.ptr());
  } else if (precon_type_ == PRECON_SMART_EWC) {
    MPCCoupledCells::precon(u,Pu);

    // make sure we can back-calc face corrections that preserve residuals on faces
    Teuchos::RCP<TreeVector> res0 = Teuchos::rcp(new TreeVector(*u));
    res0->PutScalar(0.);
    Teuchos::RCP<TreeVector> Pu_std = Teuchos::rcp(new TreeVector(*Pu));
    *Pu_std = *Pu;

    // call EWC, which does Pu_p <-- Pu_p_std + dPu_p
    ewc_->precon(u, Pu);

    // calculate dPu_lambda from dPu_p
    Pu_std->Update(1.0, *Pu, -1.0);
    mfd_preconditioner_->UpdateConsistentFaceCorrection(*res0, Pu_std.ptr());

    // update Pu_lambda <-- Pu_lambda_std + dPu_lambda
    Pu->SubVector(0)->Data()->ViewComponent("face",false)->Update(1.,
            *Pu_std->SubVector(0)->Data()->ViewComponent("face",false), 1.);
    Pu->SubVector(1)->Data()->ViewComponent("face",false)->Update(1.,
            *Pu_std->SubVector(1)->Data()->ViewComponent("face",false), 1.);

    // write residuals
    if (vo_->os_OK(Teuchos::VERB_HIGH)) *vo_->os() << "EWC Precon Correction:" << std::endl;
    std::vector<std::string> vnames;
    vnames.push_back("  PC*r_p"); vnames.push_back("  PC*r_T"); 
    std::vector< Teuchos::Ptr<const CompositeVector> > vecs;
    vecs.push_back(Pu->SubVector(0)->Data().ptr()); 
    vecs.push_back(Pu->SubVector(1)->Data().ptr()); 
    db_->WriteVectors(vnames, vecs, true);
  }
}


bool MPCSubsurface::modify_predictor_heuristic_(double h, Teuchos::RCP<TreeVector> up) {
  Errors::Message message("MPCSubsurface: heuristic no longer implemented");
  Exceptions::amanzi_throw(message);
  return false;
}

bool MPCSubsurface::modify_correction(double h,
        Teuchos::RCP<const TreeVector> res,
        Teuchos::RCP<const TreeVector> u,
        Teuchos::RCP<TreeVector> du) {

  // write residuals
  if (vo_->os_OK(Teuchos::VERB_HIGH)) *vo_->os() << "NKA'd Correction:" << std::endl;
  std::vector<std::string> vnames;
  vnames.push_back("  NKA*PC*r_p"); vnames.push_back("  NKA*PC*r_T"); 
  std::vector< Teuchos::Ptr<const CompositeVector> > vecs;
  vecs.push_back(du->SubVector(0)->Data().ptr()); 
  vecs.push_back(du->SubVector(1)->Data().ptr()); 
  db_->WriteVectors(vnames, vecs, true);

  return false;
}

} // namespace
