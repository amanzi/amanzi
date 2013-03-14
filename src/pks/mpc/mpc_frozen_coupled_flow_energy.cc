/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Derived MPC for flow and energy.  This couples using a block-diagonal coupler.
------------------------------------------------------------------------- */

#include "primary_variable_field_evaluator.hh"
#include "permafrost_model.hh"
#include "wrm_permafrost_evaluator.hh"
#include "pc_ice_evaluator.hh"
#include "pc_liquid_evaluator.hh"
#include "eos_evaluator.hh"
#include "iem_evaluator.hh"
#include "iem_water_vapor_evaluator.hh"
#include "molar_fraction_gas_evaluator.hh"

#include "richards.hh"
#include "two_phase.hh"

#include "mpc_frozen_coupled_flow_energy.hh"

namespace Amanzi {

#define DEBUG_FLAG 1

RegisteredPKFactory<MPCFrozenCoupledFlowEnergy> MPCFrozenCoupledFlowEnergy::reg_("frozen energy-flow preconditioner coupled");


// -------------------------------------------------------------
// Setup data
// -------------------------------------------------------------
void MPCFrozenCoupledFlowEnergy::setup(const Teuchos::Ptr<State>& S) {
  // off-diagonal terms needed by MPCCoupledCells
  plist_.set("conserved quantity A", "water_content");
  plist_.set("conserved quantity B", "energy");
  plist_.set("primary variable A", "pressure");
  plist_.set("primary variable B", "temperature");

  plist_.set("mesh key", "domain");
  MPCCoupledCells::setup(S);

  // stash the PKs
  flow_pk_ = Teuchos::rcp_dynamic_cast<Amanzi::Flow::Richards>(sub_pks_[0]);
  energy_pk_ = Teuchos::rcp_dynamic_cast<Amanzi::Energy::TwoPhase>(sub_pks_[1]);

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
  }

  // select the method used for nonlinear prediction
  std::string predictor_string = plist_.get<std::string>("predictor type", "none");
  if (predictor_string == "none") {
    predictor_type_ = PREDICTOR_NONE;
  } else if (predictor_string == "heuristic") {
    predictor_type_ = PREDICTOR_HEURISTIC;
    modify_thaw_to_prev_ = plist_.get<bool>("modify thawing cells to previous temp",
            true);
  } else if (predictor_string == "ewc") {
    predictor_type_ = PREDICTOR_EWC;
  } else if (predictor_string == "smart ewc") {
    predictor_type_ = PREDICTOR_SMART_EWC;
  } else {
    Errors::Message message(std::string("Invalid predictor type ")+predictor_string);
    Exceptions::amanzi_throw(message);
  }

  // Need old values for the EWC predictor
  if (predictor_type_ == PREDICTOR_EWC || predictor_type_ == PREDICTOR_SMART_EWC) {
    S->RequireField("prev_water_content",name_)->SetMesh(S->GetMesh())->SetGhosted(false)
      ->SetComponent("cell",AmanziMesh::CELL,1);
    S->RequireField("prev_energy",name_)->SetMesh(S->GetMesh())->SetGhosted(false)
      ->SetComponent("cell",AmanziMesh::CELL,1);
  }

  // Smart EWC uses a heuristic to guess when we need the EWC instead of using
  // it blindly.
  if (predictor_type_ == PREDICTOR_SMART_EWC || precon_type_ == PRECON_SMART_EWC) {
    if (plist_.isParameter("cusp distance in T")) {
      cusp_size_T_freezing_ = plist_.get<double>("cusp distance in T");
      cusp_size_T_thawing_ = cusp_size_T_freezing_;
    } else {
      cusp_size_T_freezing_ = plist_.get<double>("cusp distance in T, freezing", 0.005);
      cusp_size_T_thawing_ = plist_.get<double>("cusp distance in T, thawing", 0.005);
    }
  }
}


// -------------------------------------------------------------
// -- Initialize owned (dependent) variables.
// -------------------------------------------------------------
void MPCFrozenCoupledFlowEnergy::initialize(const Teuchos::Ptr<State>& S) {
  // This is an ugly interpolation hack to start from frozen without actually
  // running the freezing simulation.
  if (plist_.get<bool>("initialize from frozen column", false)) {
    InitialConditionFromFrozenColumn_(S);
  }

  MPCCoupledCells::initialize(S);

  // initialize the EWC model if needed
  if (predictor_type_ == PREDICTOR_EWC ||
      predictor_type_ == PREDICTOR_SMART_EWC ||
      precon_type_ == PRECON_EWC ||
      precon_type_ == PRECON_SMART_EWC) {
    InitializeModels_(S);
  }

  if (precon_type_ == PRECON_EWC || precon_type_ == PRECON_SMART_EWC) {
    int ncells = S->GetFieldData("pressure")->size("cell",false);
    jac_.resize(ncells, WhetStone::Tensor(2,2));
  }

  if (predictor_type_ == PREDICTOR_EWC || predictor_type_ == PREDICTOR_SMART_EWC) {
    *S->GetFieldData("prev_water_content",name_) = *S->GetFieldData("water_content");
    S->GetField("prev_water_content",name_)->set_initialized();
    *S->GetFieldData("prev_energy",name_) = *S->GetFieldData("energy");
    S->GetField("prev_energy",name_)->set_initialized();

    // work is done in S_inter -- get the T and p evaluators
    Teuchos::RCP<FieldEvaluator> Teval_fe = S->GetFieldEvaluator("temperature");
    Teuchos::RCP<FieldEvaluator> peval_fe = S->GetFieldEvaluator("pressure");
    Teuchos::RCP<PrimaryVariableFieldEvaluator> Teval =
      Teuchos::rcp_static_cast<PrimaryVariableFieldEvaluator>(Teval_fe);
    Teuchos::RCP<PrimaryVariableFieldEvaluator> peval =
      Teuchos::rcp_static_cast<PrimaryVariableFieldEvaluator>(peval_fe);
    Teval->SetFieldAsChanged(S);
    peval->SetFieldAsChanged(S);
  }
}


void MPCFrozenCoupledFlowEnergy::set_states(const Teuchos::RCP<const State>& S,
          const Teuchos::RCP<State>& S_inter,
        const Teuchos::RCP<State>& S_next) {
  MPCCoupledCells::set_states(S,S_inter,S_next);

  // update derivatives -- this ensures the derivatives exist in all states, not just S_next
  S_next_->GetFieldEvaluator("water_content")
    ->HasFieldDerivativeChanged(S_next_.ptr(),name_,"temperature");
  S_next_->GetFieldEvaluator("water_content")
    ->HasFieldDerivativeChanged(S_next_.ptr(),name_,"pressure");
  S_next_->GetFieldEvaluator("energy")
    ->HasFieldDerivativeChanged(S_next_.ptr(),name_,"temperature");
  S_next_->GetFieldEvaluator("energy")
    ->HasFieldDerivativeChanged(S_next_.ptr(),name_,"pressure");

  if (predictor_type_ == PREDICTOR_EWC || predictor_type_ == PREDICTOR_SMART_EWC) {
    S_work_ = Teuchos::rcp(new State(*S_next_));
    *S_work_ = *S_next_;
  }
}


// -----------------------------------------------------------------------------
// This commit state simply saves some meta-data for the EWC predictor.
// -----------------------------------------------------------------------------
void MPCFrozenCoupledFlowEnergy::commit_state(double dt, const Teuchos::RCP<State>& S) {
  MPCCoupledCells::commit_state(dt,S);

  if ((predictor_type_ == PREDICTOR_EWC || predictor_type_ == PREDICTOR_SMART_EWC)
      && S_inter_ != Teuchos::null) {
    // stash water content and energy in S_work.
    *S_work_->GetFieldData("prev_water_content",name_) =
      *S_inter_->GetFieldData("water_content");
    *S_work_->GetFieldData("prev_energy",name_) =
      *S_inter_->GetFieldData("energy");
    S_work_->set_time(S_inter_->time());
  }

  if (predictor_type_ == PREDICTOR_EWC || predictor_type_ == PREDICTOR_SMART_EWC ||
      precon_type_ == PRECON_EWC || precon_type_ == PRECON_SMART_EWC) {
    // set model's value scalars and check it is setup
    model_->set_p_atm(*S->GetScalarData("atmospheric_pressure"));
    model_->set_rock_density(*S->GetScalarData("density_rock"));
    ASSERT(model_->IsSetUp());
  }
}


// -----------------------------------------------------------------------------
// Wrapper to call the requested modify_predictor().
// -----------------------------------------------------------------------------
bool MPCFrozenCoupledFlowEnergy::modify_predictor(double h, Teuchos::RCP<TreeVector> u) {
  bool modified = false;
  if (predictor_type_ == PREDICTOR_HEURISTIC) {
    modified = modify_predictor_heuristic_(h,u);
  } else if (predictor_type_ == PREDICTOR_EWC) {
    if (S_->time() == 0.) { // this needs to be fixed!
      return false;
    } else {
      modified = modify_predictor_ewc_(h,u);
    }
  } else if (predictor_type_ == PREDICTOR_SMART_EWC) {
    if (S_->time() == 0.) { // this needs to be fixed!
      return false;
    } else {
      modified = modify_predictor_smart_ewc_(h,u);
    }
  }

  if (modified) {
    // calculate consistent faces
    if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_EXTREME, true)) {
      *out_ << "Calculating consistent faces to go along with modified cells." << std::endl;
    }

    Teuchos::RCP<TreeVector> temp_guess = u->SubVector("energy");
    Teuchos::RCP<TreeVector> pres_guess = u->SubVector("flow");
    flow_pk_->CalculateConsistentFaces(pres_guess->data().ptr());
    energy_pk_->CalculateConsistentFaces(temp_guess->data().ptr());
  }

  return modified;
};


// -----------------------------------------------------------------------------
// Modify the initial guess for the next step via EWC.
// -----------------------------------------------------------------------------
bool MPCFrozenCoupledFlowEnergy::modify_predictor_ewc_(double h, Teuchos::RCP<TreeVector> u) {

  Teuchos::OSTab tab = getOSTab();

  Teuchos::RCP<CompositeVector> temp_guess = u->SubVector("energy")->data();
  Teuchos::RCP<CompositeVector> pres_guess = u->SubVector("flow")->data();
  Epetra_MultiVector& temp_guess_c = *temp_guess->ViewComponent("cell",false);
  Epetra_MultiVector& pres_guess_c = *pres_guess->ViewComponent("cell",false);

  const Epetra_MultiVector& T1 = *S_inter_->GetFieldData("temperature")
      ->ViewComponent("cell",false);
  const Epetra_MultiVector& p1 = *S_inter_->GetFieldData("pressure")
      ->ViewComponent("cell",false);

  const Epetra_MultiVector& cv = *S_inter_->GetFieldData("cell_volume")
      ->ViewComponent("cell",false);

  // project energy and water content
  double dt_next = S_next_->time() - S_inter_->time();
  double dt_prev = S_inter_->time() - S_work_->time();

  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_EXTREME, true)) {
    *out_ << "Projecting Energy and WC, dt_prev = " << dt_prev
          << ", dt_next = " << dt_next << std::endl;
  }

  // -- get wc and energy data
  const Epetra_MultiVector& wc0 = *S_work_->GetFieldData("prev_water_content")
      ->ViewComponent("cell",false);
  const Epetra_MultiVector& wc1 = *S_inter_->GetFieldData("water_content")
      ->ViewComponent("cell",false);
  Epetra_MultiVector& wc2 = *S_work_->GetFieldData("water_content", "water_content")
      ->ViewComponent("cell",false);

  const Epetra_MultiVector& e0 = *S_work_->GetFieldData("prev_energy")
      ->ViewComponent("cell",false);
  const Epetra_MultiVector& e1 = *S_inter_->GetFieldData("energy")
      ->ViewComponent("cell",false);
  Epetra_MultiVector& e2 = *S_work_->GetFieldData("energy", "energy")
      ->ViewComponent("cell",false);

  // -- project
  wc2 = wc0;
  e2 = e0;
  double dt_ratio = (dt_next + dt_prev) / dt_prev;
  wc2.Update(dt_ratio, wc1, 1. - dt_ratio);
  e2.Update(dt_ratio, e1, 1. - dt_ratio);

  // Check and update
  const Epetra_MultiVector& poro = *S_next_->GetFieldData("porosity")
      ->ViewComponent("cell",false);

  double wc_scale = 10.;
  double e_scale = 10000.;

  int ncells = wc0.MyLength();
  for (int c=0; c!=ncells; ++c) {
    // determine how to project based upon min change in e,wc
    AmanziGeometry::Point result(2);
    int ierr = 0;

    double p = p1[0][c];
    double T = T1[0][c];

    // uses intensive forms, so must divide by cell volume.
#if DEBUG_FLAG
    std::cout << "Inverting: c = " << c << std::endl;
    std::cout << "   p,T  = " << p << ", " << T << std::endl;
    std::cout << "   wc,e = " << wc1[0][c] << ", " << e1[0][c] << std::endl;
    std::cout << "   goal = " << wc2[0][c] << ", " << e2[0][c] << std::endl;
#endif
    ierr = model_->InverseEvaluate(e2[0][c]/cv[0][c], wc2[0][c]/cv[0][c], poro[0][c], T, p);

    if (!ierr) { // valid solution, no zero determinates, etc
      temp_guess_c[0][c] = T;
      pres_guess_c[0][c] = p;
    }
  }
  return true;
}


// -----------------------------------------------------------------------------
// Modify the initial guess for the next step via selective use of EWC.
// -----------------------------------------------------------------------------
bool MPCFrozenCoupledFlowEnergy::modify_predictor_smart_ewc_(double h, Teuchos::RCP<TreeVector> u) {

  Teuchos::OSTab tab = getOSTab();

  Teuchos::RCP<CompositeVector> temp_guess = u->SubVector("energy")->data();
  Teuchos::RCP<CompositeVector> pres_guess = u->SubVector("flow")->data();
  Epetra_MultiVector& temp_guess_c = *temp_guess->ViewComponent("cell",false);
  Epetra_MultiVector& pres_guess_c = *pres_guess->ViewComponent("cell",false);

  const Epetra_MultiVector& T1 = *S_inter_->GetFieldData("temperature")
      ->ViewComponent("cell",false);
  const Epetra_MultiVector& p1 = *S_inter_->GetFieldData("pressure")
      ->ViewComponent("cell",false);

  const Epetra_MultiVector& cv = *S_inter_->GetFieldData("cell_volume")
      ->ViewComponent("cell",false);

  // DEBUG CRUFT
#if DEBUG_FLAG
  for (std::vector<int>::const_iterator lcv = cells_to_track_.begin();
       lcv!=cells_to_track_.end(); ++lcv) {
    std::cout << "Check on Prev guesses: c = " << *lcv << std::endl;
    std::cout << "  true p,T = " << p1[0][*lcv] << ", " << T1[0][*lcv] << std::endl;
  }
  cells_to_track_.clear();
#endif

  // project energy and water content
  double dt_next = S_next_->time() - S_inter_->time();
  double dt_prev = S_inter_->time() - S_work_->time();

  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_EXTREME, true)) {
    *out_ << "Projecting Energy and WC, dt_prev = " << dt_prev
          << ", dt_next = " << dt_next << std::endl;
  }

  // -- get wc and energy data
  const Epetra_MultiVector& wc0 = *S_work_->GetFieldData("prev_water_content")
      ->ViewComponent("cell",false);
  const Epetra_MultiVector& wc1 = *S_inter_->GetFieldData("water_content")
      ->ViewComponent("cell",false);
  Epetra_MultiVector& wc2 = *S_work_->GetFieldData("water_content", "water_content")
      ->ViewComponent("cell",false);

  const Epetra_MultiVector& e0 = *S_work_->GetFieldData("prev_energy")
      ->ViewComponent("cell",false);
  const Epetra_MultiVector& e1 = *S_inter_->GetFieldData("energy")
      ->ViewComponent("cell",false);
  Epetra_MultiVector& e2 = *S_work_->GetFieldData("energy", "energy")
      ->ViewComponent("cell",false);

  // -- project
  wc2 = wc0;
  e2 = e0;
  double dt_ratio = (dt_next + dt_prev) / dt_prev;
  wc2.Update(dt_ratio, wc1, 1. - dt_ratio);
  e2.Update(dt_ratio, e1, 1. - dt_ratio);

  // Check and update
  const Epetra_MultiVector& poro = *S_next_->GetFieldData("porosity")
      ->ViewComponent("cell",false);

  double wc_scale = 10.;
  double e_scale = 10000.;

  int ncells = wc0.MyLength();
  for (int c=0; c!=ncells; ++c) {
    bool ewc_predictor = false;
    AmanziGeometry::Point result(2);
    int ierr = 0;

    // Determine whether to use T/p or E/WC as the projection
    double T_prev = T1[0][c];
    double T_guess = temp_guess_c[0][c];
    // invert the T-projection to get what it was projected from
    double T_prev2 = (T_guess - dt_ratio*T_prev) / (1. - dt_ratio);

    // Curve looks like __|--    with discontinuity at 0.
    if (T_guess < T_prev) { // getting colder, so the upper cusp is the problem
      if (T_guess > 273.15) {
        ewc_predictor = false;  // both above freezing, in the upper branch
      } else if (T_prev2 < 273.15 - cusp_size_T_freezing_) {
        ewc_predictor = false;  // both well below freezing, in lower branch
      } else if (T_prev2 - T_guess < cusp_size_T_freezing_) {
        ewc_predictor = false;  // likely past the upper cusp
      } else {
        ewc_predictor = true;
      }
    } else if (T_prev < T_guess) { // getting warmer, so the lower cusp is the problem
      if (T_guess < 273.15 - cusp_size_T_thawing_) {
        ewc_predictor = false;  // both below freezing, in the lower branch
      } else if (T_prev2 > 273.15) {
        ewc_predictor = false;  // both well above freezing, in upper branch
      } else if (T_guess - T_prev2 < cusp_size_T_thawing_) {
        ewc_predictor = false;  // likely past the lowercusp
      } else {
        ewc_predictor = true;
      }
    } else {
      ewc_predictor = false;
    }

    if (ewc_predictor) {
      double p = p1[0][c];
      double T = T1[0][c];

#if DEBUG_FLAG
      std::cout << "Inverting: c = " << c << std::endl;
      std::cout << "   based upon h_old = " << dt_prev << ", h_next = " << dt_next << std::endl;
      std::cout << "   Prev wc,e: " << wc1[0][c] << ", " << e1[0][c] << std::endl;
      std::cout << "   Prev p,T: " << p << ", " << T << std::endl;
      std::cout << "   Desired wc,e: " << wc2[0][c] << ", " << e2[0][c] << std::endl;
#endif

      // uses intensive forms, so must divide by cell volume.
      ierr = model_->InverseEvaluate(e2[0][c]/cv[0][c], wc2[0][c]/cv[0][c], poro[0][c], T, p);

#if DEBUG_FLAG
      std::cout << "   p_ewc,T_ewc  = " << p << ", " << T << std::endl;
      cells_to_track_.push_back(c);
#endif

      if (!ierr) { // valid solution, no zero determinates, etc
        // Return the smaller change
        if (std::abs(T - T_prev) < std::abs(T_guess - T_prev)) {
#if DEBUG_FLAG
          std::cout << "   EWC Accepted" << std::endl;
#endif
          temp_guess_c[0][c] = T;
          pres_guess_c[0][c] = p;
        } else {
#if DEBUG_FLAG
          std::cout << "   EWC Rejected" << std::endl;
#endif
        }
      }
    }
  }
  return true;
}


// -----------------------------------------------------------------------------
// Modify the initial guess for the next step via a heuristic around freezing.
// -----------------------------------------------------------------------------
bool MPCFrozenCoupledFlowEnergy::modify_predictor_heuristic_(double h,
        Teuchos::RCP<TreeVector> u) {
  Teuchos::RCP<CompositeVector> temp_guess = u->SubVector("energy")->data();
  Epetra_MultiVector& temp_guess_c = *temp_guess->ViewComponent("cell",false);
  Epetra_MultiVector& guess_faces = *temp_guess->ViewComponent("face",false);

  Teuchos::RCP<const CompositeVector> temp = S_inter_->GetFieldData("temperature");
  Teuchos::RCP<const CompositeVector> pres = S_inter_->GetFieldData("pressure");
  const Epetra_MultiVector& temp_c = *temp->ViewComponent("cell",false);
  const Epetra_MultiVector& temp_faces = *temp->ViewComponent("face",false);

  bool update_faces(false);


  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_EXTREME, true)) {
    double minT,maxT,minTL,maxTL;
    temp_guess_c.MinValue(&minT);
    temp_guess_c.MaxValue(&maxT);
    guess_faces.MinValue(&minTL);
    guess_faces.MaxValue(&maxTL);

    *out_ << "--- Modifying Guess: ---" << std::endl;
    *out_ << "  T_min = " << minT << "     T_max = " << maxT << std::endl;
    *out_ << "  Lambda_min = " << minTL << "     Lambda_max = " << maxTL << std::endl;
  }


  int ncells = temp->size("cell",false);
  std::vector<bool> changed(ncells, false);

  for (int c=0; c!=ncells; ++c) {
    if (temp_c[0][c] >= 273.15 && temp_guess_c[0][c] < 273.15) {
      // freezing
#if DEBUG_FLAG
      std::cout << "Freezing cell " << c << std::endl;
      std::cout << "   T_prev = " << temp_c[0][c]
               << ",  T_guess = " << temp_guess_c[0][c]
               << ",  T_corrected = " << 273.15 - 1.e-3 << std::endl;
#endif
      temp_guess_c[0][c] = 273.15 - 1.e-3;
      changed[c] = true;
      update_faces = true;

    } else if (273.15 > temp_c[0][c] &&
               temp_c[0][c] >= 273.1 &&
               (temp_c[0][c] - temp_guess_c[0][c]) > 1.e-2) {
      // catch the 2nd step in freezing -- after the 2nd step the
      // extrapolation should be ok?
#if DEBUG_FLAG
        std::cout << "2nd Freezing step cell " << c << std::endl;
        std::cout << "   T_prev = " << temp_c[0][c]
              << ",  T_guess = " << temp_guess_c[0][c]
              << ",  T_corrected = " << temp_c[0][c] << std::endl;
#endif
      temp_guess_c[0][c] = temp_c[0][c];
      changed[c] = true;
      update_faces = true;

    } else if (temp_c[0][c] <= 273.15 && temp_guess_c[0][c] > 273.15) {
      if (modify_thaw_to_prev_) {
        // thawing
#if DEBUG_FLAG
        std::cout << "Thawing cell " << c << std::endl;
        std::cout << "   T_prev = " << temp_c[0][c]
                  << ",  T_guess = " << temp_guess_c[0][c]
                  << ",  T_corrected = " << temp_c[0][c] << std::endl;
#endif
        temp_guess_c[0][c] = temp_c[0][c];
        changed[c] = true;
        update_faces = true;

      } else {
        // thawing, 2nd step
#if DEBUG_FLAG
        std::cout << "Thawing cell " << c << std::endl;
        std::cout << "   T_prev = " << temp_c[0][c]
                  << ",  T_guess = " << temp_guess_c[0][c]
                  << ",  T_corrected = " << 273.15 - 1.e-3 << std::endl;
#endif
        temp_guess_c[0][c] = std::min((temp_guess_c[0][c] + temp_c[0][c]) / 2., 273.14);
        changed[c] = true;
        update_faces = true;
      }
    }
  }
  return true;
}


// -----------------------------------------------------------------------------
// Wrapper to call the requested update_precon().
// -----------------------------------------------------------------------------
void MPCFrozenCoupledFlowEnergy::update_precon(double t,
        Teuchos::RCP<const TreeVector> up, double h) {
  if (precon_type_ == PRECON_NONE) {
    return;
  } else if (precon_type_ == PRECON_BLOCK_DIAGONAL) {
    return StrongMPC::update_precon(t,up,h);
  } else if (precon_type_ == PRECON_PICARD) {
    return MPCCoupledCells::update_precon(t,up,h);
  } else if (precon_type_ == PRECON_EWC) {
    return update_precon_ewc_(t,up,h);
  } else if (precon_type_ == PRECON_SMART_EWC) {
    return update_precon_ewc_(t,up,h);
  } else {
    ASSERT(0);
  }
}


// -----------------------------------------------------------------------------
// Update extra parts for the EWC precon
// -----------------------------------------------------------------------------
void MPCFrozenCoupledFlowEnergy::update_precon_ewc_(double t,
        Teuchos::RCP<const TreeVector> up, double h) {
  MPCCoupledCells::update_precon(t,up,h);

  const Epetra_MultiVector& dedT = *S_next_->GetFieldData("denergy_dtemperature")
      ->ViewComponent("cell",false);
  const Epetra_MultiVector& dedp = *S_next_->GetFieldData("denergy_dpressure")
      ->ViewComponent("cell",false);
  const Epetra_MultiVector& dwcdT = *S_next_->GetFieldData("dwater_content_dtemperature")
      ->ViewComponent("cell",false);
  const Epetra_MultiVector& dwcdp = *S_next_->GetFieldData("dwater_content_dpressure")
      ->ViewComponent("cell",false);

  int ncells = up->SubVector("energy")->data()->size("cell", false);
  for (int c=0; c!=ncells; ++c) {
    jac_[c](0,0) = dedT[0][c];
    jac_[c](0,1) = dedp[0][c];
    jac_[c](1,0) = dwcdT[0][c];
    jac_[c](1,1) = dwcdp[0][c];
  }
}


// -----------------------------------------------------------------------------
// Wrapper to call the requested preconditioner.
// -----------------------------------------------------------------------------
void MPCFrozenCoupledFlowEnergy::precon(Teuchos::RCP<const TreeVector> u,
        Teuchos::RCP<TreeVector> Pu) {
  if (precon_type_ == PRECON_NONE) {
    *Pu = *u;
    return;
  } else if (precon_type_ == PRECON_BLOCK_DIAGONAL) {
    return StrongMPC::precon(u,Pu);
  } else if (precon_type_ == PRECON_PICARD) {
    return MPCCoupledCells::precon(u,Pu);
  } else if (precon_type_ == PRECON_EWC) {
    return precon_ewc_(u, Pu);
  } else if (precon_type_ == PRECON_SMART_EWC) {
    return precon_smart_ewc_(u, Pu);
  } else {
    ASSERT(0);
  }
}

// -----------------------------------------------------------------------------
// The E-WC preconditioner is composed by:
//   -- using precon() to calculate a correction, (dT,dp)
//   -- calculating a correction in the new vars, (de,dwc) = J(dT,dp)
//       where J is the Jacobian of the E-WC transformation.
//   -- calculating the corrected (e',wc') = (e_prev - de, wc_prev - dwc)
//   -- inverting for corrected (T',p') = E-WC^inv(e',wc')
//   -- Pu = (T' - T_prev, p' - p_prev)
// -----------------------------------------------------------------------------
void MPCFrozenCoupledFlowEnergy::precon_ewc_(Teuchos::RCP<const TreeVector> u,
        Teuchos::RCP<TreeVector> Pu) {
  // using precon() to calculate a correction, (dT,dp)
  MPCCoupledCells::precon(u,Pu);

  Teuchos::RCP<TreeVector> dtemp = Pu->SubVector("energy");
  Teuchos::RCP<TreeVector> dpres = Pu->SubVector("flow");
  Epetra_MultiVector& dtemp_c = *dtemp->data()->ViewComponent("cell",false);
  Epetra_MultiVector& dpres_c = *dpres->data()->ViewComponent("cell",false);

  const Epetra_MultiVector& temp_prev = *S_next_->GetFieldData("temperature")
      ->ViewComponent("cell",false);
  const Epetra_MultiVector& pres_prev = *S_next_->GetFieldData("pressure")
      ->ViewComponent("cell",false);
  const Epetra_MultiVector& poro = *S_next_->GetFieldData("porosity")
      ->ViewComponent("cell",false);
  const Epetra_MultiVector& e_prev = *S_next_->GetFieldData("energy")
      ->ViewComponent("cell",false);
  const Epetra_MultiVector& wc_prev = *S_next_->GetFieldData("water_content")
      ->ViewComponent("cell",false);
  const Epetra_MultiVector& cell_volume = *S_next_->GetFieldData("cell_volume")
      ->ViewComponent("cell",false);

  Teuchos::RCP<TreeVector> temp_new = Teuchos::rcp(new TreeVector(*dtemp));
  Teuchos::RCP<TreeVector> pres_new = Teuchos::rcp(new TreeVector(*dpres));
  Epetra_MultiVector& temp_new_c = *temp_new->data()->ViewComponent("cell",false);
  Epetra_MultiVector& pres_new_c = *pres_new->data()->ViewComponent("cell",false);

  AmanziGeometry::Point dTp(2), dewc(2);

  int ncells = dtemp_c.MyLength();
  for (int c=0; c!=ncells; ++c) {

    dTp[0] = dtemp_c[0][c]; dTp[1] = dpres_c[0][c];
    dewc = jac_[c] * dTp;

    // calculating the corrected (e',wc') = (e_prev - de, wc_prev - dwc)
    //   note _int suffix imples intensive, per unit volume.
    double e_new_int = e_prev[0][c]/cell_volume[0][c] - dewc[0];
    double wc_new_int = wc_prev[0][c]/cell_volume[0][c] - dewc[1];

    // inverting for corrected (T',p') = E-WC^inv(e',wc')
    double T_new = temp_prev[0][c];
    double p_new = pres_prev[0][c];

    // using intensive values for inversion
    int ierr = model_->InverseEvaluate(e_new_int/cell_volume[0][c],
            wc_new_int/cell_volume[0][c], poro[0][c], T_new, p_new);

    if (!ierr) { // valid solution, no zero determinates, etc
      temp_new_c[0][c] = T_new;
      pres_new_c[0][c] = p_new;
    } else {
      // fall back on the standard correction
      temp_new_c[0][c] = temp_prev[0][c] - dtemp_c[0][c];
      pres_new_c[0][c] = pres_prev[0][c] - dpres_c[0][c];
    }
  }

  // Calculate consistent faces
  double h = S_next_->time() - S_inter_->time();
  energy_pk_->CalculateConsistentFaces(temp_new->data().ptr());
  flow_pk_->CalculateConsistentFaces(pres_new->data().ptr());

  // calculate the correction
  *dtemp->data() = *S_next_->GetFieldData("temperature");
  dtemp->data()->Update(-1., *temp_new->data(), 1.);
  *dpres->data() = *S_next_->GetFieldData("pressure");
  dpres->data()->Update(-1., *pres_new->data(), 1.);
}


// -----------------------------------------------------------------------------
// The E-WC preconditioner is composed by:
//   -- using precon() to calculate a correction, (dT,dp)
//   -- calculating a correction in the new vars, (de,dwc) = J(dT,dp)
//       where J is the Jacobian of the E-WC transformation.
//   -- calculating the corrected (e',wc') = (e_prev - de, wc_prev - dwc)
//   -- inverting for corrected (T',p') = E-WC^inv(e',wc')
//   -- Pu = (T' - T_prev, p' - p_prev)
//
// The smart EW-C preconditioner is done by only doing the E-WC preconditioner
// when it is useful, i.e. temperature is decreasing into the upper cusp of
// the freezing point or increasing into the lower cusp of the thawing point.
// -----------------------------------------------------------------------------
void MPCFrozenCoupledFlowEnergy::precon_smart_ewc_(Teuchos::RCP<const TreeVector> u,
        Teuchos::RCP<TreeVector> Pu) {
  // using precon() to calculate a correction, (dT,dp)
  MPCCoupledCells::precon(u,Pu);

  // debugging
  const Epetra_MultiVector& res_temp_c = *u->SubVector("energy")->data()
      ->ViewComponent("cell",false);
  const Epetra_MultiVector& res_pres_c = *u->SubVector("flow")->data()
      ->ViewComponent("cell",false);

  Teuchos::RCP<TreeVector> dtemp = Pu->SubVector("energy");
  Teuchos::RCP<TreeVector> dpres = Pu->SubVector("flow");
  Epetra_MultiVector& dtemp_c = *dtemp->data()->ViewComponent("cell",false);
  Epetra_MultiVector& dpres_c = *dpres->data()->ViewComponent("cell",false);

  const Epetra_MultiVector& temp_prev = *S_next_->GetFieldData("temperature")
      ->ViewComponent("cell",false);
  const Epetra_MultiVector& pres_prev = *S_next_->GetFieldData("pressure")
      ->ViewComponent("cell",false);
  const Epetra_MultiVector& poro = *S_next_->GetFieldData("porosity")
      ->ViewComponent("cell",false);
  const Epetra_MultiVector& e_prev = *S_next_->GetFieldData("energy")
      ->ViewComponent("cell",false);
  const Epetra_MultiVector& wc_prev = *S_next_->GetFieldData("water_content")
      ->ViewComponent("cell",false);
  const Epetra_MultiVector& cell_volume = *S_next_->GetFieldData("cell_volume")
      ->ViewComponent("cell",false);

  Teuchos::RCP<TreeVector> temp_new = Teuchos::rcp(new TreeVector(*dtemp));
  Teuchos::RCP<TreeVector> pres_new = Teuchos::rcp(new TreeVector(*dpres));
  Epetra_MultiVector& temp_new_c = *temp_new->data()->ViewComponent("cell",false);
  Epetra_MultiVector& pres_new_c = *pres_new->data()->ViewComponent("cell",false);

  AmanziGeometry::Point dTp(2), dewc(2);

  bool ewc_any = false;
  int ncells = dtemp_c.MyLength();
  for (int c=0; c!=ncells; ++c) {
    // determine which precon to use
    bool ewc_precon = true;

    // Determine whether to use P or P * J (where P is the Picard precon and J
    // is the Jacobian of the E-WC transform).
    double T_prev = temp_prev[0][c];
    double T_picard = T_prev - dtemp_c[0][c];

    //                     __
    // Curve looks like __|    with discontinuity at 0.
    if (T_picard < T_prev) { // getting colder, so the upper cusp is the problem
      if (T_picard > 273.15) {
        ewc_precon = false;  // both above freezing, in the upper branch
      } else if (T_prev < 273.15 - cusp_size_T_freezing_) {
        ewc_precon = false;  // both well below freezing, in lower branch
      } else if (T_prev - T_picard < cusp_size_T_freezing_) {
        ewc_precon = false;  // likely past the upper cusp
      } else {
        ewc_precon = true;
      }
    } else if (T_prev < T_picard) { // getting warmer, so the lower cusp is the problem
      if (T_picard < 273.15 - cusp_size_T_thawing_) {
        ewc_precon = false;  // both below freezing, in the lower branch
      } else if (T_prev > 273.15) {
        ewc_precon = false;  // both well above freezing, in upper branch
      } else if (T_picard - T_prev < cusp_size_T_thawing_) {
        ewc_precon = false;  // likely past the lowercusp
      } else {
        ewc_precon = true;
      }
    } else {
      ewc_precon = false;
    }

    if (ewc_precon) {
      ewc_any = true;

#if DEBUG_FLAG
      std::cout << "-----------------------------------------" << std::endl;
      std::cout << "Using EWC Precon on cell " << c << std::endl;
      std::cout << " with residual (T,p) = " << res_temp_c[0][c] << ", "
                << res_pres_c[0][c] << std::endl;
      std::cout << "  T,p_old = " << temp_prev[0][c] << ", " << pres_prev[0][c] << std::endl;
      std::cout << "  dT,dp = " << dtemp_c[0][c] << ", " << dpres_c[0][c] << std::endl;
      std::cout << "  T,p_new = " << temp_prev[0][c] - dtemp_c[0][c]
                << ", " << pres_prev[0][c] - dpres_c[0][c] << std::endl;
#endif

      dTp[0] = dtemp_c[0][c]; dTp[1] = dpres_c[0][c];
      dewc = jac_[c] * dTp;

      // calculating the corrected (e',wc') = (e_prev - de, wc_prev - dwc)
      //   note _int suffix imples intensive, per unit volume.
      double e_new_int = e_prev[0][c] - dewc[0];
      double wc_new_int = wc_prev[0][c] - dewc[1];

#if DEBUG_FLAG
      std::cout << "  --- " << std::endl;
      std::cout << " Jac = [ " << jac_[c](0,0) << ", " << jac_[c](0,1) << "]" << std::endl;
      std::cout << "       [ " << jac_[c](1,0) << ", " << jac_[c](1,1) << "]" << std::endl;
      std::cout << "  e,wc_old = " << e_prev[0][c] << ", " << wc_prev[0][c] << std::endl;
      std::cout << "  de,dwc = " << dewc[0] << ", " << dewc[1] << std::endl;
      std::cout << "  e,wc_new = " << e_new_int << ", " << wc_new_int << std::endl;
#endif

      // inverting for corrected (T',p') = E-WC^inv(e',wc')
      double T_new = temp_prev[0][c];
      double p_new = pres_prev[0][c];

      // using intensive values for inversion
      int ierr = model_->InverseEvaluate(e_new_int/cell_volume[0][c],
              wc_new_int/cell_volume[0][c], poro[0][c], T_new, p_new);

#if DEBUG_FLAG
      std::cout << "  T,p_new = " << T_new << ", " << p_new << std::endl;
#endif

      //      ASSERT(!ierr); // remove me for robustness --etc
      if (!ierr) { // valid solution, no zero determinates, etc
        temp_new_c[0][c] = T_new;
        pres_new_c[0][c] = p_new;
      } else {
        // fall back on the standard correction
        temp_new_c[0][c] = temp_prev[0][c] - dtemp_c[0][c];
        pres_new_c[0][c] = pres_prev[0][c] - dpres_c[0][c];
      }

    } else {
      // the standard correction
      temp_new_c[0][c] = temp_prev[0][c] - dtemp_c[0][c];
      pres_new_c[0][c] = pres_prev[0][c] - dpres_c[0][c];
    }
  }

  if (ewc_any) {
    // update with consistent faces
    double h = S_next_->time() - S_inter_->time();
    energy_pk_->CalculateConsistentFaces(temp_new->data().ptr());
    flow_pk_->CalculateConsistentFaces(pres_new->data().ptr());

    // calculate the correction
    *dtemp->data() = *S_next_->GetFieldData("temperature");
    dtemp->data()->Update(-1., *temp_new->data(), 1.);
    *dpres->data() = *S_next_->GetFieldData("pressure");
    dpres->data()->Update(-1., *pres_new->data(), 1.);
  }
}


// -----------------------------------------------------------------------------
// This permafrost model is needed by EWC to invert energy and water content
// for pressure and temperature.
// -----------------------------------------------------------------------------
void MPCFrozenCoupledFlowEnergy::InitializeModels_(const Teuchos::Ptr<State>& S) {
  // get the WRM models and their regions
  Teuchos::RCP<FieldEvaluator> me = S->GetFieldEvaluator("saturation_gas");
  Teuchos::RCP<Flow::FlowRelations::WRMPermafrostEvaluator> wrm_me =
      Teuchos::rcp_dynamic_cast<Flow::FlowRelations::WRMPermafrostEvaluator>(me);
  ASSERT(wrm_me != Teuchos::null);
  Teuchos::RCP<Flow::FlowRelations::WRMPermafrostModelPartition> wrms =
      wrm_me->get_WRMPermafrostModels();

  // this needs fixed eventually, but for now assuming one WRM, and therefore
  // one model --etc
  ASSERT(wrms->second.size() == 1);

  // -- WRMs
  model_ = Teuchos::rcp(new PermafrostModel());
  model_->set_WRM(wrms->second[0]);

  // -- liquid EOS
  me = S->GetFieldEvaluator("molar_density_liquid");
  Teuchos::RCP<Relations::EOSEvaluator> eos_liquid_me =
      Teuchos::rcp_dynamic_cast<Relations::EOSEvaluator>(me);
  ASSERT(eos_liquid_me != Teuchos::null);
  Teuchos::RCP<Relations::EOS> eos = eos_liquid_me->get_EOS();
  model_->set_liquid_EOS(eos);

  // -- ice EOS
  me = S->GetFieldEvaluator("molar_density_ice");
  Teuchos::RCP<Relations::EOSEvaluator> eos_ice_me =
      Teuchos::rcp_dynamic_cast<Relations::EOSEvaluator>(me);
  ASSERT(eos_ice_me != Teuchos::null);
  eos = eos_ice_me->get_EOS();
  model_->set_ice_EOS(eos);

  // -- gas EOS
  me = S->GetFieldEvaluator("molar_density_gas");
  Teuchos::RCP<Relations::EOSEvaluator> eos_gas_me =
      Teuchos::rcp_dynamic_cast<Relations::EOSEvaluator>(me);
  ASSERT(eos_gas_me != Teuchos::null);
  eos = eos_gas_me->get_EOS();
  model_->set_gas_EOS(eos);

  // -- gas vapor pressure
  me = S->GetFieldEvaluator("mol_frac_gas");
  Teuchos::RCP<Relations::MolarFractionGasEvaluator> mol_frac_me =
      Teuchos::rcp_dynamic_cast<Relations::MolarFractionGasEvaluator>(me);
  ASSERT(mol_frac_me != Teuchos::null);
  Teuchos::RCP<Relations::VaporPressureRelation> vpr = mol_frac_me->get_VaporPressureRelation();
  model_->set_vapor_pressure_relation(vpr);

  // -- capillary pressure for ice/water
  me = S->GetFieldEvaluator("capillary_pressure_liq_ice");
  Teuchos::RCP<Flow::FlowRelations::PCIceEvaluator> pc_ice_me =
      Teuchos::rcp_dynamic_cast<Flow::FlowRelations::PCIceEvaluator>(me);
  ASSERT(pc_ice_me != Teuchos::null);
  Teuchos::RCP<Flow::FlowRelations::PCIceWater> pc_ice = pc_ice_me->get_PCIceWater();
  model_->set_pc_ice_water(pc_ice);

  // -- capillary pressure for liq/gas
  me = S->GetFieldEvaluator("capillary_pressure_gas_liq");
  Teuchos::RCP<Flow::FlowRelations::PCLiquidEvaluator> pc_liq_me =
      Teuchos::rcp_dynamic_cast<Flow::FlowRelations::PCLiquidEvaluator>(me);
  ASSERT(pc_liq_me != Teuchos::null);
  Teuchos::RCP<Flow::FlowRelations::PCLiqAtm> pc_liq = pc_liq_me->get_PCLiqAtm();
  model_->set_pc_liq_gas(pc_liq);

  // -- iem for liquid
  me = S->GetFieldEvaluator("internal_energy_liquid");
  Teuchos::RCP<Energy::EnergyRelations::IEMEvaluator> iem_liquid_me =
      Teuchos::rcp_dynamic_cast<Energy::EnergyRelations::IEMEvaluator>(me);
  ASSERT(iem_liquid_me != Teuchos::null);
  Teuchos::RCP<Energy::EnergyRelations::IEM> iem = iem_liquid_me->get_IEM();
  model_->set_liquid_IEM(iem);

  // -- iem for ice
  me = S->GetFieldEvaluator("internal_energy_ice");
  Teuchos::RCP<Energy::EnergyRelations::IEMEvaluator> iem_ice_me =
      Teuchos::rcp_dynamic_cast<Energy::EnergyRelations::IEMEvaluator>(me);
  ASSERT(iem_ice_me != Teuchos::null);
  iem = iem_ice_me->get_IEM();
  model_->set_ice_IEM(iem);

  // -- iem for gas
  me = S->GetFieldEvaluator("internal_energy_gas");
  Teuchos::RCP<Energy::EnergyRelations::IEMWaterVaporEvaluator> iem_gas_me =
      Teuchos::rcp_dynamic_cast<Energy::EnergyRelations::IEMWaterVaporEvaluator>(me);
  ASSERT(iem_gas_me != Teuchos::null);
  Teuchos::RCP<Energy::EnergyRelations::IEMWaterVapor> iem_gas = iem_gas_me->get_IEM();
  model_->set_gas_IEM(iem_gas);

  // -- iem for rock
  me = S->GetFieldEvaluator("internal_energy_rock");
  Teuchos::RCP<Energy::EnergyRelations::IEMEvaluator> iem_rock_me =
      Teuchos::rcp_dynamic_cast<Energy::EnergyRelations::IEMEvaluator>(me);
  ASSERT(iem_rock_me != Teuchos::null);
  iem = iem_rock_me->get_IEM();
  model_->set_rock_IEM(iem);

}


// -----------------------------------------------------------------------------
// The freezing problem can be a pain to run, this initializes by
// interpolating a solution from freezing a single column.
// -----------------------------------------------------------------------------
void MPCFrozenCoupledFlowEnergy::InitialConditionFromFrozenColumn_(
    const Teuchos::Ptr<State>& S) {

  Teuchos::RCP<CompositeVector> pres = S->GetFieldData("pressure", flow_pk_->name());
  Teuchos::RCP<CompositeVector> temp = S->GetFieldData("temperature", energy_pk_->name());

  int npoints = 100;
  double ref_T[] = {262.7091428,262.641515665,262.58610644,262.54456626,262.516834608,262.503140657,262.496062716,262.488999243,262.48192875,262.474849863,262.467762477,262.460666527,262.453561996,262.446448892,262.439327226,262.432196991,262.425058094,262.417910659,262.410754807,262.403590438,262.396417523,262.38923607,262.382046106,262.374847613,262.367640513,262.360424583,262.353200074,262.345967205,262.338725906,262.331476326,262.324218348,262.316951874,262.309676754,262.302392788,262.295100307,262.287799762,262.280490954,262.273173738,262.26584795,262.258513003,262.25116961,262.243818573,262.236459336,262.229091793,262.221715879,262.214331552,262.206938887,262.199537835,262.192127993,262.184709455,262.177282338,262.169845193,262.162392238,262.154903957,262.147215752,262.138831145,262.127574119,262.106832867,262.070124703,262.01912569,261.95878179,261.891743368,261.819681263,261.743550148,261.664006281,261.581541241,261.496529897,261.409277942,261.320038524,261.229029146,261.136443739,261.042456713,260.947225061,260.850890294,260.753580116,260.655409681,260.556482544,260.456891469,260.356719229,260.256039443,260.154917462,260.053411276,259.951572408,259.849446765,259.747075431,259.644495381,259.541740117,259.438840219,259.335823827,259.232717041,259.129544246,259.026328383,258.92309113,258.819853031,258.716633546,258.613451025,258.510322629,258.407264182,258.304290022,258.201412884};
  double ref_p[] = {-91144793.8313,-79198250.2123,-64553378.2237,-53817268.5748,-40376600.31,-9186846.56641,-6245241.77642,-5844385.55992,-5808523.22231,-5797105.81683,-5795148.67225,-5807088.92269,-5810167.68087,-5821627.31323,-5820860.37694,-5836780.45965,-5868415.98829,-5840987.06814,-5838227.17932,-5843736.56489,-5852397.98775,-5856537.87728,-5857887.33366,-5869968.3054,-5899256.57528,-5985661.00795,-5927751.02996,-5956801.49746,-5915396.93489,-5903438.58908,-5900368.18448,-5923657.7528,-5971066.03108,-6058214.67597,-5996882.73853,-5952227.64174,-5956081.90559,-5958021.93473,-6019441.59966,-6202372.07811,-6030838.97887,-5995013.54658,-5986821.9308,-5986523.53408,-6003380.78559,-6018426.01063,-6014029.63473,-6050617.29614,-6174214.01679,-6175627.87735,-6271705.74182,-6657225.50566,-7802907.37546,-9925062.27734,-15631485.3766,-21637300.6575,-32255189.5224,-47929250.3074,-63007803.9699,-72942365.2584,-81113258.8592,-87994540.9417,-94323588.0452,-100411902.799,-106347449.951,-112258287.406,-118190275.792,-124185038.048,-130274367.97,-136457417.992,-142720313.132,-149043914.067,-155403433.988,-161767751.951,-168100565.584,-174362689.812,-180514544.613,-186518257.197,-192339144.137,-197946557.031,-203314195.973,-208420043.01,-213246066.64,-217777821.856,-222004036.177,-225916241.135,-229508484.983,-232777146.221,-235720857.447,-238340542.772,-240639567.369,-242623992.387,-244302919.818,-245688897.123,-246798326.804,-247651788.05,-248274123.31,-248694073.748,-248943176.096,-249053590.406};

  double dz = 0.1;
  int ref_water_table_k = 53;

  double water_table = plist_.get<double>("water table height");

  Epetra_MultiVector& temp_c = *temp->ViewComponent("cell",false);
  Epetra_MultiVector& pres_c = *pres->ViewComponent("cell",false);

  int ncells = temp->size("cell",false);
  for (int c=0; c!=ncells; ++c) {
    double c_z = pres->mesh()->cell_centroid(c)[2];
    double z_over_wt = c_z - water_table;
    int dk = std::floor(z_over_wt / dz);
    int k = dk + ref_water_table_k;

    if (k < 0 || k > npoints - 2) {
      Errors::Message message("Initialization of frozen column is off the charts");
      Exceptions::amanzi_throw(message);
    }

    double p0 = ref_p[k];
    double p1 = ref_p[k+1];
    double T0 = ref_T[k];
    double T1 = ref_T[k+1];

    double param = (z_over_wt - dz*dk) / dz;

    temp_c[0][c] = T0*param + T1*(1.-param);
    pres_c[0][c] = p0*param + p1*(1.-param);

  }
}


} // namespace
