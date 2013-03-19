/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Interface for EWC, a helper class that does projections and preconditioners in
energy/water-content space instead of temperature/pressure space.
------------------------------------------------------------------------- */

#include "ewc_model.hh"
#include "mpc_delegate_ewc.hh"

namespace Amanzi {

#define DEBUG_FLAG 1

// -----------------------------------------------------------------------------
// Constructor
// -----------------------------------------------------------------------------
MPCDelegateEWC::MPCDelegateEWC(Teuchos::ParameterList& plist) :
    plist_(plist) {}

// -----------------------------------------------------------------------------
// Allocate any data or models required.
// -----------------------------------------------------------------------------
void MPCDelegateEWC::setup(const Teuchos::Ptr<State>& S) {
  // Get the mesh
  Key domain = plist_.get<std::string>("domain key", "");
  if (domain.size() != 0) {
    mesh_ = S->GetMesh(domain);
  } else {
    mesh_ = S->GetMesh("domain");
  }

  // Process the parameter list for data Keys
  if (domain.size() > 0) domain.append(1, '_');
  pres_key_ = plist_.get<std::string>("pressure key", domain+std::string("pressure"));
  temp_key_ = plist_.get<std::string>("temperature key",
          domain+std::string("temperature"));
  e_key_ = plist_.get<std::string>("energy key", domain+std::string("energy"));
  wc_key_ = plist_.get<std::string>("water content key",
          domain+std::string("water_content"));
  cv_key_ = plist_.get<std::string>("cell volume key",
          domain+std::string("cell_volume"));
  Key poro_default;
  if (domain.size() == 0) {
    poro_default = "porosity";
  } else {
    poro_default = "";
  }
  poro_key_ = plist_.get<std::string>("porosity key", poro_default);

  // Process the parameter list for methods
  std::string precon_string = plist_.get<std::string>("preconditioner type", "none");
  if (precon_string == "none") {
    precon_type_ = PRECON_NONE;
  } else if (precon_string == "ewc") {
    precon_type_ = PRECON_EWC;
  } else if (precon_string == "smart ewc") {
    precon_type_ = PRECON_SMART_EWC;
  }

  // select the method used for nonlinear prediction
  std::string predictor_string = plist_.get<std::string>("predictor type", "none");
  if (predictor_string == "none") {
    predictor_type_ = PREDICTOR_NONE;
  } else if (predictor_string == "ewc") {
    predictor_type_ = PREDICTOR_EWC;
  } else if (predictor_string == "smart ewc") {
    predictor_type_ = PREDICTOR_SMART_EWC;
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


// -----------------------------------------------------------------------------
// Initialize any data required in setup.
// -----------------------------------------------------------------------------
void MPCDelegateEWC::initialize(const Teuchos::Ptr<State>& S) {
  // Create and initialize old stored data for previous steps.
  if (predictor_type_ == PREDICTOR_EWC || predictor_type_ == PREDICTOR_SMART_EWC) {
    const Epetra_MultiVector& wc = *S->GetFieldData("water_content")
        ->ViewComponent("cell",false);
    const Epetra_MultiVector& e = *S->GetFieldData("energy")
        ->ViewComponent("cell",false);

    wc_prev2_ = Teuchos::rcp(new Epetra_MultiVector(wc));
    e_prev2_ = Teuchos::rcp(new Epetra_MultiVector(e));
    wc_prev2_->PutScalar(0.);
    e_prev2_->PutScalar(0.);

    time_prev2_ = S->time();
  }

  // initialize the Jacobian
  if (precon_type_ == PRECON_EWC || precon_type_ == PRECON_SMART_EWC) {
    int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
    jac_.resize(ncells, WhetStone::Tensor(2,2));
  }

  // initialize the model, which grabs all needed models from state
  model_->InitializeModel(S);
}


// -----------------------------------------------------------------------------
// Set state pointers.
// -----------------------------------------------------------------------------
void MPCDelegateEWC::set_states(const Teuchos::RCP<const State>& S,
        const Teuchos::RCP<State>& S_inter,
        const Teuchos::RCP<State>& S_next) {
  S_= S;
  S_inter_ = S_inter;
  S_next_ = S_next;
}


// -----------------------------------------------------------------------------
// Save info from previous iterations if needed.
// -----------------------------------------------------------------------------
void MPCDelegateEWC::commit_state(double dt, const Teuchos::RCP<State>& S) {
  if ((predictor_type_ == PREDICTOR_EWC || predictor_type_ == PREDICTOR_SMART_EWC)
      && S_inter_ != Teuchos::null) {
    // stash water content and energy in S_work.
    *wc_prev2_ = *S_inter_->GetFieldData("water_content")->ViewComponent("cell",false);
    *e_prev2_ = *S_inter_->GetFieldData("energy")->ViewComponent("cell",false);
    time_prev2_ = S_inter_->time();
  }

  if (predictor_type_ == PREDICTOR_EWC || predictor_type_ == PREDICTOR_SMART_EWC ||
      precon_type_ == PRECON_EWC || precon_type_ == PRECON_SMART_EWC) {
    // set model's value scalars and check it is setup
    model_->UpdateModel(S.ptr());
  }
}


// -----------------------------------------------------------------------------
// Modify the prediction from linearization of the time integration.
// -----------------------------------------------------------------------------
bool MPCDelegateEWC::modify_predictor(double h, Teuchos::RCP<TreeVector> up) {
  bool modified = false;
  if (predictor_type_ == PREDICTOR_EWC) {
    if (S_->time() == 0.) { // this needs to be fixed!
      return false;
    } else {
      modified = modify_predictor_ewc_(h,up);
    }
  } else if (predictor_type_ == PREDICTOR_SMART_EWC) {
    if (S_->time() == 0.) { // this needs to be fixed!
      return false;
    } else {
      modified = modify_predictor_smart_ewc_(h,up);
    }
  }
  return modified;
}


// -----------------------------------------------------------------------------
// Update of the preconditioner
// -----------------------------------------------------------------------------
void MPCDelegateEWC::update_precon(double t, Teuchos::RCP<const TreeVector> up, double h) {
  if (precon_type_ == PRECON_EWC || precon_type_ == PRECON_SMART_EWC) {
    update_precon_ewc_(t,up,h);
  }
}


// -----------------------------------------------------------------------------
// Application of the preconditioner.
// -----------------------------------------------------------------------------
void MPCDelegateEWC::precon(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu) {
  if (precon_type_ == PRECON_EWC) {
    precon_ewc_(u,Pu);
  } else if (precon_type_ == PRECON_SMART_EWC) {
    precon_smart_ewc_(u,Pu);
  }
}


bool MPCDelegateEWC::modify_predictor_ewc_(double h, Teuchos::RCP<TreeVector> up) {
  // projected guesses for T and p
  Teuchos::RCP<CompositeVector> temp_guess = up->SubVector(1)->data();
  Epetra_MultiVector& temp_guess_c = *temp_guess->ViewComponent("cell",false);
  Teuchos::RCP<CompositeVector> pres_guess = up->SubVector(0)->data();
  Epetra_MultiVector& pres_guess_c = *pres_guess->ViewComponent("cell",false);

  // T, p at the previous step
  const Epetra_MultiVector& T1 = *S_inter_->GetFieldData("temperature")
      ->ViewComponent("cell",false);
  const Epetra_MultiVector& p1 = *S_inter_->GetFieldData("pressure")
      ->ViewComponent("cell",false);

  // project energy and water content
  double dt_next = S_next_->time() - S_inter_->time();
  double dt_prev = S_inter_->time() - time_prev2_;

  // -- get wc and energy data
  const Epetra_MultiVector& wc0 = *wc_prev2_;
  const Epetra_MultiVector& wc1 = *S_inter_->GetFieldData("water_content")
      ->ViewComponent("cell",false);
  Epetra_MultiVector& wc2 = *S_next_->GetFieldData("water_content", "water_content")
      ->ViewComponent("cell",false);

  const Epetra_MultiVector& e0 = *e_prev2_;
  const Epetra_MultiVector& e1 = *S_inter_->GetFieldData("energy")
      ->ViewComponent("cell",false);
  Epetra_MultiVector& e2 = *S_next_->GetFieldData("energy", "energy")
      ->ViewComponent("cell",false);

  // -- project
  wc2 = wc0;
  e2 = e0;
  double dt_ratio = (dt_next + dt_prev) / dt_prev;
  wc2.Update(dt_ratio, wc1, 1. - dt_ratio);
  e2.Update(dt_ratio, e1, 1. - dt_ratio);


  // -- extra data
  const Epetra_MultiVector& poro = *S_next_->GetFieldData("porosity")
      ->ViewComponent("cell",false);
  const Epetra_MultiVector& cv = *S_next_->GetFieldData("cell_volume")
      ->ViewComponent("cell",false);

  double wc_scale = 10.;
  double e_scale = 10000.;

  int ncells = wc0.MyLength();
  for (int c=0; c!=ncells; ++c) {
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


bool MPCDelegateEWC::modify_predictor_smart_ewc_(double h, Teuchos::RCP<TreeVector> up) {
  // projected guesses for T and p
  Teuchos::RCP<CompositeVector> temp_guess = up->SubVector(1)->data();
  Epetra_MultiVector& temp_guess_c = *temp_guess->ViewComponent("cell",false);
  Teuchos::RCP<CompositeVector> pres_guess = up->SubVector(0)->data();
  Epetra_MultiVector& pres_guess_c = *pres_guess->ViewComponent("cell",false);

  // T, p at the previous step
  const Epetra_MultiVector& T1 = *S_inter_->GetFieldData("temperature")
      ->ViewComponent("cell",false);
  const Epetra_MultiVector& p1 = *S_inter_->GetFieldData("pressure")
      ->ViewComponent("cell",false);

  // project energy and water content
  double dt_next = S_next_->time() - S_inter_->time();
  double dt_prev = S_inter_->time() - time_prev2_;

  // -- get wc and energy data
  const Epetra_MultiVector& wc0 = *wc_prev2_;
  const Epetra_MultiVector& wc1 = *S_inter_->GetFieldData("water_content")
      ->ViewComponent("cell",false);
  Epetra_MultiVector& wc2 = *S_next_->GetFieldData("water_content", "water_content")
      ->ViewComponent("cell",false);

  const Epetra_MultiVector& e0 = *e_prev2_;
  const Epetra_MultiVector& e1 = *S_inter_->GetFieldData("energy")
      ->ViewComponent("cell",false);
  Epetra_MultiVector& e2 = *S_next_->GetFieldData("energy", "energy")
      ->ViewComponent("cell",false);

  // -- project
  wc2 = wc0;
  e2 = e0;
  double dt_ratio = (dt_next + dt_prev) / dt_prev;
  wc2.Update(dt_ratio, wc1, 1. - dt_ratio);
  e2.Update(dt_ratio, e1, 1. - dt_ratio);


  // -- extra data
  const Epetra_MultiVector& poro = *S_next_->GetFieldData("porosity")
      ->ViewComponent("cell",false);
  const Epetra_MultiVector& cv = *S_next_->GetFieldData("cell_volume")
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
      ierr = model_->InverseEvaluate(e2[0][c]/cv[0][c], wc2[0][c]/cv[0][c],
              poro[0][c], T, p);

#if DEBUG_FLAG
      std::cout << "   p_ewc,T_ewc  = " << p << ", " << T << std::endl;
#endif

      if (!ierr) { // valid solution, no zero determinates, etc
#if DEBUG_FLAG
        std::cout << "   EWC Accepted" << std::endl;
#endif
        temp_guess_c[0][c] = T;
        pres_guess_c[0][c] = p;
      }
    }
  }
  return true;
}


void MPCDelegateEWC::precon_smart_ewc_(Teuchos::RCP<const TreeVector> u,
        Teuchos::RCP<TreeVector> Pu) {

}


void MPCDelegateEWC::precon_ewc_(Teuchos::RCP<const TreeVector> u,
        Teuchos::RCP<TreeVector> Pu) {

}


void MPCDelegateEWC::update_precon_ewc_(double t, Teuchos::RCP<const TreeVector> up, double h) {
  Key dedT_key = std::string("d")+e_key_+std::string("_d")+temp_key_;
  const Epetra_MultiVector& dedT = *S_next_->GetFieldData(dedT_key)
      ->ViewComponent("cell",false);

  Key dedp_key = std::string("d")+e_key_+std::string("_d")+pres_key_;
  const Epetra_MultiVector& dedp = *S_next_->GetFieldData(dedp_key)
      ->ViewComponent("cell",false);

  Key dwcdT_key = std::string("d")+wc_key_+std::string("_d")+temp_key_;
  const Epetra_MultiVector& dwcdT = *S_next_->GetFieldData(dwcdT_key)
      ->ViewComponent("cell",false);

  Key dwcdp_key = std::string("d")+wc_key_+std::string("_d")+pres_key_;
  const Epetra_MultiVector& dwcdp = *S_next_->GetFieldData(dwcdp_key)
      ->ViewComponent("cell",false);

  int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  for (int c=0; c!=ncells; ++c) {
    jac_[c](0,0) = dedT[0][c];
    jac_[c](0,1) = dedp[0][c];
    jac_[c](1,0) = dwcdT[0][c];
    jac_[c](1,1) = dwcdp[0][c];
  }
}

} // namespace
