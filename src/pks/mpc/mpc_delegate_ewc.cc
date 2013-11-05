/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Interface for EWC, a helper class that does projections and preconditioners in
energy/water-content space instead of temperature/pressure space.
------------------------------------------------------------------------- */
#include "global_verbosity.hh"
#include "FieldEvaluator.hh"
#include "ewc_model.hh"
#include "mpc_delegate_ewc.hh"

namespace Amanzi {

#define DEBUG_FLAG 1

// -----------------------------------------------------------------------------
// Constructor
// -----------------------------------------------------------------------------
MPCDelegateEWC::MPCDelegateEWC(Teuchos::ParameterList& plist) :
    plist_(Teuchos::rcpFromRef(plist)) {
  // set up the VerboseObject
  std::string name = plist_->get<std::string>("PK name")+std::string(" EWC");
  vo_ = Teuchos::rcp(new VerboseObject(name, *plist_));
}

// -----------------------------------------------------------------------------
// Allocate any data or models required.
// -----------------------------------------------------------------------------
void MPCDelegateEWC::setup(const Teuchos::Ptr<State>& S) {
  // Verbosity
  std::string name = plist_->get<std::string>("PK name")+std::string(" EWC");

  // Get the mesh
  Key domain = plist_->get<std::string>("domain key", "");
  if (domain.size() != 0) {
    mesh_ = S->GetMesh(domain);
  } else {
    mesh_ = S->GetMesh("domain");
  }

  // set up a debugger
  db_ = Teuchos::rcp(new Debugger(mesh_, name, *plist_));

  // Process the parameter list for data Keys
  if (domain.size() > 0) domain.append(1, '_');
  pres_key_ = plist_->get<std::string>("pressure key", domain+std::string("pressure"));
  temp_key_ = plist_->get<std::string>("temperature key",
          domain+std::string("temperature"));
  e_key_ = plist_->get<std::string>("energy key", domain+std::string("energy"));
  wc_key_ = plist_->get<std::string>("water content key",
          domain+std::string("water_content"));
  cv_key_ = plist_->get<std::string>("cell volume key",
          domain+std::string("cell_volume"));

  // Process the parameter list for methods
  std::string precon_string = plist_->get<std::string>("preconditioner type", "none");
  if (precon_string == "none") {
    precon_type_ = PRECON_NONE;
  } else if (precon_string == "ewc") {
    precon_type_ = PRECON_EWC;
  } else if (precon_string == "smart ewc") {
    precon_type_ = PRECON_SMART_EWC;
  }

  // select the method used for nonlinear prediction
  std::string predictor_string = plist_->get<std::string>("predictor type", "none");
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
    if (plist_->isParameter("cusp distance in T")) {
      cusp_size_T_freezing_ = plist_->get<double>("cusp distance in T");
      cusp_size_T_thawing_ = cusp_size_T_freezing_;
    } else {
      cusp_size_T_freezing_ = plist_->get<double>("cusp distance in T, freezing", 0.005);
      cusp_size_T_thawing_ = plist_->get<double>("cusp distance in T, thawing", 0.005);
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
}


// -----------------------------------------------------------------------------
// Modify the prediction from linearization of the time integration.
// -----------------------------------------------------------------------------
bool MPCDelegateEWC::modify_predictor(double h, Teuchos::RCP<TreeVector> up) {
  bool modified = false;
  double dt_prev = S_inter_->time() - time_prev2_;

  if (predictor_type_ == PREDICTOR_EWC) {
    if (dt_prev > 0.) {
      modified = modify_predictor_ewc_(h,up);
    } else {
      return false;
    }
  } else if (predictor_type_ == PREDICTOR_SMART_EWC) {
    if (dt_prev > 0.) {
      modified = modify_predictor_smart_ewc_(h,up);
    } else {
      return false;
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
  if ((precon_type_ == PRECON_EWC) || (precon_type_ == PRECON_SMART_EWC)) {
    precon_ewc_(u,Pu);
  }
}


bool MPCDelegateEWC::modify_predictor_ewc_(double h, Teuchos::RCP<TreeVector> up) {
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_HIGH))
    *vo_->os() << "  Modifying predictor using EWC algorithm" << std::endl;

  // projected guesses for T and p
  Teuchos::RCP<CompositeVector> temp_guess = up->SubVector(1)->Data();
  Epetra_MultiVector& temp_guess_c = *temp_guess->ViewComponent("cell",false);
  Teuchos::RCP<CompositeVector> pres_guess = up->SubVector(0)->Data();
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
  const Epetra_MultiVector& cv = *S_next_->GetFieldData(cv_key_)
      ->ViewComponent("cell",false);

  int ncells = wc0.MyLength();
  for (int c=0; c!=ncells; ++c) {
    AmanziGeometry::Point result(2);
    int ierr = 0;

    double p = p1[0][c];
    double T = T1[0][c];

    // uses intensive forms, so must divide by cell volume.
#if DEBUG_FLAG
    if (vo_->os_OK(Teuchos::VERB_EXTREME)) {
      int rank = mesh_->get_comm()->MyPID();
      Teuchos::RCP<VerboseObject> dcvo = db_->GetVerboseObject(c, rank);
      if (dcvo != Teuchos::null && dcvo->os_OK(Teuchos::VERB_EXTREME)) {
        Teuchos::OSTab tab = dcvo->getOSTab();
        *dcvo->os() << "Inverting: c = " << c << std::endl;
        *dcvo->os() << "   p,T  = " << p << ", " << T << std::endl;
        *dcvo->os() << "   wc,e = " << wc1[0][c] << ", " << e1[0][c] << std::endl;
        *dcvo->os() << "   goal = " << wc2[0][c] << ", " << e2[0][c] << std::endl;
      }
    }
#endif
    model_->UpdateModel(S_next_.ptr(), c);
    ierr = model_->InverseEvaluate(e2[0][c]/cv[0][c], wc2[0][c]/cv[0][c], T, p);

    if (!ierr) { // valid solution, no zero determinates, etc
      temp_guess_c[0][c] = T;
      pres_guess_c[0][c] = p;
    }
  }
  return true;
}


bool MPCDelegateEWC::modify_predictor_smart_ewc_(double h, Teuchos::RCP<TreeVector> up) {
  Teuchos::OSTab tab = vo_->getOSTab();
  // projected guesses for T and p
  Teuchos::RCP<CompositeVector> temp_guess = up->SubVector(1)->Data();
  Epetra_MultiVector& temp_guess_c = *temp_guess->ViewComponent("cell",false);
  Teuchos::RCP<CompositeVector> pres_guess = up->SubVector(0)->Data();
  Epetra_MultiVector& pres_guess_c = *pres_guess->ViewComponent("cell",false);

  if (vo_->os_OK(Teuchos::VERB_HIGH)) {
    *vo_->os() << "  Modifying predictor using SmartEWC algorithm" << std::endl;
    std::vector<std::string> vnames;
    vnames.push_back("p_extrap"); vnames.push_back("T_extrap");
    std::vector< Teuchos::Ptr<const CompositeVector> > vecs;
    vecs.push_back(pres_guess.ptr()); vecs.push_back(temp_guess.ptr());
    db_->WriteVectors(vnames, vecs, true);
  }

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
  const Epetra_MultiVector& cv = *S_next_->GetFieldData(cv_key_)
      ->ViewComponent("cell",false);

  int rank = mesh_->get_comm()->MyPID();
  int ncells = wc0.MyLength();
  for (int c=0; c!=ncells; ++c) {
    Teuchos::RCP<VerboseObject> dcvo = Teuchos::null;
    if (vo_->os_OK(Teuchos::VERB_EXTREME))
      dcvo = db_->GetVerboseObject(c, rank);
    Teuchos::OSTab dctab = dcvo == Teuchos::null ? vo_->getOSTab() : dcvo->getOSTab();

    AmanziGeometry::Point result(2);
    int ierr = 0;

    double T_guess = temp_guess_c[0][c];
    double T_prev = T1[0][c];
    double T_prev2 = (T_guess - dt_ratio*T_prev) / (1. - dt_ratio);

    double p = p1[0][c];
    double T = T1[0][c];

    double wc_tmp(0.), e_tmp(0.);
    if (dcvo != Teuchos::null && dcvo->os_OK(Teuchos::VERB_EXTREME))
      *dcvo->os() << "Predicting: c = " << c << std::endl
                  << "   based upon h_old = " << dt_prev << ", h_next = " << dt_next << std::endl
                  << "   -------------" << std::endl
                  << "   Prev wc,e: " << wc1[0][c] << ", " << e1[0][c] << std::endl
                  << "   Prev p,T: " << p << ", " << T << std::endl
                  << "   -------------" << std::endl
                  << "   Extrap wc,e: " << wc2[0][c] << ", " << e2[0][c] << std::endl
                  << "   Extrap p,T: " << pres_guess_c[0][c] << ", " << T_guess << std::endl;

    model_->UpdateModel(S_next_.ptr(), c);
    ierr = model_->Evaluate(T_guess, pres_guess_c[0][c],
                            e_tmp, wc_tmp);
    ASSERT(!ierr);

    if (dcvo != Teuchos::null && dcvo->os_OK(Teuchos::VERB_EXTREME))
      *dcvo->os() << "   Calc wc,e of extrap: " << wc_tmp*cv[0][c] << ", " << e_tmp*cv[0][c] << std::endl
                  << "   -------------" << std::endl;

    if (T_guess - T < 0.) {  // decreasing, freezing
      if (dcvo != Teuchos::null && dcvo->os_OK(Teuchos::VERB_EXTREME))
          *dcvo->os() << "   decreasing temps..." << std::endl;

      if (T_guess >= 273.149) {
        if (dcvo != Teuchos::null && dcvo->os_OK(Teuchos::VERB_EXTREME))
          *dcvo->os() << "   above freezing, keep T,p projections" << std::endl;
        // pass, guesses are good
      } else {
        // -- invert for T,p at the projected ewc
        ierr = model_->InverseEvaluate(e2[0][c]/cv[0][c], wc2[0][c]/cv[0][c], T, p);
        double s_g, s_l, s_i;
        ierr |= model_->EvaluateSaturations(T,p,s_g,s_l,s_i);
        if (ierr) {
          if (dcvo != Teuchos::null && dcvo->os_OK(Teuchos::VERB_EXTREME))
              *dcvo->os() << "FAILED EWC PREDICTOR" << std::endl;

        } else {
          if (dcvo != Teuchos::null && dcvo->os_OK(Teuchos::VERB_EXTREME))
            *dcvo->os() << "   saturations new ewc values (g,l,i): " << s_g << ", " << s_l << ", " << s_i << std::endl
                        << "   frozen_fraction at the new ewc values: " << (s_i/(s_l + s_i)) << std::endl;

          // Check if our energy projection is before the 2nd inflection point
          if ( s_i / (s_l + s_i) < 0.99) {
            if (dcvo != Teuchos::null && dcvo->os_OK(Teuchos::VERB_EXTREME))
              *dcvo->os() << "     kept within the transition zone." << std::endl
                          << "   p,T = " << p << ", " << T << std::endl;

            // in the transition zone of latent heat exchange
            if (T > 200.) {
              temp_guess_c[0][c] = T;
              pres_guess_c[0][c] = p;
            } else {
              if (dcvo != Teuchos::null && dcvo->os_OK(Teuchos::VERB_EXTREME))
                *dcvo->os() << "       not admissible!" << std::endl;
            }
          } else {
            // pass, past the transition and we are just chilling ice
            if (dcvo != Teuchos::null && dcvo->os_OK(Teuchos::VERB_EXTREME))
              *dcvo->os() << "     outside the transition zone." << std::endl;
          }
        }
      }

    } else { // increasing, thawing
      if (dcvo != Teuchos::null && dcvo->os_OK(Teuchos::VERB_EXTREME))
        *dcvo->os() << "   increasing temps..." << std::endl;

      if (T >= 273.149) {
        if (dcvo != Teuchos::null && dcvo->os_OK(Teuchos::VERB_EXTREME))
          *dcvo->os() << "   above freezing, keep T,p projections" << std::endl;

        // pass, guesses are good
      } else {
        double s_g, s_l, s_i;
        ierr = model_->EvaluateSaturations(temp_guess_c[0][c],pres_guess_c[0][c],
                s_g,s_l,s_i);
        ASSERT(!ierr);
        if (dcvo != Teuchos::null && dcvo->os_OK(Teuchos::VERB_EXTREME))
            *dcvo->os() << "   saturations extrap pT values (g,l,i): " << s_g << ", " << s_l << ", " << s_i << std::endl
                        << "   frozen_fraction at the new pT values: " << (s_i/(s_l + s_i)) << std::endl;

        if (( s_i / (s_l + s_i) > 0.99)) {
          // pass, warming ice but still outside of the transition region
          if (dcvo != Teuchos::null && dcvo->os_OK(Teuchos::VERB_EXTREME))
              *dcvo->os() << "     outside the transition zone." << std::endl;
        } else {
          // in the transition zone of latent heat exchange
          ierr = model_->InverseEvaluate(e2[0][c]/cv[0][c], wc2[0][c]/cv[0][c], T, p);
          if (dcvo != Teuchos::null && dcvo->os_OK(Teuchos::VERB_EXTREME))
              *dcvo->os() << "     kept within the transition zone." << std::endl
                          << "   p,T = " << p << ", " << T << std::endl;

          if (!ierr && T < 273.149) {
            temp_guess_c[0][c] = T;
            pres_guess_c[0][c] = p;
          } else {
            // pass, past the transition zone and onto the upper branch
          }
        }
      }
    }
  }
  return true;
}



bool MPCDelegateEWC::modify_predictor_energy_(double h, Teuchos::RCP<TreeVector> up) {
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_HIGH))
    *vo_->os() << "  Modifying predictor using SmartEWC algorithm" << std::endl;

  // projected guesses for T and p
  Teuchos::RCP<CompositeVector> temp_guess = up->SubVector(1)->Data();
  Epetra_MultiVector& temp_guess_c = *temp_guess->ViewComponent("cell",false);
  Teuchos::RCP<CompositeVector> pres_guess = up->SubVector(0)->Data();
  Epetra_MultiVector& pres_guess_c = *pres_guess->ViewComponent("cell",false);

  // T, p at the previous step
  const Epetra_MultiVector& T1 = *S_inter_->GetFieldData("temperature")
      ->ViewComponent("cell",false);

  // project energy
  double dt_next = S_next_->time() - S_inter_->time();
  double dt_prev = S_inter_->time() - time_prev2_;

  // -- get energy data
  const Epetra_MultiVector& e0 = *e_prev2_;
  const Epetra_MultiVector& e1 = *S_inter_->GetFieldData("energy")
      ->ViewComponent("cell",false);
  Epetra_MultiVector& e2 = *S_next_->GetFieldData("energy", "energy")
      ->ViewComponent("cell",false);

  // -- project
  e2 = e0;
  double dt_ratio = (dt_next + dt_prev) / dt_prev;
  e2.Update(dt_ratio, e1, 1. - dt_ratio);


  // -- extra data
  const Epetra_MultiVector& cv = *S_next_->GetFieldData(cv_key_)
      ->ViewComponent("cell",false);

  double e_scale = 10000.;

  int rank = mesh_->get_comm()->MyPID();
  int ncells = e0.MyLength();
  for (int c=0; c!=ncells; ++c) {
    Teuchos::RCP<VerboseObject> dcvo = Teuchos::null;
    if (vo_->os_OK(Teuchos::VERB_EXTREME))
      dcvo = db_->GetVerboseObject(c, rank);
    Teuchos::OSTab dctab = dcvo == Teuchos::null ? vo_->getOSTab() : dcvo->getOSTab();

    bool ewc_predictor = false;
    AmanziGeometry::Point result(2);
    int ierr = 0;

    // Determine whether to use T/p or E/WC as the projection
    double T_prev = T1[0][c];
    double T_guess = temp_guess_c[0][c];
    double T_prev2 = (T_guess - dt_ratio*T_prev) / (1. - dt_ratio);

    double T = T1[0][c];

    double wc_tmp(0.), e_tmp(0.);
    if (dcvo != Teuchos::null && dcvo->os_OK(Teuchos::VERB_EXTREME))
      *dcvo->os() << "Inverting: c = " << c << std::endl
                  << "   based upon h_old = " << dt_prev << ", h_next = " << dt_next << std::endl
                  << "   -------------" << std::endl
                  << "   2Prev e: " << e0[0][c] << std::endl
                  << "   2Prev T: " << T_prev2 << std::endl
                  << "   -------------" << std::endl
                  << "   Prev e: " << e1[0][c] << std::endl
                  << "   Prev T: " << T << std::endl
                  << "   -------------" << std::endl
                  << "   Extrap p,T: " << pres_guess_c[0][c] << ", " << temp_guess_c[0][c] << std::endl;

    model_->UpdateModel(S_next_.ptr(), c);
    ierr = model_->Evaluate(temp_guess_c[0][c], pres_guess_c[0][c],
                            e_tmp, wc_tmp);
    ASSERT(!ierr);

    if (dcvo != Teuchos::null && dcvo->os_OK(Teuchos::VERB_EXTREME))
      *dcvo->os() << "   Calc e of extrap: " << e_tmp*cv[0][c] << std::endl
                  << "   Extrap e: " << e2[0][c] << std::endl;


    if (std::abs(e_tmp - e1[0][c]) > std::abs(e2[0][c] - e1[0][c])) { // extrap energy change is smaller
      // uses intensive forms, so must divide by cell volume.
      ierr = model_->InverseEvaluateEnergy(e2[0][c]/cv[0][c], pres_guess_c[0][c], T);

      if (dcvo != Teuchos::null && dcvo->os_OK(Teuchos::VERB_EXTREME))
        *dcvo->os() << "   Inverted T = " << T << std::endl;

      if (!ierr) { // valid solution, no zero determinates, etc
        if (std::abs(temp_guess_c[0][c] - T1[0][c]) > std::abs(T - T1[0][c])) {
          if (dcvo != Teuchos::null && dcvo->os_OK(Teuchos::VERB_EXTREME))
            *dcvo->os() << "   EWC Accepted" << std::endl;
          temp_guess_c[0][c] = T;
        } else {
          if (dcvo != Teuchos::null && dcvo->os_OK(Teuchos::VERB_EXTREME))
            *dcvo->os() << "   EWC NOT Accepted" << std::endl;
        }
      }
    } else {
      if (dcvo != Teuchos::null && dcvo->os_OK(Teuchos::VERB_EXTREME))
        *dcvo->os() << "   EWC energy change not smaller!" << std::endl;
    }
  }
  return true;
}



bool MPCDelegateEWC::modify_predictor_smart_energy_(double h, Teuchos::RCP<TreeVector> up) {
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_HIGH))
    *vo_->os() << "  Modifying predictor using SmartEWC algorithm" << std::endl;

  // projected guesses for T and p
  Teuchos::RCP<CompositeVector> temp_guess = up->SubVector(1)->Data();
  Epetra_MultiVector& temp_guess_c = *temp_guess->ViewComponent("cell",false);
  Teuchos::RCP<CompositeVector> pres_guess = up->SubVector(0)->Data();
  Epetra_MultiVector& pres_guess_c = *pres_guess->ViewComponent("cell",false);

  // T, p at the previous step
  const Epetra_MultiVector& T1 = *S_inter_->GetFieldData("temperature")
      ->ViewComponent("cell",false);

  // project energy
  double dt_next = S_next_->time() - S_inter_->time();
  double dt_prev = S_inter_->time() - time_prev2_;

  // -- get energy data
  const Epetra_MultiVector& e0 = *e_prev2_;
  const Epetra_MultiVector& e1 = *S_inter_->GetFieldData("energy")
      ->ViewComponent("cell",false);
  Epetra_MultiVector& e2 = *S_next_->GetFieldData("energy", "energy")
      ->ViewComponent("cell",false);

  // -- project
  e2 = e0;
  double dt_ratio = (dt_next + dt_prev) / dt_prev;
  e2.Update(dt_ratio, e1, 1. - dt_ratio);


  // -- extra data
  const Epetra_MultiVector& cv = *S_next_->GetFieldData(cv_key_)
      ->ViewComponent("cell",false);

  double e_scale = 10000.;

  int rank = mesh_->get_comm()->MyPID();
  int ncells = e0.MyLength();
  for (int c=0; c!=ncells; ++c) {
    // debugger
    Teuchos::RCP<VerboseObject> dcvo = Teuchos::null;
    if (vo_->os_OK(Teuchos::VERB_EXTREME))
      dcvo = db_->GetVerboseObject(c, rank);
    Teuchos::OSTab dctab = dcvo == Teuchos::null ? vo_->getOSTab() : dcvo->getOSTab();

    bool ewc_predictor = false;
    AmanziGeometry::Point result(2);
    int ierr = 0;

    // Determine whether to use T/p or E/WC as the projection
    double T_prev = T1[0][c];
    double T_guess = temp_guess_c[0][c];
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
      double T = T1[0][c];

      double wc_tmp(0.), e_tmp(0.);
      if (dcvo != Teuchos::null && dcvo->os_OK(Teuchos::VERB_EXTREME))
        *dcvo->os() << "Inverting: c = " << c << std::endl
                    << "   based upon h_old = " << dt_prev << ", h_next = " << dt_next << std::endl
                    << "   -------------" << std::endl
                    << "   2Prev e: " << e0[0][c] << std::endl
                    << "   2Prev T: " << T_prev2 << std::endl
                    << "   -------------" << std::endl
                    << "   Prev e: " << e1[0][c] << std::endl
                    << "   Prev T: " << T << std::endl
                    << "   -------------" << std::endl
                    << "   Extrap p,T: " << pres_guess_c[0][c] << ", " << temp_guess_c[0][c] << std::endl;

      model_->UpdateModel(S_next_.ptr(), c);
      ierr = model_->Evaluate(temp_guess_c[0][c], pres_guess_c[0][c],
              e_tmp, wc_tmp);
      ASSERT(!ierr);

      if (dcvo != Teuchos::null && dcvo->os_OK(Teuchos::VERB_EXTREME))
        *dcvo->os() << "   Calc e of extrap: " << e_tmp*cv[0][c] << std::endl
                    << "   Extrap e: " << e2[0][c] << std::endl;


      if (std::abs(e_tmp - e1[0][c]) > std::abs(e2[0][c] - e1[0][c])) { // extrap energy change is smaller
        // uses intensive forms, so must divide by cell volume.
        ierr = model_->InverseEvaluateEnergy(e2[0][c]/cv[0][c], pres_guess_c[0][c], T);

        if (dcvo != Teuchos::null && dcvo->os_OK(Teuchos::VERB_EXTREME))
          *dcvo->os() << "   Inverted T = " << T << std::endl;

        if (!ierr) { // valid solution, no zero determinates, etc
          if (std::abs(temp_guess_c[0][c] - T1[0][c]) > std::abs(T - T1[0][c])) {
            if (dcvo != Teuchos::null && dcvo->os_OK(Teuchos::VERB_EXTREME))
              *dcvo->os() << "   EWC Accepted" << std::endl;
            temp_guess_c[0][c] = T;
          } else {
            if (dcvo != Teuchos::null && dcvo->os_OK(Teuchos::VERB_EXTREME))
              *dcvo->os() << "   EWC NOT Accepted" << std::endl;
          }
        }
      } else {
        if (dcvo != Teuchos::null && dcvo->os_OK(Teuchos::VERB_EXTREME))
          *dcvo->os() << "   EWC energy change not smaller!" << std::endl;
      }
    }
  }
  return true;
}


void MPCDelegateEWC::precon_ewc_(Teuchos::RCP<const TreeVector> u,
        Teuchos::RCP<TreeVector> Pu) {
  // Pu currently stores (dp_std,dT_std).  The approach here is:
  // -- calculate (wc_ewc,e_ewc), a wc and energy-basis correction:
  //        wc_ewc = wc_0 - [ dwc/dp  dwc/dT ]( dp_std )
  //         e_ewc =  e_0 - [  de/dp   de/dT ]( dT_std )
  // -- invert the relation, calculating (p_ewc, T_ewc) = ewc^-1(wc_ewc, e_ewc)
  // -- use the smaller update of ||(p_ewc,T_ewc) - (p_0,T_0)|| and ||(dp_std, dT_std)||.

  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_HIGH))
    *vo_->os() << "  Preconditioning using SmartEWC algorithm" << std::endl;

  // projected guesses for T and p
  Epetra_MultiVector& dp_std = *Pu->SubVector(0)->Data()->ViewComponent("cell",false);
  Epetra_MultiVector& dT_std = *Pu->SubVector(1)->Data()->ViewComponent("cell",false);

  // additional data required
  const Epetra_MultiVector& cv = *S_next_->GetFieldData("cell_volume")->ViewComponent("cell",false);

  // old values
  const Epetra_MultiVector& p_old = *S_next_->GetFieldData("pressure")->ViewComponent("cell",false);
  const Epetra_MultiVector& T_old = *S_next_->GetFieldData("temperature")->ViewComponent("cell",false);
  const Epetra_MultiVector& wc_old = *S_next_->GetFieldData("water_content")->ViewComponent("cell",false);
  const Epetra_MultiVector& e_old = *S_next_->GetFieldData("energy")->ViewComponent("cell",false);

  // min change values... ewc is not useful near convergence
  double dT_min = 0.01;
  double dp_min = 100.;

  int rank = mesh_->get_comm()->MyPID();
  int ncells = cv.MyLength();
  for (int c=0; c!=ncells; ++c) {

    // debugger
    Teuchos::RCP<VerboseObject> dcvo = Teuchos::null;
    if (vo_->os_OK(Teuchos::VERB_EXTREME))
      dcvo = db_->GetVerboseObject(c, rank);
    Teuchos::OSTab dctab = dcvo == Teuchos::null ? vo_->getOSTab() : dcvo->getOSTab();

    double T_prev = T_old[0][c];
    double T_std = T_prev - dT_std[0][c];
    bool precon_ewc = false;

    if (dcvo != Teuchos::null && dcvo->os_OK(Teuchos::VERB_EXTREME))
      *dcvo->os() << "Precon: c = " << c << std::endl;

    // only do EWC if corrections are large, ie not clearly converging
    model_->UpdateModel(S_next_.ptr(), c);
    if (std::abs(dT_std[0][c]) > dT_min || std::abs(dp_std[0][c]) > dp_min) {
      if (-dT_std[0][c] < 0.) {  // decreasing, freezing
        if (dcvo != Teuchos::null && dcvo->os_OK(Teuchos::VERB_EXTREME))
          *dcvo->os() << "   decreasing temps..." << std::endl;

        if (T_std >= 273.15) {
          if (dcvo != Teuchos::null && dcvo->os_OK(Teuchos::VERB_EXTREME))
            *dcvo->os() << "   above freezing, keep std correction" << std::endl;
          // pass, guesses are good

        } else {
          // calculate the correction in ewc
          double wc_ewc = wc_old[0][c]
              - (jac_[c](0,0) * dp_std[0][c] + jac_[c](0,1) * dT_std[0][c]);
          double e_ewc = e_old[0][c]
              - (jac_[c](1,0) * dp_std[0][c] + jac_[c](1,1) * dT_std[0][c]);
          if (dcvo != Teuchos::null && dcvo->os_OK(Teuchos::VERB_EXTREME))
            *dcvo->os() << "   Prev p,T: " << p_old[0][c] << ", " << T_old[0][c] << std::endl
                        << "   Prev wc,e: " << wc_old[0][c] << ", " << e_old[0][c] << std::endl
                        << "   -------------" << std::endl
                        << "   Std correction dp,dT: " << dp_std[0][c] << ", " << dT_std[0][c] << std::endl
                        << "   applying the EWC precon:" << std::endl
                        << "     wc,e_ewc = " << wc_ewc << ", " << e_ewc << std::endl;

          // -- invert for T,p at the projected ewc
          double T(T_prev), p(p_old[0][c]);
          int ierr = model_->InverseEvaluate(e_ewc/cv[0][c], wc_ewc/cv[0][c],
                  T, p);
          double s_g, s_l, s_i;
          ierr |= model_->EvaluateSaturations(T,p,s_g,s_l,s_i);

          if (ierr) {
            if (dcvo != Teuchos::null && dcvo->os_OK(Teuchos::VERB_EXTREME))
              *dcvo->os() << "FAILED EWC PRECON" << std::endl;

          } else {
            if (dcvo != Teuchos::null && dcvo->os_OK(Teuchos::VERB_EXTREME))
              *dcvo->os() << "   frozen_fraction at the new ewc corrected values: " << (s_i/(s_l + s_i)) << std::endl;

            // Check if our energy projection is before the 2nd inflection point
            if ( s_i / (s_l + s_i) < 0.99) {
              double dT_ewc = T_old[0][c] - T;
              double dp_ewc = p_old[0][c] - p;

              if (dcvo != Teuchos::null && dcvo->os_OK(Teuchos::VERB_EXTREME)) {
                *dcvo->os() << "     kept within the transition zone." << std::endl
                            << "   p,T_ewc = " << p << ", " << T << std::endl
                            << "   dp,dT_ewc = " << dp_ewc << ", " << dT_ewc << std::endl;
                if (std::abs(dT_ewc) > dT_min || std::abs(dp_ewc) > dp_min) {
                  *dcvo->os() << "  sufficient change" << std::endl;
                } else {
                  *dcvo->os() << "  insufficient change" << std::endl;
                }
              }

              if (std::abs(dT_ewc) > dT_min || std::abs(dp_ewc) > dp_min) {
                dT_std[0][c] = dT_ewc;
                dp_std[0][c] = dp_ewc;
              }
            } else {
              // pass, past the transition and we are just chilling ice
              if (dcvo != Teuchos::null && dcvo->os_OK(Teuchos::VERB_EXTREME))
                *dcvo->os() << "     outside the transition zone." << std::endl;
            }
          }
        }
      } else { // increasing, thawing
        if (dcvo != Teuchos::null && dcvo->os_OK(Teuchos::VERB_EXTREME))
          *dcvo->os() << "   increasing temps..." << std::endl;

        if (T_std >= 273.15) {
          // pass, update is are good
        } else {
          double s_g, s_l, s_i;
          double p_std = p_old[0][c] - dT_std[0][c];
          int ierr = model_->EvaluateSaturations(T_std,p_std,s_g,s_l,s_i);
          ASSERT(!ierr);

          if (( s_i / (s_l + s_i) > 0.99)) {
            // pass, warming ice but still outside of the transition region
          } else {
            // calculate the correction in ewc
            double wc_ewc = wc_old[0][c]
                - (jac_[c](0,0) * dp_std[0][c] + jac_[c](0,1) * dT_std[0][c]);
            double e_ewc = e_old[0][c]
                - (jac_[c](1,0) * dp_std[0][c] + jac_[c](1,1) * dT_std[0][c]);
            if (dcvo != Teuchos::null && dcvo->os_OK(Teuchos::VERB_EXTREME))
              *dcvo->os() << "   Prev p,T: " << p_old[0][c] << ", " << T_old[0][c] << std::endl
                          << "   Prev wc,e: " << wc_old[0][c] << ", " << e_old[0][c] << std::endl
                          << "   -------------" << std::endl
                          << "   Std correction dp,dT: " << dp_std[0][c] << ", " << dT_std[0][c] << std::endl
                          << "   applying the EWC precon:" << std::endl
                          << "     wc,e_ewc = " << wc_ewc << ", " << e_ewc << std::endl;

            // -- invert for T,p at the projected ewc
            double T(T_prev), p(p_old[0][c]);
            int ierr = model_->InverseEvaluate(e_ewc/cv[0][c], wc_ewc/cv[0][c],
                    T, p);

            if (!ierr && T < 273.15) {
              double dT_ewc = T_old[0][c] - T;
              double dp_ewc = p_old[0][c] - p;

              if (dcvo != Teuchos::null && dcvo->os_OK(Teuchos::VERB_EXTREME)) {
                *dcvo->os() << "     kept within the transition zone." << std::endl
                            << "   p,T_ewc = " << p << ", " << T << std::endl
                            << "   dp,dT_ewc = " << dp_ewc << ", " << dT_ewc << std::endl;
                if (std::abs(dT_ewc) > dT_min || std::abs(dp_ewc) > dp_min) {
                  *dcvo->os() << "  sufficient change" << std::endl;
                } else {
                  *dcvo->os() << "  insufficient change" << std::endl;
                }
              }

              if (std::abs(dT_ewc) > dT_min || std::abs(dp_ewc) > dp_min) {
                dT_std[0][c] = dT_ewc;
                dp_std[0][c] = dp_ewc;
              }

            } else {
              // pass, past the transition zone and onto the upper branch
              if (dcvo != Teuchos::null && dcvo->os_OK(Teuchos::VERB_EXTREME))
                *dcvo->os() << "     outside the transition zone." << std::endl;
            }
          }
        }
      }
    }
  }
}



void MPCDelegateEWC::update_precon_ewc_(double t, Teuchos::RCP<const TreeVector> up, double h) {
  Key dedT_key = std::string("d")+e_key_+std::string("_d")+temp_key_;
  S_next_->GetFieldEvaluator(e_key_)
      ->HasFieldDerivativeChanged(S_next_.ptr(), "ewc", temp_key_);
  const Epetra_MultiVector& dedT = *S_next_->GetFieldData(dedT_key)
      ->ViewComponent("cell",false);

  Key dedp_key = std::string("d")+e_key_+std::string("_d")+pres_key_;
  S_next_->GetFieldEvaluator(e_key_)
      ->HasFieldDerivativeChanged(S_next_.ptr(), "ewc", pres_key_);
  const Epetra_MultiVector& dedp = *S_next_->GetFieldData(dedp_key)
      ->ViewComponent("cell",false);

  Key dwcdT_key = std::string("d")+wc_key_+std::string("_d")+temp_key_;
  S_next_->GetFieldEvaluator(wc_key_)
      ->HasFieldDerivativeChanged(S_next_.ptr(), "ewc", temp_key_);
  const Epetra_MultiVector& dwcdT = *S_next_->GetFieldData(dwcdT_key)
      ->ViewComponent("cell",false);

  Key dwcdp_key = std::string("d")+wc_key_+std::string("_d")+pres_key_;
  S_next_->GetFieldEvaluator(wc_key_)
      ->HasFieldDerivativeChanged(S_next_.ptr(), "ewc", pres_key_);
  const Epetra_MultiVector& dwcdp = *S_next_->GetFieldData(dwcdp_key)
      ->ViewComponent("cell",false);

  int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  for (int c=0; c!=ncells; ++c) {
    jac_[c](0,0) = dwcdp[0][c];
    jac_[c](0,1) = dwcdT[0][c];
    jac_[c](1,0) = dedp[0][c];
    jac_[c](1,1) = dedT[0][c];
  }
}

} // namespace
