/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Default base with default implementations of methods for a PK integrated using
BDF.
------------------------------------------------------------------------- */

#include "Teuchos_TimeMonitor.hpp"
#include "bdf1_time_integrator.hh"
#include "pk_bdf_base.hh"

namespace Amanzi {


// -----------------------------------------------------------------------------
// Setup
// -----------------------------------------------------------------------------
void PKBDFBase::setup(const Teuchos::Ptr<State>& S) {
  PKDefaultBase::setup(S);

  // initial timestep
  dt_ = plist_.get<double>("initial time step", 1.);
};


// -----------------------------------------------------------------------------
// Initialization of timestepper.
// -----------------------------------------------------------------------------
void PKBDFBase::initialize(const Teuchos::Ptr<State>& S) {
  // set up the timestepping algorithm
  if (!plist_.get<bool>("strongly coupled PK", false)) {
    // -- get the timestepper plist and grab what we need
    Teuchos::RCP<Teuchos::ParameterList> bdf_plist_p =
      Teuchos::rcp(new Teuchos::ParameterList(plist_.sublist("time integrator")));
    backtracking_ = (bdf_plist_p->get<int>("max backtrack count",0) > 0);

    // -- instantiate time stepper
    time_stepper_ = Teuchos::rcp(new BDF1TimeIntegrator(this, bdf_plist_p, solution_));

    // -- initialize time derivative
    Teuchos::RCP<TreeVector> solution_dot = Teuchos::rcp(new TreeVector(*solution_));
    solution_dot->PutScalar(0.0);

    // -- set initial state
    time_stepper_->set_initial_state(S->time(), solution_, solution_dot);
  }

  // set up the wallclock timer
  step_walltime_ = Teuchos::TimeMonitor::getNewCounter(name_+std::string(" BDF step timer"));
};


// -----------------------------------------------------------------------------
// Initialization of timestepper.
// -----------------------------------------------------------------------------
double PKBDFBase::get_dt() { return dt_; }


// -----------------------------------------------------------------------------
// Advance from state S to state S_next at time S.time + dt.
// -----------------------------------------------------------------------------
bool PKBDFBase::advance(double dt) {
  state_to_solution(S_next_, solution_);

  // take a bdf timestep
  double dt_solver;
  bool fail;
  if (true) { // this is here simply to create a context for timer,
              // which stops the clock when it is destroyed at the
              // closing brace.
    Teuchos::TimeMonitor timer(*step_walltime_);
    fail = time_stepper_->time_step(dt, dt_solver, solution_);
  }

  // VerboseObject stuff.
  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_HIGH, true)) {
    Teuchos::OSTab tab = getOSTab();
    Teuchos::TimeMonitor::summarize();
  }

  if (!fail) {
    // commit the step as successful
    time_stepper_->commit_solution(dt, solution_);
    solution_to_state(solution_, S_next_);
    commit_state(dt, S_next_);

    // update the timestep size
    if (dt_solver < dt_ && dt_solver >= dt) {
      // We took a smaller step than we recommended, and it worked fine (not
      // suprisingly).  Likely this was due to constraints from other PKs or
      // vis.  Do not reduce our recommendation.
    } else {
      dt_ = dt_solver;
    }
  } else {
    // take the decreased timestep size
    dt_ = dt_solver;
  }

  return fail;
};


// -----------------------------------------------------------------------------
// Allows a PK to reject solutions, which forces a timestep cut, or invokes
// backtracking.
// -----------------------------------------------------------------------------
bool PKBDFBase::is_admissible(Teuchos::RCP<const TreeVector> up) {
  if (!backtracking_) return true;

  // const issues with BDF1 need cleaned up...
  Teuchos::RCP<TreeVector> up_nc = Teuchos::rcp_const_cast<TreeVector>(up);

  // evaluate the residual
  Teuchos::RCP<TreeVector> g = Teuchos::rcp(new TreeVector(*up));
  fun(S_inter_->time(), S_next_->time(), Teuchos::null, up_nc, g);
  double norm(0.);
  int ierr = g->NormInf(&norm);

  // ensure the new residual is smaller than the old residual
  if (ierr || norm > residual_norm_) return false;

  residual_norm_ = norm;
  return true;
}


// -----------------------------------------------------------------------------
// Allows a PK to modify the initial guess.
// -----------------------------------------------------------------------------
bool PKBDFBase::modify_predictor(double h, Teuchos::RCP<TreeVector> up) {
  if (backtracking_) {
    // evaluate the residual
    Teuchos::RCP<TreeVector> g = Teuchos::rcp(new TreeVector(*up));
    fun(S_inter_->time(), S_next_->time(), Teuchos::null, up, g);
    g->NormInf(&residual_norm_);
  }
  return false;
}


} // namespace
