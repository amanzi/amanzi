/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Default base with default implementations of methods for a PK integrated using
BDF.
------------------------------------------------------------------------- */

#include "Teuchos_TimeMonitor.hpp"
#include "BDF1_TI.hh"
#include "pk_bdf_base.hh"

namespace Amanzi {


// -----------------------------------------------------------------------------
// Setup
// -----------------------------------------------------------------------------
void PKBDFBase::setup(const Teuchos::Ptr<State>& S) {
  PKDefaultBase::setup(S);

  // initial timestep
  dt_ = plist_.get<double>("initial time step", 1.);

  // precon assembly
  assemble_preconditioner_ = plist_.get<bool>("assemble preconditioner", true);
};


// -----------------------------------------------------------------------------
// Initialization of timestepper.
// -----------------------------------------------------------------------------
void PKBDFBase::initialize(const Teuchos::Ptr<State>& S) {
  // set up the timestepping algorithm
  if (!plist_.get<bool>("strongly coupled PK", false)) {
    // -- get the timestepper plist and grab backtracking control
    Teuchos::ParameterList bdf_plist = plist_.sublist("time integrator");

    // NOT CURRENTLY USED, would need to be added to modify_corrector()
    backtracking_iterations_ = bdf_plist.get<int>("max backtrack count",0);
    backtracking_ = (backtracking_iterations_ > 0);

    // -- instantiate time stepper
    bdf_plist.set("initial time", S->time());
    time_stepper_ = Teuchos::rcp(new BDF1_TI<TreeVector>(*this, bdf_plist, solution_));

    // -- initialize time derivative
    Teuchos::RCP<TreeVector> solution_dot = Teuchos::rcp(new TreeVector(*solution_));
    solution_dot->PutScalar(0.0);

    // -- set initial state
    time_stepper_->set_initial_state(S->time(), solution_, solution_dot);
  }

  // set up the wallclock timer
  //  step_walltime_ = Teuchos::TimeMonitor::getNewCounter(name_+std::string(" BDF step timer"));
};


// -----------------------------------------------------------------------------
// Initialization of timestepper.
// -----------------------------------------------------------------------------
double PKBDFBase::get_dt() { return dt_; }


// -----------------------------------------------------------------------------
// Apply the preconditioner (default application).
// -----------------------------------------------------------------------------
void PKBDFBase::precon(Teuchos::RCP<const TreeVector> u,
        Teuchos::RCP<TreeVector> Pu) {
  preconditioner_->ApplyInverse(*u, Pu.ptr());
}


// -----------------------------------------------------------------------------
// Advance from state S to state S_next at time S.time + dt.
// -----------------------------------------------------------------------------
bool PKBDFBase::advance(double dt) {
  state_to_solution(S_next_, solution_);
  backtracking_count_ = 0;

  // take a bdf timestep
  double dt_solver;
  bool fail;
  if (true) { // this is here simply to create a context for timer,
              // which stops the clock when it is destroyed at the
              // closing brace.
    //    Teuchos::TimeMonitor timer(*step_walltime_);
    fail = time_stepper_->time_step(dt, dt_solver, solution_);
  }

  // VerboseObject stuff.
  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_HIGH, true)) {
    Teuchos::OSTab tab = getOSTab();
    //    Teuchos::TimeMonitor::summarize();
  }

  if (!fail) {
    // commit the step as successful
    time_stepper_->commit_solution(dt, solution_);
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


} // namespace
