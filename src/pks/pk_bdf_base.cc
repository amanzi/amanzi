/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Default base with default implementations of methods for a PK integrated using
BDF.
------------------------------------------------------------------------- */

#include "bdf1_time_integrator.hh"
#include "pk_bdf_base.hh"

namespace Amanzi {


// -----------------------------------------------------------------------------
// Setup
// -----------------------------------------------------------------------------
void PKBDFBase::setup(const Teuchos::Ptr<State>& S) {
  PKDefaultBase::setup(S);

  // timestep evolution
  dt_ = plist_.get<double>("initial time step", 1.);
  time_step_reduction_factor_ =
    plist_.get<double>("time step reduction factor", 1.);

  // the default is no backtracking
  max_backtrack_count_ = plist_.get<int>("max backtrack count", 0);
  backtrack_damping_ = plist_.get<double>("backtrack damping", 1.0);

  // debugging option -- The BDF advance() method catches all errors, and then
  // re-throws ones it doesn't recognize.  Putting "false" in the PK PList
  // removes the catch, making debugging and debugger usage easier.
  catch_errors_ = plist_.get<bool>("catch thrown errors", true);

};


// -----------------------------------------------------------------------------
// Initialization of timestepper.
// -----------------------------------------------------------------------------
void PKBDFBase::initialize(const Teuchos::Ptr<State>& S) {
  // set up the timestepping algorithm
  if (!plist_.get<bool>("strongly coupled PK", false)) {
    // -- instantiate time stepper
    Teuchos::RCP<Teuchos::ParameterList> bdf_plist_p =
      Teuchos::rcp(new Teuchos::ParameterList(plist_.sublist("time integrator")));
    time_stepper_ = Teuchos::rcp(new BDF1TimeIntegrator(this, bdf_plist_p, solution_));

    // -- initialize time derivative
    Teuchos::RCP<TreeVector> solution_dot = Teuchos::rcp(new TreeVector(*solution_));
    solution_dot->PutScalar(0.0);

    // -- set initial state
    time_stepper_->set_initial_state(S->time(), solution_, solution_dot);
  }
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
  double h = dt;
  double dt_solver;
  if (catch_errors_) {
    try {
      dt_solver = time_stepper_->time_step(h, solution_);
    } catch (Exceptions::Amanzi_exception &error) {
      // fix me! --etc
      std::cout << "Timestepper called error: " << error.what() << std::endl;

      if (error.what() == std::string("BDF time step failed") ||
          error.what() == std::string("Cut time step")) {
        // try cutting the timestep
        dt_ = h*time_step_reduction_factor_;
        return true;
      } else {
        throw error;
      }
    }
  } else {
    dt_solver = time_stepper_->time_step(h, solution_);
  }

  // commit the step as successful
  time_stepper_->commit_solution(h, solution_);
  solution_to_state(solution_, S_next_);
  commit_state(h, S_next_);

  // update the timestep size
  if (dt_solver < dt_ && dt_solver >= h) {
    // We took a smaller step than we recommended, and it worked fine (not
    // suprisingly).  Likely this was due to constraints from other PKs or
    // vis.  Do not reduce our recommendation.
  } else {
    dt_ = dt_solver;
  }

  return false;
};


// -----------------------------------------------------------------------------
// Allows a PK to reject solutions, which forces a timestep cut.
// -----------------------------------------------------------------------------
  bool PKBDFBase::is_admissible(Teuchos::RCP<const TreeVector> up) { return true; }


// -----------------------------------------------------------------------------
// Allows a PK to modify the initial guess.
// -----------------------------------------------------------------------------
bool PKBDFBase::modify_predictor(double h, Teuchos::RCP<TreeVector> up) { return false; }


} // namespace
