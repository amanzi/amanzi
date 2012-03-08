/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
   ATS

   License: see $ATS_DIR/COPYRIGHT
   Author: Ethan Coon

   Implementation for the NullEnergy PK.  This PK simply provides a constant
   temperature, and is provided for testing with other PKs that depend upon an
   energy equation.  This could easily be provided by the state as an independent
   variable, but this is nice for testing the full hierarchy with a simple PK.

   Example usage:

   <ParameterList name="energy">
   <Parameter name="PK model" type="string" value="Constant Temperature"/>
   <Parameter name="Constant Temperature" type="double" value="290.0"/>
   </ParameterList>

   ------------------------------------------------------------------------- */

#include "composite_vector.hh"
#include "bdf2_time_integrator.hh"
#include "constant_temperature.hh"

namespace Amanzi {
namespace Energy {

ConstantTemperature::ConstantTemperature(Teuchos::ParameterList& energy_plist,
        Teuchos::RCP<State>& S, Teuchos::RCP<TreeVector>& solution) :
  energy_plist_(energy_plist) {

  solution_ = solution;

  // require fields for the state and solution
  S->RequireField("temperature", "energy", AmanziMesh::CELL);
  S->GetRecord("temperature","energy")->set_io_vis(true);
  Teuchos::RCP<CompositeVector> temp = S->GetFieldData("temperature", "energy");
  solution_->set_data(temp);

  S->RequireField("temperature_dot", "energy", AmanziMesh::CELL);

  // check if we need to make a time integrator
  if (!energy_plist_.get<bool>("Strongly Coupled PK", false)) {
    Teuchos::RCP<Teuchos::ParameterList> bdf2_plist_p =
      Teuchos::rcp(new Teuchos::ParameterList(energy_plist_.sublist("Time integrator")));
    time_stepper_ = Teuchos::rcp(new BDF2TimeIntegrator(this, bdf2_plist_p, solution_));
  }
};

// initialize ICs
void ConstantTemperature::initialize(const Teuchos::RCP<State>& S) {
  // constant initial temperature
  T_ = energy_plist_.get<double>("Constant temperature", 290.0);
  S->GetFieldData("temperature", "energy")->PutScalar(T_);
  S->GetRecord("temperature", "energy")->set_initialized();

  S->GetFieldData("temperature_dot", "energy")->PutScalar(0.0);
  S->GetRecord("temperature_dot", "energy")->set_initialized();

  if (!energy_plist_.get<bool>("Strongly Coupled PK", false)) {
    // model evaluator params
    // -- tolerances
    atol_ = energy_plist_.get<double>("Absolute error tolerance",1.0);
    rtol_ = energy_plist_.get<double>("Relative error tolerance",1e-5);

    // -- initialize time derivative
    Teuchos::RCP<TreeVector> solution_dot = Teuchos::rcp( new TreeVector(*solution_));
    solution_dot->PutScalar(0.0);

    // -- set initial state
    time_stepper_->set_initial_state(S->time(), solution_, solution_dot);
  }
};


// Pointer copy of state to solution
void ConstantTemperature::state_to_solution(const Teuchos::RCP<State>& S,
        const Teuchos::RCP<TreeVector>& soln) {
  soln->set_data(S->GetFieldData("temperature", "energy"));
};

// Pointer copy temperature fields from solution vector into state.  Used within
// compute_f() of strong couplers to set the current iterate in the state for
// use with other PKs.
void ConstantTemperature::solution_to_state(const Teuchos::RCP<TreeVector>& soln,
        const Teuchos::RCP<State>& S) {
  S->SetData("temperature", "energy", soln->data());
};

// Advance methods calculate the constant value
// -- advance using the analytic value
bool ConstantTemperature::advance_analytic_(double dt) {
  solution_->PutScalar(T_);
  return false;
};

// -- advance using the BDF integrator
bool ConstantTemperature::advance_bdf_(double dt) {
  state_to_solution(S_next_, solution_);

  // take the bdf timestep
  double h = dt;
  time_stepper_->time_step(h, solution_);
  time_stepper_->commit_solution(h, solution_);

  // In the case where this is a leaf, and therefore advancing itself (and not
  // within a strongly coupled solver), we will call the local residual
  // function only.  This local residual function need NOT copy the guess for
  // u/u_dot into the state (as it is a leaf and therefore already has access.
  // Therefore, the S_next's temperature pointer was not overwritten, and we
  // need not copy it back into S_next, like is required in StrongMPC.
  return false;
};

// -- call your favorite
bool ConstantTemperature::advance(double dt) {
  return advance_analytic_(dt);
};

// Methods for the BDF integrator
// -- residual
void ConstantTemperature::fun(double t_old, double t_new, Teuchos::RCP<TreeVector> u_old,
        Teuchos::RCP<TreeVector> u_new, Teuchos::RCP<TreeVector> f) {
  *f = *u_new;
  f->Shift(-T_);
};

// -- preconditioning (currently none)
void ConstantTemperature::precon(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu) {
  *Pu = *u;
};

// computes a norm on u-du and returns the result
double ConstantTemperature::enorm(Teuchos::RCP<const TreeVector> u,
        Teuchos::RCP<const TreeVector> du) {
  double enorm_val = 0.0;
  Teuchos::RCP<const Epetra_MultiVector> temp_vec = u->data()->ViewComponent("temperature", false);
  Teuchos::RCP<const Epetra_MultiVector> temp_dot_vec = du->data()->ViewComponent("temperature_dot", false);

  for (unsigned int lcv=0; lcv != temp_vec->MyLength(); ++lcv) {
    double tmp = abs((*(*temp_dot_vec)(0))[lcv])/(atol_ + rtol_*abs((*(*temp_vec)(0))[lcv]));
    enorm_val = std::max<double>(enorm_val, tmp);
  }

#ifdef HAVE_MPI
  double buf = enorm_val;
  MPI_Allreduce(&buf, &enorm_val, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#endif

  return enorm_val;
};

} // namespace
} // namespace
