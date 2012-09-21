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
#include "composite_vector_factory.hh"
#include "bdf2_time_integrator.hh"
#include "constant_temperature.hh"

namespace Amanzi {
namespace Energy {

RegisteredPKFactory<ConstantTemperature> ConstantTemperature::reg_("constant temperature energy");

ConstantTemperature::ConstantTemperature(Teuchos::ParameterList& energy_plist,
        const Teuchos::RCP<State>& S, const Teuchos::RCP<TreeVector>& solution) :
  energy_plist_(energy_plist) {

  solution_ = solution;

  // require fields for the state and solution
  Teuchos::RCP<CompositeVectorFactory> factory =
    S->RequireField("temperature", "energy");

  // Set up the data structure.
  // Since there is only one field, we'll do this manually.
  factory->SetMesh(S->GetMesh());
  factory->SetComponent("cell", AmanziMesh::CELL, 1);
  factory->SetGhosted(true);

  // Note that this the above lines are equivalent to the fancier/more concise
  // version:
  //  S->RequireField("temperature", "energy")->SetMesh(S->GetMesh())->
  //            SetComponent("cell", AmanziMesh::CELL, 1)->SetGhosted(true);

  S->GetField("temperature","energy")->set_io_vis(true);
  Teuchos::RCP<CompositeVector> temp = S->GetFieldData("temperature", "energy");
};

// initialize ICs
void ConstantTemperature::initialize(const Teuchos::RCP<State>& S) {
  // This pk provides a constant temperature, given by the intial temp.
  // Therefore we store the initial temp to evaluate changes.
  Teuchos::RCP<CompositeVector> temp = S->GetFieldData("temperature", "energy");
  temp0_ = Teuchos::rcp(new CompositeVector(*temp));
  *temp0_ = *temp;

  // initialize the time integrator and nonlinear solver
  solution_->set_data(temp);
  atol_ = energy_plist_.get<double>("Absolute error tolerance",1.0);
  rtol_ = energy_plist_.get<double>("Relative error tolerance",1.0);

  // check if we need to make a time integrator
  if (!energy_plist_.get<bool>("Strongly Coupled PK", false)) {
    // -- instantiate the time integrator
    Teuchos::RCP<Teuchos::ParameterList> bdf2_plist_p =
      Teuchos::rcp(new Teuchos::ParameterList(energy_plist_.sublist("Time integrator")));
    time_stepper_ = Teuchos::rcp(new BDF2TimeIntegrator(this, bdf2_plist_p, solution_));

    // -- initialize time derivative
    Teuchos::RCP<TreeVector> solution_dot = Teuchos::rcp(new TreeVector(*solution_));
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
  *solution_->data() = *temp0_;
  return false;
};

// -- advance using the BDF integrator
bool ConstantTemperature::advance_bdf_(double dt) {
  state_to_solution(S_next_, solution_);

  // take the bdf timestep
  double h = dt;
  time_stepper_->time_step(h, solution_);

  // commit the step as successful
  time_stepper_->commit_solution(h, solution_);
  solution_to_state(solution_, S_next_);
  commit_state(h, S_next_);

  return false;
};

// -- call your favorite
bool ConstantTemperature::advance(double dt) {
  return advance_bdf_(dt);
};

// Methods for the BDF integrator
// -- residual
void ConstantTemperature::fun(double t_old, double t_new, Teuchos::RCP<TreeVector> u_old,
        Teuchos::RCP<TreeVector> u_new, Teuchos::RCP<TreeVector> f) {
  *f = *u_new;
  f->data()->Update(-1.0, *temp0_, 1.0); // T - T0
};

// -- preconditioning (the identity matrix)
void ConstantTemperature::precon(Teuchos::RCP<const TreeVector> u,
        Teuchos::RCP<TreeVector> Pu) {
  *Pu = *u;
};


// -- update the preconditioner (no need to do anything)
void ConstantTemperature::update_precon(double t,
        Teuchos::RCP<const TreeVector> up, double h) {};


// -- computes a norm on (u, du) and returns the result
double ConstantTemperature::enorm(Teuchos::RCP<const TreeVector> u,
        Teuchos::RCP<const TreeVector> du) {
  double enorm_val = 0.0;
  Teuchos::RCP<const Epetra_MultiVector> temp_vec = u->data()->ViewComponent("cell", false);
  Teuchos::RCP<const Epetra_MultiVector> temp_dot_vec = du->data()->ViewComponent("cell", false);

  for (int lcv=0; lcv != temp_vec->MyLength(); ++lcv) {
    double tmp = abs((*(*temp_dot_vec)(0))[lcv]) /
      (atol_ + rtol_*abs((*(*temp_vec)(0))[lcv]));
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
