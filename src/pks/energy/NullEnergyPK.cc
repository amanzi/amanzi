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

#include "NullEnergyPK.hh"

namespace Amanzi {

NullEnergyPK::NullEnergyPK(Teuchos::ParameterList& energy_plist,
                           Teuchos::RCP<State>& S,
                           Teuchos::RCP<TreeVector>& solution) :
    energy_plist_(energy_plist) {

  solution_ = solution;

  // require fields for the state and solution
  S->require_field("temperature", FIELD_LOCATION_CELL, "energy");
  S->get_field_record("temperature")->set_io_vis(true);
  Teuchos::RCP<Epetra_MultiVector> temp_ptr = S->get_field("temperature", "energy");
  solution_->PushBack(temp_ptr);

  S->require_field("temperature_dot", FIELD_LOCATION_CELL, "energy");

  T_ = energy_plist.get<double>("Constant temperature", 290.0);
};

// initialize ICs
void NullEnergyPK::initialize(Teuchos::RCP<State>& S) {
  S->set_field("temperature", "energy", T_);
  S->get_field_record("temperature")->set_initialized();

  S->set_field("temperature_dot", "energy", 0.0);
  S->get_field_record("temperature_dot")->set_initialized();

  // model evaluator params
  atol_ = energy_plist_.get<double>("Absolute error tolerance",1.0);
  rtol_ = energy_plist_.get<double>("Relative error tolerance",1e-5);
};


// Pointer copy of state to solution
void NullEnergyPK::state_to_solution(State& S,
                                     Teuchos::RCP<TreeVector>& soln) {
  ((*soln)[0]) = S.get_field("temperature", "energy");
};

// Pointer copy temperature fields from state into solution vector.
void NullEnergyPK::state_to_solution(State& S, Teuchos::RCP<TreeVector>& soln,
                                     Teuchos::RCP<TreeVector>& soln_dot) {
  (*soln)[0] = S.get_field("temperature", "energy");
  (*soln_dot)[0] = S.get_field("temperature_dot", "energy");
};

// Pointer copy temperature fields from solution vector into state.  Used within
// compute_f() of strong couplers to set the current iterate in the state for
// use with other PKs.
void NullEnergyPK::solution_to_state(TreeVector& soln, Teuchos::RCP<State>& S) {
  Teuchos::RCP<Epetra_MultiVector> temp_ptr = soln[0];
  S->set_field_pointer("temperature", "energy", temp_ptr);
};

void NullEnergyPK::solution_to_state(TreeVector& soln, TreeVector& soln_dot,
                                     Teuchos::RCP<State>& S) {
  Teuchos::RCP<Epetra_MultiVector> temp_ptr = soln[0];
  S->set_field_pointer("temperature", "energy", temp_ptr);
  Teuchos::RCP<Epetra_MultiVector> temp_dot_ptr = soln_dot[0];
  S->set_field_pointer("temperature_dot", "energy", temp_dot_ptr);
};

// Advance methods calculate the constant value
// -- advance using the analytic value
bool NullEnergyPK::advance_analytic(double dt) {
  *solution_ = T_;
  return false;
};

// -- advance using the BDF integrator
bool NullEnergyPK::advance_bdf(double dt) {
};

// -- call your favorite
bool NullEnergyPK::advance(double dt) {
  return advance_analytic(dt);
};

// overwrite the state pointers, and make sure the solutions T pointer points to S_next's T
void NullEnergyPK::set_states(Teuchos::RCP<const State>& S, Teuchos::RCP<State>& S_next) {
  S_ = S;
  S_next_ = S_next;
  // pointer copy the state's temperature field to be the solution's temperature field
  (*solution_)[0] = S_next->get_field("temperature", "energy");
};

// Methods for the BDF integrator
// -- residual
void NullEnergyPK::fun(const double t, const TreeVector& soln, const TreeVector& udot,
                       TreeVector& f) {
  f = soln;
  f.Shift(-T_);
};

// -- preconditioning (currently none)
void NullEnergyPK::precon(const TreeVector& u, TreeVector& Pu) {
  Pu = u;
};

// computes a norm on u-du and returns the result
double NullEnergyPK::enorm(const TreeVector& u, const TreeVector& du) {
  double enorm_val = 0.0;
  Epetra_Vector temp_vec = *(*u[0])(0);
  Epetra_Vector temp_dot_vec = *(*du[0])(0);
  for (unsigned int lcv=0; lcv != u[0]->MyLength(); ++lcv) {
    double tmp = abs(temp_dot_vec[lcv])/(atol_ + rtol_*abs(temp_vec[lcv]));
    enorm_val = std::max<double>(enorm_val, tmp);
  }

#ifdef HAVE_MPI
  double buf = enorm_val;
  MPI_Allreduce(&buf, &enorm_val, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#endif

  return enorm_val;
};

} // namespace
