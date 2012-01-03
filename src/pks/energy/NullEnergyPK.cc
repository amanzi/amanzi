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
                           Teuchos::RCP<TreeVector>& soln) :
    energy_plist_(energy_plist) {

  // require fields for the state and solution
  S->require_field("temperature", FIELD_LOCATION_CELL, "energy");
  S->get_field_record("temperature")->set_io_vis(true);
  Teuchos::RCP<Epetra_MultiVector> soln_temp = Teuchos::rcp(new Epetra_MultiVector(*S->get_field("temperature")));
  soln->PushBack(soln_temp);

  S->require_field("temperature_dot", FIELD_LOCATION_CELL, "energy");

  T_ = energy_plist.get<double>("Constant temperature", 290.0);
};

// initialize ICs
void NullEnergyPK::initialize(Teuchos::RCP<State>& S,
                              Teuchos::RCP<TreeVector>& soln) {
  S->set_field("temperature", "energy", T_);
  S->get_field_record("temperature")->set_initialized();

  S->set_field("temperature_dot", "energy", 0.0);
  S->get_field_record("temperature_dot")->set_initialized();

  state_to_solution(*S, soln);
};

// Copy temperature field from state into solution vector.
// Used to reset timestep and initialize the solution vector.
  void NullEnergyPK::state_to_solution(const State& S,
                                     Teuchos::RCP<TreeVector>& soln) {
  *((*soln)[0]) = *S.get_field("temperature");
};

// Copy temperature fields from state into solution vector.
void NullEnergyPK::state_to_solution(const State& S,
                                     Teuchos::RCP<TreeVector>& soln,
                                     Teuchos::RCP<TreeVector>& soln_dot) {
  *(*soln)[0] = *S.get_field("temperature");
  *(*soln_dot)[0] = *S.get_field("temperature_dot");
};

// Copy temperature fields from solution vector into state.  Used within
// compute_f() of strong couplers to set the current iterate in the state for
// use with other PKs.
void NullEnergyPK::solution_to_state(const TreeVector& soln,
                                     Teuchos::RCP<State>& S) {
  S->set_field("temperature", "energy", *soln[0]);
};

void NullEnergyPK::solution_to_state(const TreeVector& soln,
                                     const TreeVector& soln_dot,
                                     Teuchos::RCP<State>& S) {
  S->set_field("temperature", "energy", *soln[0]);
  S->set_field("temperature_dot", "energy", *soln_dot[0]);
};

bool NullEnergyPK::advance(double dt, Teuchos::RCP<TreeVector> &soln) {
  *soln = T_;
  solution_to_state(*soln, S_next_);
  return false;
};

void NullEnergyPK::compute_f(const double t, const Vector& u,
                             const Vector& udot, Vector& f) {
  f = u;
  f.Shift(-T_);
};
} // namespace
