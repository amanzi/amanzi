/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

#include "NullEnergyPK.hh"

namespace Amanzi {

NullEnergyPK::NullEnergyPK(Teuchos::ParameterList& energy_plist,
                           Teuchos::RCP<State>& S,
                           Teuchos::RCP<TreeVector>& soln) :
    energy_plist_(energy_plist), S_(S) {

  // require fields for the state and solution
  S->require_field("temperature", Amanzi::AmanziMesh::CELL, "energy");
  S->get_field_record("temperature")->set_io_vis(true);
  Teuchos::RCP<Epetra_MultiVector> soln_temp = Teuchos::rcp(new Epetra_MultiVector(*S->get_field("temperature", "energy")));
  soln->PushBack(soln_temp);

  S->require_field("dtemperature", Amanzi::AmanziMesh::CELL, "energy");

  T_ = energy_plist.get<double>("Constant temperature", 290.0);
};

void NullEnergyPK::initialize(Teuchos::RCP<State>& S,
                              Teuchos::RCP<TreeVector>& soln) {
  S_->set_field("temperature", "energy", T_);
  S_->set_field("temperature", "energy", T_);
  S_->set_field("dtemperature", "energy", 0.0);
  state_to_solution(S, soln);
};

// transfer operators
void NullEnergyPK::state_to_solution(Teuchos::RCP<const State>& S,
                                     Teuchos::RCP<TreeVector>& soln) {
  *(*soln[0]) = S->get_field("temperature");
};

void NullEnergyPK::state_to_solution(Teuchos::RCP<const State>& S,
                                     Teuchos::RCP<TreeVector>& soln,
                                     Teuchos::RCP<TreeVector>& soln_dot) {
  *(*soln)[0] = S->get_field("temperature");
  *(*soln_dot)[0] = S->get_field("dtemperature");
};

void NullEnergyPK::solution_to_state(Teuchos::RCP<const TreeVector>& soln,
                                     Teuchos::RCP<State>& S) {
  S->set_field("temperature", "energy", (*soln)[0]);
};

void NullEnergyPK::solution_to_state(Teuchos::RCP<const TreeVector>& soln,
                                     Teuchos::RCP<const TreeVector>& soln_dot,
                                     Teuchos::RCP<State>& S) {
  S->set_field("temperature", "energy", (*soln)[0]);
  S->set_field("dtemperature", "energy", (*soln_dot)[0]);
};

bool NullEnergyPK::advance(double dt, Teuchos::RCP<State> &S0,
          Teuchos::RCP<State> &S1, Teuchos::RCP<TreeVector> &soln) {
  *soln = T_;
  solution_to_state(soln, S1);
  return false;
};

void NullEnergyPK::compute_f(const double t, const Vector& u,
                             const Vector& udot, Vector& f) {
  f = u;
  f.Shift(-T_);
};
} // namespace
