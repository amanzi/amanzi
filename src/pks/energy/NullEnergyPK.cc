/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

#include "NullEnergyPK.hh"

//namespace Amanzi {

NullEnergyPK::NullEnergyPK(Teuchos::ParameterList&, Teuchos::RCP<State> S) :
    energy_plist_(plist), S_(S) {
  // require fields for the state
  S->require_field("temperature", Amanzi::AmanziMesh::CELL, "energy");
  S->get_field_record("temperature")->set_io_vis(true);

  S->require_field("dtemperature", Amanzi::AmanziMesh::CELL, "energy");

  T_ = plist.get<double>("Constant temperature", 290.0);
};

NullEnergyPK::initialize() {
  S->set_field("temperature", "energy", T_);
  S->set_field("temperature", "energy", T_);
};

bool NullEnergyPK::advance(double dt, const Teuchos::RCP<State> &S0,
          Teuchos::RCP<State> &S1, Teuchos::RCP<TreeVector> &solution) {
  *solution = T_;
  
  return false;
};

void NullEnergyPK::compute_f(const double t, const Vector& u,
                             const Vector& udot, Vector& f) {
  f = u;
  Vector work(u);
  work = 1.0;
  f.Update(-T_, work, 1.0);
};
    
void NullEnergyPK::solution_to_state(const TreeVector& u,
                                     const TreeVector& udot) {
  solution_to_state(u, S_);
  Teuchos::RCP<const Epetra_MultiVector> dtempdata = udot[0];
  S_->set_field("dtemperature", "energy", dtempdata);
};

void NullEnergyPK::solution_to_state(const TreeVector& u,
                                     Teuchos::RCP<State> &S) {
  Teuchos::RCP<const Epetra_MultiVector> tempdata = u[0];
  S->set_field("temperature", "energy", tempdata);
};
