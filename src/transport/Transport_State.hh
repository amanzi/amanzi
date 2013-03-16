/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
Amanzi Flow

License: see COPYRIGHT
Author: Ethan Coon

Interface layer between Flow and State, this is a harness for
accessing the new state-dev from the old Flow PK.

 ------------------------------------------------------------------------- */

#ifndef AMANZI_FLOW_STATE_NEW_HH_
#define AMANZI_FLOW_STATE_NEW_HH_

#include "PK_State.hh"

namespace Amanzi {
namespace AmanziTransport {

class Transport_State : public PK_State {

  explicit Transport_State(Teuchos::RCP<AmanziMesh::Mesh> mesh);
  explicit Transport_State(Teuchos::RCP<State> S);
  explicit Transport_State(State& S);
  Transport_State(Transport_State& other, PKStateConstructMode mode);

  void Initialize();

  // access methods for state variables
  Teuchos::RCP<Epetra_MultiVector> total_component_concentration() {
    return S_->GetFieldData("total_component_concentration", name_)->ViewComponent("cell", true); }

  Teuchos::RCP<Epetra_Vector> porosity() {
    return Teuchos::rcpFromRef(*(*S_->GetFieldData("porosity", name_)->ViewComponent("cell", ghosted_))(0)); }
  Teuchos::RCP<Epetra_Vector> water_saturation() {
    return Teuchos::rcpFromRef(*(*S_->GetFieldData("water_saturation", name_)->ViewComponent("cell", ghosted_))(0)); }
  Teuchos::RCP<Epetra_Vector> prev_water_saturation() {
    return Teuchos::rcpFromRef(*(*S_->GetFieldData("prev_water_saturation", name_)->ViewComponent("cell", ghosted_))(0)); }
  Teuchos::RCP<Epetra_Vector> darcy_flux() {
    return Teuchos::rcpFromRef(*(*S_->GetFieldData("darcy_flux", name_)->ViewComponent("face", ghosted_))(0)); }
  Teuchos::RCP<Epetra_Vector> water_density() {
    return S_->GetScalarData("fluid_density", name_); }

  Epetra_MultiVector& ref_total_component_concentration() { return *total_component_concentration(); }
  Epetra_Vector& ref_porosity() { return *porosity(); }
  Epetra_Vector& ref_water_saturation() { return *water_saturation(); }
  Epetra_Vector& ref_prev_water_saturation() { return *prev_water_saturation(); }
  Epetra_Vector& ref_darcy_flux() { return *darcy_flux(); }
  Epetra_Vector& ref_water_density() { return *water_density(); }


  // component name and number accessors
  int get_component_number(const std::string component_name);
  std::string get_component_name(const int component_number);

protected:
  void Construct_();

protected:
  std::map<std::string,int> comp_numbers_;
  std::vector<std::string> comp_names_;

};
} // namespace

#endif
