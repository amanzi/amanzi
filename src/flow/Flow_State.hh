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
namespace AmanziFlow {

class Flow_State : public PK_State {

public:

  explicit Flow_State(Teuchos::RCP<AmanziMesh::Mesh> mesh);
  explicit Flow_State(Teuchos::RCP<State> S);
  explicit Flow_State(State& S);
  Flow_State(Flow_State& other, PKStateConstructMode mode);

  void Initialize();

  // access methods
  Teuchos::RCP<AmanziGeometry::Point> gravity();

  Teuchos::RCP<double> fluid_density() { return S_->GetScalarData("fluid_density", name_); }
  Teuchos::RCP<double> fluid_viscosity() {
    return S_->GetScalarData("fluid_viscosity", name_); }
  Teuchos::RCP<Epetra_Vector> pressure() {
    return S_->GetFieldData("pressure", name_)->ViewComponent("cell", ghosted_)(0); }
  Teuchos::RCP<Epetra_Vector> lambda() {
    return S_->GetFieldData("pressure", name_)->ViewComponent("face", ghosted_)(0); }
  Teuchos::RCP<Epetra_Vector> darcy_flux() {
    return S_->GetFieldData("darcy_flux", name_)->ViewComponent("face", ghosted_)(0); }

  Teuchos::RCP<Epetra_Vector> vertical_permeability() {
    return S_->GetFieldData("permeability", name_)->ViewComponent("cell", ghosted_)(0); }
  Teuchos::RCP<Epetra_Vector> horizontal_permeability() {
    return S_->GetFieldData("permeability", name_)->ViewComponent("cell", ghosted_)(1); }
  Teuchos::RCP<Epetra_Vector> porosity() {
    return S_->GetFieldData("porosity", name_)->ViewComponent("cell", ghosted_)(0); }
  Teuchos::RCP<Epetra_Vector> water_saturation() {
    return S_->GetFieldData("water_saturation", name_)->ViewComponent("cell", ghosted_)(0); }
  Teuchos::RCP<Epetra_Vector> prev_water_saturation() {
    return S_->GetFieldData("prev_water_saturation", name_)->ViewComponent("cell", ghosted_)(0); }

  Teuchos::RCP<Epetra_Vector> specific_storage() {
    return S_->GetFieldData("specific_storage", name_)->ViewComponent("cell", ghosted_)(0); }
  Teuchos::RCP<Epetra_Vector> specific_yield() {
    return S_->GetFieldData("specific_yield", name_)->ViewComponent("cell", ghosted_)(0); }

  double ref_fluid_density() { return *fluid_density(); }
  double ref_fluid_viscosity() { return *fluid_viscosity(); }
  Epetra_Vector& ref_pressure() { return *pressure(); }
  Epetra_Vector& ref_lambda() { return *lambda(); }
  Epetra_Vector& ref_darcy_flux() { return *darcy_flux(); }
  Epetra_MultiVector& ref_darcy_velocity() { return *darcy_velocity(); }
  const AmanziGeometry::Point& ref_gravity() { return *gravity(); }

  Epetra_Vector& ref_vertical_permeability() { return *vertical_permeability(); }
  Epetra_Vector& ref_horizontal_permeability() { return *horizontal_permeability(); }
  Epetra_Vector& ref_porosity() { return *porosity(); }
  Epetra_Vector& ref_water_saturation() { return *water_saturation(); }
  Epetra_Vector& ref_prev_water_saturation() { return *prev_water_saturation(); }

  Epetra_Vector& ref_specific_storage() { return *specific_storage(); }
  Epetra_Vector& ref_specific_yield() { return *specific_yield(); }

  // miscaleneous
  double get_time() { return (S_ == Teuchos::null) ? -1.0 : S_->time(); }

  // debug routines
  void set_fluid_density(double rho);
  void set_fluid_viscosity(double mu);
  void set_porosity(double phi);
  void set_pressure_hydrostatic(double z0, double p0);
  void set_permeability(double Kh, double Kv);
  void set_permeability(double Kh, double Kv, const string region);
  void set_gravity(double g);
  void set_specific_storage(double ss);

protected:
  void Construct_();


};
} // namespace

#endif
