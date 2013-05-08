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

  typedef AmanziGeometry::Point f_flux_t(const AmanziGeometry::Point&, const double);
  typedef double f_conc_t(const AmanziGeometry::Point&, const double);


 public:

  explicit Transport_State(Teuchos::RCP<AmanziMesh::Mesh> mesh, const int);
  explicit Transport_State(Teuchos::RCP<State> S);
  explicit Transport_State(State& S);
  Transport_State(Transport_State& other, PKStateConstructMode mode);

  void Initialize();

  // access methods for state variables
  // const methods
  Teuchos::RCP<const Epetra_MultiVector> total_component_concentration() const {
    return S_->GetFieldData("total_component_concentration", name_)->ViewComponent("cell", true); }

  Teuchos::RCP<const Epetra_Vector> porosity() const {
    return Teuchos::rcpFromRef(*(*S_->GetFieldData("porosity", name_)->ViewComponent("cell", ghosted_))(0)); }
  Teuchos::RCP<const Epetra_Vector> water_saturation() const {
    return Teuchos::rcpFromRef(*(*S_->GetFieldData("water_saturation", name_)->ViewComponent("cell", ghosted_))(0)); }
  Teuchos::RCP<const Epetra_Vector> prev_water_saturation() const {
    return Teuchos::rcpFromRef(*(*S_->GetFieldData("prev_water_saturation", name_)->ViewComponent("cell", ghosted_))(0)); }
  Teuchos::RCP<const Epetra_Vector> darcy_flux() const {
    return Teuchos::rcpFromRef(*(*S_->GetFieldData("darcy_flux", name_)->ViewComponent("face", ghosted_))(0)); }
  Teuchos::RCP<const double> water_density() const {
    return S_->GetScalarData("fluid_density", name_); }

  // non const access methods
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
  Teuchos::RCP<double> water_density() {
    return S_->GetScalarData("fluid_density", name_); }

  Epetra_MultiVector& ref_total_component_concentration() { return *total_component_concentration(); }
  Epetra_Vector& ref_porosity() { return *porosity(); }
  Epetra_Vector& ref_water_saturation() { return *water_saturation(); }
  Epetra_Vector& ref_prev_water_saturation() { return *prev_water_saturation(); }
  Epetra_Vector& ref_darcy_flux() { return *darcy_flux(); }
  double& ref_water_density() { return *water_density(); }


  // component name and number accessors
  int get_component_number(const std::string component_name);
  std::string get_component_name(const int component_number);

  // time
  double initial_time() const { return S_->time(); }
  double intermediate_time() const { return S_->intermediate_time(); }
  double final_time() const { return S_->final_time(); }

  // routines that may not belong in state
  void InterpolateCellVector(const Epetra_Vector& v0,
                             const Epetra_Vector& v1,
                             double dT_int, double dT,
                             Epetra_Vector& v_int);

  // debug routines
  void set_porosity(const double p = 0.2);
  void set_water_saturation(const double ws = 1.0);
  void set_prev_water_saturation(const double ws = 0.0);
  void set_water_density(const double wd = 1000.0);
  void set_total_component_concentration(f_conc_t f, const double t=0.0);
  void set_total_component_concentration(const double tcc);  
  void set_darcy_flux(f_flux_t f, const double t=0.0);
  void set_darcy_flux(const AmanziGeometry::Point& u);  
  void error_total_component_concentration(f_conc_t f, double t, double* L1, double* L2);

 protected:
  void Construct_();

 private:
  // not implemented
  Transport_State(const Transport_State& other);

 protected:
  std::map<std::string,int> comp_numbers_;
  std::vector<std::string> comp_names_;

};
} // namespace
} // namespace

#endif
