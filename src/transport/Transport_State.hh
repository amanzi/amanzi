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

class Transport_State_New : public PK_State {

  explicit Transport_State_New(Teuchos::RCP<AmanziMesh::Mesh> mesh) :
    PK_State(std::string("state"), mesh) { Construct_(); }

  explicit Transport_State_New(Teuchos::RCP<State> S) :
    PK_State(std::string("state"), S) { Construct_(); }

  explicit Transport_State_New(State& S) :
    PK_State(std::string("state"), S) {}

  Transport_State_New(Transport_State_New& other, PKStateConstructMode mode) :
      PK_State(other) {
    if (mode == PK_STATE_CONSTRUCT_MODE_VIEW_DATA) {
      ghosted_ = false;

    } else if (mode == PK_STATE_CONSTRUCT_MODE_VIEW_DATA_GHOSTED) {
      ASSERT(0);
    } else if (mode == PK_STATE_CONSTRUCT_MODE_COPY_DATA) {
      ASSERT(0);
    } else if (mode == PK_STATE_CONSTRUCT_MODE_COPY_DATA_GHOSTED) {
      ghosted_ = true;

      CompositeVectorFactory fac_tcc();
      fac_tcc.SetMesh(mesh_);
      fac_tcc.SetComponent("cell", AmanziMesh::CELL, comp_names_.size());
      Teuchos::RCP<CompositeVector> tcc = fac_tcc.CreateVector(true);
      *tcc->ViewComponent("cell",false) = *other->total_component_concentration();
      tcc->ScatterMasterToGhosted();
      S_->SetData("total_component_concentration", name_, flux);

      CompositeVectorFactory fac();
      fac.SetMesh(mesh_);
      fac.SetComponent("face", AmanziMesh::FACE, 1);
      Teuchos::RCP<CompositeVector> flux = fac.CreateVector(true);
      flux->PutScalar(0.);
      *flux->ViewComponent("face",false) = *other.darcy_flux();
      flux->ScatterMasterToGhosted();
      S_->SetData("darcy_flux", name_, flux);
    }
  }

  void Initialize();


  // access methods for state variables
  Teuchos::RCP<Epetra_MultiVector> total_component_concentration() {
    return S_->GetFieldData("total_component_concentration", name_)->ViewComponent("cell", true); }

  Teuchos::RCP<Epetra_Vector> porosity() {
    return S_->GetFieldData("porosity", name_)->ViewComponent("cell", true)(0); }
  Teuchos::RCP<Epetra_Vector> water_saturation() {
    return S_->GetFieldData("water_saturation", name_)->ViewComponent("cell", true)(0); }

  Teuchos::RCP<Epetra_Vector> prev_water_saturation() {
    return S_->GetFieldData("prev_water_saturation", name_)->ViewComponent("cell", true)(0); }

  Teuchos::RCP<Epetra_Vector> darcy_flux() {
    return S_->GetFieldData("darcy_flux", name_)->ViewComponent("face", true)(0); }

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
