
#include "Transport_State.hh"

namespace Amanzi {
namespace AmanziTransport {

Transport_State::Transport_State(Teuchos::RCP<AmanziMesh::Mesh> mesh, const int ncomp) :
    PK_State(std::string("state"), mesh) {

  // create dummy component names
  for (int i=0; i<ncomp; ++i) {
    std::stringstream ss;
    ss << "Component " << i; 
    comp_names_.push_back(ss.str());
    comp_numbers_[ss.str()] = i;
  }

  Construct_();
}

Transport_State::Transport_State(Teuchos::RCP<State> S) :
    PK_State(std::string("state"), S) {
  Construct_();
}

Transport_State::Transport_State(State& S) :
    PK_State(std::string("state"), S) {
  Construct_();
}

Transport_State::Transport_State(Transport_State& other,
        PKStateConstructMode mode) :
    PK_State(other, STATE_CONSTRUCT_MODE_COPY_POINTERS),
    comp_numbers_(other.comp_numbers_),
    comp_names_(other.comp_names_) {

  if (mode == CONSTRUCT_MODE_VIEW_DATA) {
    ghosted_ = false;  // non-ghosted views
  } else if (mode == CONSTRUCT_MODE_VIEW_DATA_GHOSTED) {
    ghosted_ = true; // no guarantees -- if other is not ghosted, this is not
                     // ghosted either!
  } else if (mode == CONSTRUCT_MODE_COPY_DATA) {
    // Not currently needed by Transport?
    ASSERT(0);
  } else if (mode == CONSTRUCT_MODE_COPY_DATA_GHOSTED) {
    // Copy, with ghosted TCC and flux
    ghosted_ = true;

    CompositeVectorFactory fac_tcc;
    fac_tcc.SetMesh(mesh_);
    fac_tcc.SetComponent("cell", AmanziMesh::CELL, comp_names_.size());
    Teuchos::RCP<CompositeVector> tcc = fac_tcc.CreateVector(true);
    tcc->CreateData();
    *tcc->ViewComponent("cell",false) = *other.total_component_concentration();
    tcc->ScatterMasterToGhosted();
    S_->SetData("total_component_concentration", name_, tcc);

    CompositeVectorFactory fac;
    fac.SetMesh(mesh_);
    fac.SetComponent("face", AmanziMesh::FACE, 1);
    Teuchos::RCP<CompositeVector> flux = fac.CreateVector(true);
    flux->CreateData();
    *flux->ViewComponent("face",false) = *other.darcy_flux();
    flux->ScatterMasterToGhosted();
    S_->SetData("darcy_flux", name_, flux);
  }
}


void Transport_State::Construct_() {
  // Somehow much set up comp_names and comp_numbers!

  // Require data, all owned by "state" to cheat the system.
  S_->RequireScalar("fluid_density", name_);
  S_->RequireField("porosity", name_)->SetMesh(mesh_)->SetGhosted(false)
    ->SetComponent("cell", AmanziMesh::CELL, 1);
  S_->RequireField("water_saturation", name_)->SetMesh(mesh_)->SetGhosted(false)
    ->SetComponent("cell", AmanziMesh::CELL, 1);
  S_->RequireField("prev_water_saturation", name_)->SetMesh(mesh_)->SetGhosted(false)
    ->SetComponent("cell", AmanziMesh::CELL, 1);
  S_->RequireField("darcy_flux", name_)->SetMesh(mesh_)->SetGhosted(false)
    ->SetComponent("face", AmanziMesh::FACE, 1);

  S_->RequireField("total_component_concentration", name_)->SetMesh(mesh_)
    ->SetGhosted(false)
    ->SetComponent("cell", AmanziMesh::CELL, comp_names_.size());
}


void Transport_State::Initialize() {
  if (standalone_mode_) {
    S_->Setup();
    S_->GetField("total_component_concentration",name_)->set_initialized();
    S_->GetField("fluid_density",name_)->set_initialized();
    S_->GetField("porosity",name_)->set_initialized();
    S_->GetField("water_saturation",name_)->set_initialized();
    S_->GetField("prev_water_saturation",name_)->set_initialized();
    S_->GetField("darcy_flux",name_)->set_initialized();
    S_->Initialize();
  } else {
    // BEGIN REMOVE ME once flow tests pass --etc
    S_->GetFieldData("total_component_concentration", name_)->PutScalar(-1.);
    S_->GetFieldData("porosity", name_)->PutScalar(-1.);
    S_->GetFieldData("water_saturation", name_)->PutScalar(-1.);
    S_->GetFieldData("prev_water_saturation", name_)->PutScalar(-1.);
    S_->GetFieldData("darcy_flux", name_)->PutScalar(-1.);
    // END REMOVE ME
  }
}

int Transport_State::get_component_number(const std::string component_name) {
  std::map<std::string,int>::const_iterator lb =
    comp_numbers_.lower_bound(component_name);
  if (lb != comp_numbers_.end() && !(comp_numbers_.key_comp()(component_name, lb->first))) {
    return lb->second;
  } else {
    return -1;
  }
}

std::string Transport_State::get_component_name(const int component_number) {
  return comp_names_[component_number];
}


/* *******************************************************************
* Interpolate linearly in time between two values v0 and v1. The time 
* is measuared relative to value v0; so that v1 is at time dT. The
* interpolated data are at time dT_int.            
******************************************************************* */
void Transport_State::InterpolateCellVector(const Epetra_Vector& v0,
        const Epetra_Vector& v1, double dT_int, double dT,
        Epetra_Vector& v_int) {
  double a = dT_int / dT;
  double b = 1.0 - a;
  v_int.Update(b, v0, a, v1, 0.);
}


} // namspace
} // namspace
