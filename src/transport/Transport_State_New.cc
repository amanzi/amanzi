
#include "Transport_State_New.hh"

namespace Amanzi {
namespace AmanziTransport {

Transport_State_New::Transport_State_New(Teuchos::RCP<AmanziMesh::Mesh> mesh) :
    PK_State(std::string("state"), mesh) {
  Construct_();
}

Transport_State_New::Transport_State_New(Teuchos::RCP<State> S) :
    PK_State(std::string("state"), S) {
  Construct_();
}

Transport_State_New::Transport_State_New(State& S) :
    PK_State(std::string("state"), S) {
  Construct_();
}

Transport_State_New::Transport_State_New(Transport_State_New& other,
        PKStateConstructMode mode) :
    PK_State(other),
    comp_numbers_(other.comp_numbers_),
    comp_names_(other.comp_names_) {

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


void Transport_State_New::Construct_() {
  // Somehow much set up comp_names and comp_numbers!

  // Require data, all owned
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

void Transport_State_New::Initialize() {
  S_->GetField("total_component_concentration",name_)->set_initialized();
  S_->GetField("fluid_density",name_)->set_initialized();
  S_->GetField("porosity",name_)->set_initialized();
  S_->GetField("water_saturation",name_)->set_initialized();
  S_->GetField("prev_water_saturation",name_)->set_initialized();
  S_->GetField("darcy_flux",name_)->set_initialized();
}


int Transport_State_New::get_component_number(const std::string component_name) {
  std::map<std::string,int>::const_iterator lb =
    comp_numbers_.lower_bound(component_name);
  if (lb != comp_numbers_.end() && !(comp_numbers_.key_comp()(key, lb->first))) {
    return lb->second;
  } else {
    return -1;
  }
}

std::string get_component_name(const int component_number) {
  return comp_names_[component_number];
}


} // namspace
} // namspace
