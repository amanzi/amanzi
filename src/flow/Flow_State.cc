/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
Amanzi Flow

License: see COPYRIGHT
Author: Ethan Coon

Interface layer between Flow and State, this is a harness for
accessing the new state-dev from the old Flow PK.

------------------------------------------------------------------------- */

#include "Flow_State.hh"
#include "state_defs.hh"

namespace Amanzi {
namespace AmanziFlow {

/* ****************************************************************
* Constructors.
**************************************************************** */
Flow_State::Flow_State(Teuchos::RCP<AmanziMesh::Mesh> mesh) :
    PK_State(std::string("state"), mesh) {
  Construct_();
}

Flow_State::Flow_State(Teuchos::RCP<State> S) :
    PK_State(std::string("state"), S) {
  Construct_();
}

Flow_State::Flow_State(State& S) :
    PK_State(std::string("state"), S) {
  Construct_();
}


/* ****************************************************************
* Constructors.
**************************************************************** */
Flow_State::Flow_State(const Flow_State& other,
    PKStateConstructMode mode) :
    PK_State(other, STATE_CONSTRUCT_MODE_COPY_POINTERS) {
  if (mode == CONSTRUCT_MODE_VIEW_DATA) {
    ghosted_ = false; // non-ghosted views
  } else if (mode == CONSTRUCT_MODE_VIEW_DATA_GHOSTED) {
    ghosted_ = true; // no guarantees -- if other is not ghosted, this is not
                     // ghosted either!
  } else if (mode == CONSTRUCT_MODE_COPY_DATA) {
    // Not currently needed by Flow?
    ASSERT(0);
  } else if (mode == CONSTRUCT_MODE_COPY_DATA_GHOSTED) {
    // Copy data, making new vector for Darcy flux with ghosted entries.
    ghosted_ = true;

    CompositeVectorSpace fac;
    fac.SetMesh(mesh_);
    fac.SetComponent("face", AmanziMesh::FACE, 1);
    Teuchos::RCP<CompositeVector> flux = Teuchos::rcp(new CompositeVector(fac));
    //flux->CreateData();
    *flux->ViewComponent("face",false) = *other.darcy_flux();
    flux->ScatterMasterToGhosted();
    S_->SetData("darcy_flux", name_, flux);
  }
}


/* ****************************************************************
* TBW.
**************************************************************** */
void Flow_State::Construct_() {
  // for creating fields
  std::vector<AmanziMesh::Entity_kind> locations(2);
  std::vector<std::string> names(2);
  std::vector<int> ndofs(2,1);
  locations[0] = AmanziMesh::CELL; locations[1] = AmanziMesh::FACE;
  names[0] = "cell"; names[1] = "face";

  // Require data, all owned
  S_->RequireScalar("fluid_density", name_);
  S_->RequireScalar("fluid_viscosity", name_);
  S_->RequireConstantVector("gravity", name_, mesh_->space_dimension());
  
  if (!S_->HasField("pressure")) {
    S_->RequireField("pressure", name_)->SetMesh(mesh_)->SetGhosted(false)
      ->SetComponents(names, locations, ndofs);
  }
  if (!S_->HasField("permeability")) {
    S_->RequireField("permeability", name_)->SetMesh(mesh_)->SetGhosted(false)
      ->SetComponent("cell", AmanziMesh::CELL, mesh_->space_dimension());
  }
  if (!S_->HasField("porosity")) {
    S_->RequireField("porosity", name_)->SetMesh(mesh_)->SetGhosted(false)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
  }
  if (!S_->HasField("water_saturation")) {
    S_->RequireField("water_saturation", name_)->SetMesh(mesh_)->SetGhosted(false)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
  }
  if (!S_->HasField("prev_water_saturation")) {
    S_->RequireField("prev_water_saturation", name_)->SetMesh(mesh_)->SetGhosted(false)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
  }
  if (!S_->HasField("specific_storage")) {
    S_->RequireField("specific_storage", name_)->SetMesh(mesh_)->SetGhosted(false)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
  }
  if (!S_->HasField("specific_yield")) {
    S_->RequireField("specific_yield", name_)->SetMesh(mesh_)->SetGhosted(false)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
  }
  if (!S_->HasField("particle_density")) {
    S_->RequireField("particle_density", name_)->SetMesh(mesh_)->SetGhosted(false)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
  }
  if (!S_->HasField("darcy_flux")) {
    S_->RequireField("darcy_flux", name_)->SetMesh(mesh_)->SetGhosted(false)
      ->SetComponent("face", AmanziMesh::FACE, 1);
  }
  if (!S_->HasField("darcy_velocity")) {
    S_->RequireField("darcy_velocity", name_)->SetMesh(mesh_)->SetGhosted(false)
      ->SetComponent("cell", AmanziMesh::CELL, mesh_->space_dimension());
  }
  if (!S_->HasField("hydraulic_head")) {
    S_->RequireField("hydraulic_head", name_)->SetMesh(mesh_)->SetGhosted(false)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
  }
};


/* ****************************************************************
* TBW.
**************************************************************** */
void Flow_State::Initialize() {
  if (standalone_mode_) {
    S_->Setup();
    S_->GetField("fluid_density",name_)->set_initialized();
    S_->GetField("fluid_viscosity",name_)->set_initialized();
    S_->GetField("gravity",name_)->set_initialized();
    S_->GetField("pressure",name_)->set_initialized();
    S_->GetField("permeability",name_)->set_initialized();
    S_->GetField("porosity",name_)->set_initialized();
    S_->GetField("water_saturation",name_)->set_initialized();
    S_->GetField("prev_water_saturation",name_)->set_initialized();
    S_->GetField("specific_storage",name_)->set_initialized();
    S_->GetField("specific_yield",name_)->set_initialized();
    S_->GetField("darcy_flux",name_)->set_initialized();
    S_->GetField("darcy_velocity",name_)->set_initialized();
    S_->GetField("hydraulic_head",name_)->set_initialized();
    S_->GetField("particle_density",name_)->set_initialized();
    S_->InitializeFields();
  } else {
    
    // fields that might be initialized through the input
    // file need to be tested and initialized 'by hand' here

    if (!S_->GetField("pressure",name_)->initialized()) {
      pressure()->PutScalar(0.0);
      S_->GetField("pressure",name_)->set_initialized();
    }
    if (!S_->GetField("darcy_flux",name_)->initialized()) {
      darcy_flux()->PutScalar(0.0);
      S_->GetField("darcy_flux",name_)->set_initialized();
    }
    if (!S_->GetField("water_saturation",name_)->initialized()) {
      water_saturation()->PutScalar(1.0);
      S_->GetField("water_saturation",name_)->set_initialized();
    }
    if (!S_->GetField("prev_water_saturation",name_)->initialized()) {
      prev_water_saturation()->PutScalar(1.0);
      S_->GetField("prev_water_saturation",name_)->set_initialized();
    }
    if (!S_->GetField("specific_yield",name_)->initialized()) {
      specific_yield()->PutScalar(0.0);
      S_->GetField("specific_yield",name_)->set_initialized();
    }
    if (!S_->GetField("specific_storage",name_)->initialized()) {
      specific_storage()->PutScalar(0.0);
      S_->GetField("specific_storage",name_)->set_initialized();
    }
    if (!S_->GetField("particle_density",name_)->initialized()) {
      particle_density()->PutScalar(1.0);
      S_->GetField("particle_density",name_)->set_initialized();
    }    


    // aux data that is only for vis/checkpoint/observation
    // these will never be initialized in the input file, so
    // we need to set them to initialized here
    S_->GetField("darcy_velocity",name_)->set_initialized();
    S_->GetField("hydraulic_head",name_)->set_initialized();
  }
}


/* ****************************************************************
* TBW.
**************************************************************** */
Teuchos::RCP<AmanziGeometry::Point> Flow_State::gravity() {
  Teuchos::RCP<Epetra_Vector> gvec = S_->GetConstantVectorData("gravity", name_);
  Teuchos::RCP<AmanziGeometry::Point> gpoint =
      Teuchos::rcp(new AmanziGeometry::Point(gvec->MyLength()));
  for (int i = 0; i != gvec->MyLength(); ++i) (*gpoint)[i] = (*gvec)[i];
  return gpoint;
}


/* ****************************************************************
* Debug routines
**************************************************************** */
void Flow_State::set_fluid_density(double rho) {
  *fluid_density() = rho;
}

void Flow_State::set_fluid_viscosity(double mu) {
  *fluid_viscosity() = mu;
}

void Flow_State::set_porosity(double phi) {
  porosity()->PutScalar(phi);
}

void Flow_State::set_pressure_hydrostatic(double z0, double p0) {
  int dim = mesh_->space_dimension();
  int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::USED);

  double rho = *fluid_density();
  double g = (*gravity())[dim - 1];

  Epetra_Vector& pres = ref_pressure();
  for (int c=0; c!=ncells; ++c) {
    const AmanziGeometry::Point& xc = mesh_->cell_centroid(c);
    pres[c] = p0 + rho * g * (xc[dim - 1] - z0);
  }
}

void Flow_State::set_permeability_2D(double Kx, double Ky) {
  (*permeability())(0)->PutScalar(Kx);
  (*permeability())(1)->PutScalar(Ky); 
}

void Flow_State::set_permeability_3D(double Kx, double Ky, double Kz) {
  (*permeability())(0)->PutScalar(Kx);
  (*permeability())(1)->PutScalar(Ky);
  (*permeability())(2)->PutScalar(Kz);
}

void Flow_State::set_permeability_2D(double Kx, double Ky, const string region) {
  Epetra_Vector& perm_x = *(*permeability())(0);
  Epetra_Vector& perm_y = *(*permeability())(1);

  AmanziMesh::Entity_ID_List block;
  mesh_->get_set_entities(region, AmanziMesh::CELL, AmanziMesh::OWNED, &block);
  int ncells = block.size();

  for (int i=0; i!=ncells; ++i) {
    int c = block[i];
    perm_x[c] = Kx;
    perm_y[c] = Ky;
  }
}

void Flow_State::set_permeability_3D(double Kx, double Ky, double Kz, const string region) {
  Epetra_Vector& perm_x = *(*permeability())(0);
  Epetra_Vector& perm_y = *(*permeability())(1);
  Epetra_Vector& perm_z = *(*permeability())(2);

  AmanziMesh::Entity_ID_List block;
  mesh_->get_set_entities(region, AmanziMesh::CELL, AmanziMesh::OWNED, &block);
  int ncells = block.size();

  for (int i = 0; i != ncells; ++i) {
    int c = block[i];
    perm_x[c] = Kx;
    perm_y[c] = Ky;
    perm_z[c] = Kz;
  }
}

void Flow_State::set_gravity(double g) {
  Teuchos::RCP<Epetra_Vector> gvec = S_->GetConstantVectorData("gravity",name_);
  int dim = mesh_->space_dimension();

  gvec->PutScalar(0.);
  (*gvec)[dim-1] = g;
  // (*gvec)[2] = g;  // Waiting for Markus ticket (lipnikov@lanl.gov)
}

void Flow_State::set_specific_storage(double ss) {
  specific_storage()->PutScalar(ss);
}

} // namespace
} // namespace
