/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon
------------------------------------------------------------------------- */

#include "energy_bc_factory.hh"
#include "advection_factory.hh"

#include "composite_vector_function.hh"
#include "composite_vector_function_factory.hh"

#include "bdf1_time_integrator.hh"

#include "energy_base.hh"

#define MORE_DEBUG_FLAG 0


namespace Amanzi {
namespace Energy {


// -------------------------------------------------------------
// Setup
// -------------------------------------------------------------
void EnergyBase::setup(const Teuchos::Ptr<State>& S) {
  PKPhysicalBDFBase::setup(S);
  SetupEnergy_(S);
  SetupPhysicalEvaluators_(S);
};


void EnergyBase::SetupEnergy_(const Teuchos::Ptr<State>& S) {
  // Set up keys if they were not already set.
  if (energy_key_ == std::string()) {
    energy_key_ = plist_.get<std::string>("energy key",
            domain_prefix_+std::string("energy"));
  }
  if (cell_vol_key_ == std::string()) {
    cell_vol_key_ = plist_.get<std::string>("cell volume key",
            domain_prefix_+std::string("cell_volume"));
  }
  if (enthalpy_key_ == std::string()) {
    enthalpy_key_ = plist_.get<std::string>("enthalpy key",
            domain_prefix_+std::string("enthalpy"));
  }
  if (flux_key_ == std::string()) {
    flux_key_ = plist_.get<std::string>("flux key",
            domain_prefix_+std::string("flux"));
  }
  if (flux_key_ == std::string()) {
    flux_key_ = plist_.get<std::string>("energy flux key",
            domain_prefix_+std::string("energy_flux"));
  }
  if (conductivity_key_ == std::string()) {
    conductivity_key_ = plist_.get<std::string>("conductivity key",
            domain_prefix_+std::string("thermal_conductivity"));
  }
  if (de_dT_key_ == std::string()) {
    de_dT_key_ = plist_.get<std::string>("de/dT key",
            std::string("d")+energy_key_+std::string("_d")+key_);
  }

  // Require fields and evaluators for those fields.
  // primary variable: temperature on both cells and faces, ghosted, with 1 dof
  std::vector<AmanziMesh::Entity_kind> locations2(2);
  std::vector<std::string> names2(2);
  std::vector<int> num_dofs2(2,1); // = [1, 1]
  locations2[0] = AmanziMesh::CELL;
  locations2[1] = AmanziMesh::FACE;
  names2[0] = "cell";
  names2[1] = "face";

  S->RequireField(key_, name_)->SetMesh(mesh_)
    ->SetGhosted()->SetComponents(names2, locations2, num_dofs2);

#if MORE_DEBUG_FLAG
  for (int i=1; i!=23; ++i) {
    std::stringstream namestream;
    namestream << domain_prefix_ << "energy_residual_" << i;
    std::stringstream solnstream;
    solnstream << domain_prefix_ << "energy_solution_" << i;
    S->RequireField(namestream.str(), name_)->SetMesh(mesh_)
                    ->SetComponents(names2, locations2, num_dofs2);
    S->RequireField(solnstream.str(), name_)->SetMesh(mesh_)
                    ->SetComponents(names2, locations2, num_dofs2);
  }
#endif

  // Require a field and evaluator for cell volume.
  S->RequireFieldEvaluator(cell_vol_key_);

  // Require a field for the mass flux for advection.
  S->RequireField(flux_key_)->SetMesh(mesh_)->SetGhosted()
                                ->AddComponent("face", AmanziMesh::FACE, 1);

  // Require a field for the energy flux.
  std::string updatestring = plist_.get<std::string>("update flux mode", "never");
  if (updatestring == "iteration") {
    update_flux_ = UPDATE_FLUX_ITERATION;
  } else if (updatestring == "timestep") {
    update_flux_ = UPDATE_FLUX_TIMESTEP;
  } else if (updatestring == "vis") {
    update_flux_ = UPDATE_FLUX_VIS;
  } else if (updatestring == "never") {
    update_flux_ = UPDATE_FLUX_NEVER;
  } else {
    Errors::Message message(std::string("Unknown frequence for updating the overland flux: ")+updatestring);
    Exceptions::amanzi_throw(message);
  }

  // boundary conditions
  Teuchos::ParameterList bc_plist = plist_.sublist("boundary conditions", true);
  EnergyBCFactory bc_factory(mesh_, bc_plist);
  bc_temperature_ = bc_factory.CreateTemperature();
  bc_flux_ = bc_factory.CreateEnthalpyFlux();

  // coupling terms
  // -- subsurface PK, coupled to the surface
  coupled_to_surface_ = plist_.get<bool>("coupled to surface", false);
  if (coupled_to_surface_) {
    // surface temperature used for BCs
    S->RequireField("surface_temperature");
    update_flux_ = UPDATE_FLUX_ITERATION;
  }

  // flux of energy
  if (update_flux_ != UPDATE_FLUX_NEVER) {
    S->RequireField(energy_flux_key_, name_)->SetMesh(mesh_)->SetGhosted()
        ->SetComponent("face", AmanziMesh::FACE, 1);
  }

  // operator for advection terms
  Operators::AdvectionFactory advection_factory;
  Teuchos::ParameterList advect_plist = plist_.sublist("Advection");
  advection_ = advection_factory.create(advect_plist, mesh_);
  advection_->set_num_dofs(1);

  // operator for the diffusion terms
  Teuchos::ParameterList mfd_plist = plist_.sublist("Diffusion");
  matrix_ = Teuchos::rcp(new Operators::MatrixMFD(mfd_plist, mesh_));
  matrix_->SetSymmetryProperty(true);
  matrix_->SymbolicAssembleGlobalMatrices();
  matrix_->CreateMFDmassMatrices(Teuchos::null);

  // preconditioner
  Teuchos::ParameterList mfd_pc_plist = plist_.sublist("Diffusion PC");
  Teuchos::RCP<Operators::Matrix> precon =
    Teuchos::rcp(new Operators::MatrixMFD(mfd_pc_plist, mesh_));
  set_preconditioner(precon);
  assemble_preconditioner_ = plist_.get<bool>("assemble preconditioner", true);

  // constraint on max delta T, which kicks us out of bad iterates faster?
  dT_max_ = plist_.get<double>("maximum temperature change", 10.);
};


// -------------------------------------------------------------
// Initialize PK
// -------------------------------------------------------------
void EnergyBase::initialize(const Teuchos::Ptr<State>& S) {
  // initialize BDF stuff and physical domain stuff
  PKPhysicalBDFBase::initialize(S);

#if MORE_DEBUG_FLAG
  for (int i=1; i!=23; ++i) {
    std::stringstream namestream;
    namestream << domain_prefix_ << "energy_residual_" << i;
    S->GetFieldData(namestream.str(),name_)->PutScalar(0.);
    S->GetField(namestream.str(),name_)->set_initialized();

    std::stringstream solnstream;
    solnstream << domain_prefix_ << "energy_solution_" << i;
    S->GetFieldData(solnstream.str(),name_)->PutScalar(0.);
    S->GetField(solnstream.str(),name_)->set_initialized();
  }

#endif

  // initialize boundary conditions
  int nfaces = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::USED);
  bc_markers_.resize(nfaces, Operators::MATRIX_BC_NULL);
  bc_values_.resize(nfaces, 0.0);

  // initialize flux
  if (update_flux_ != UPDATE_FLUX_NEVER) {
    S->GetFieldData("darcy_flux", name_)->PutScalar(0.0);
    S->GetField("darcy_flux", name_)->set_initialized();
  }

};


// -----------------------------------------------------------------------------
// Update any secondary (dependent) variables given a solution.
//
//   After a timestep is evaluated (or at ICs), there is no way of knowing if
//   secondary variables have been updated to be consistent with the new
//   solution.
// -----------------------------------------------------------------------------
void EnergyBase::commit_state(double dt, const Teuchos::RCP<State>& S) {
  niter_ = 0;

  bool update = S->GetFieldEvaluator(conductivity_key_)
      ->HasFieldChanged(S.ptr(), name_);

  if (update_flux_ == UPDATE_FLUX_TIMESTEP ||
      (update_flux_ == UPDATE_FLUX_ITERATION && update)) {
    Teuchos::RCP<const CompositeVector> conductivity =
        S->GetFieldData(conductivity_key_);
    matrix_->CreateMFDstiffnessMatrices(conductivity.ptr());

    Teuchos::RCP<CompositeVector> temp = S->GetFieldData(key_, name_);
    Teuchos::RCP<CompositeVector> flux = S->GetFieldData(energy_flux_key_, name_);
    matrix_->DeriveFlux(*temp, flux.ptr());
  }
};


// -----------------------------------------------------------------------------
// Evaluate boundary conditions at the current time.
// -----------------------------------------------------------------------------
void EnergyBase::UpdateBoundaryConditions_() {
  for (int n=0; n!=bc_markers_.size(); ++n) {
    bc_markers_[n] = Operators::MATRIX_BC_NULL;
    bc_values_[n] = 0.0;
  }

  // Dirichlet temperature boundary conditions
  for (Functions::BoundaryFunction::Iterator bc=bc_temperature_->begin();
       bc!=bc_temperature_->end(); ++bc) {
    int f = bc->first;
    bc_markers_[f] = Operators::MATRIX_BC_DIRICHLET;
    bc_values_[f] = bc->second;
  }

  // Neumann flux boundary conditions
  for (Functions::BoundaryFunction::Iterator bc=bc_flux_->begin();
       bc!=bc_flux_->end(); ++bc) {
    int f = bc->first;
    bc_markers_[f] = Operators::MATRIX_BC_FLUX;
    bc_values_[f] = bc->second;
  }

  // Dirichlet temperature boundary conditions from a coupled surface.
  if (coupled_to_surface_) {
    // Face is Dirichlet with value of surface temp
    Teuchos::RCP<const AmanziMesh::Mesh> surface = S_next_->GetMesh("surface");
    const Epetra_MultiVector& temp = *S_next_->GetFieldData("surface_temperature")
        ->ViewComponent("cell",false);

    int ncells_surface = temp.MyLength();
    for (int c=0; c!=ncells_surface; ++c) {
      // -- get the surface cell's equivalent subsurface face
      AmanziMesh::Entity_ID f =
        surface->entity_get_parent(AmanziMesh::CELL, c);

      // -- set that value to dirichlet
      bc_markers_[f] = Operators::MATRIX_BC_DIRICHLET;
      bc_values_[f] = temp[0][c];
    }
  }

};


// -----------------------------------------------------------------------------
// Add a boundary marker to owned faces.
// -----------------------------------------------------------------------------
void EnergyBase::ApplyBoundaryConditions_(const Teuchos::RCP<CompositeVector>& temperature) {
  int nfaces = temperature->size("face");
  for (int f=0; f!=nfaces; ++f) {
    if (bc_markers_[f] == Operators::MATRIX_BC_DIRICHLET) {
      (*temperature)("face",f) = bc_values_[f];
    }
  }
};


// -----------------------------------------------------------------------------
// Check admissibility of the solution guess.
// -----------------------------------------------------------------------------
bool EnergyBase::is_admissible(Teuchos::RCP<const TreeVector> up) {
  Teuchos::OSTab tab = getOSTab();
  // For some reason, wandering PKs break most frequently with an unreasonable
  // temperature.  This simply tries to catch that before it happens.
  Teuchos::RCP<const CompositeVector> temp = up->data();

  const Epetra_MultiVector& temp_v = *temp->ViewComponent("cell",false);
  double minT(0.), maxT(0.);
  int ierr = temp_v.MinValue(&minT);
  ierr |= temp_v.MaxValue(&maxT);

  if(out_.get() && includesVerbLevel(verbosity_,Teuchos::VERB_EXTREME,true)) {
    *out_ << "Admissible T? (min/max): " << minT << ",  " << maxT << std::endl;
  }

  if (ierr || minT < 200.0 || maxT > 300.0) {
    if(out_.get() && includesVerbLevel(verbosity_,Teuchos::VERB_HIGH,true)) {
      *out_ << " is not admissible, as it is not within bounds of constitutive models: min(T) = " << minT << ", max(T) = " << maxT << std::endl;
    }
    return false;
  }
  return true;
}


// -----------------------------------------------------------------------------
// BDF takes a prediction step -- make sure it is physical and otherwise ok.
// -----------------------------------------------------------------------------
bool EnergyBase::modify_predictor(double h, const Teuchos::RCP<TreeVector>& u) {
  if (modify_predictor_with_consistent_faces_) {
    CalculateConsistentFaces(u->data().ptr());
    return true;
  }

  return PKPhysicalBDFBase::modify_predictor(h, u);
}


// -----------------------------------------------------------------------------
// Given an arbitrary set of cell values, calculate consitent face constraints.
//
//  This is useful for prediction steps, hacky preconditioners, etc.
// -----------------------------------------------------------------------------
void EnergyBase::CalculateConsistentFaces(const Teuchos::Ptr<CompositeVector>& u) {
  // update boundary conditions
  bc_temperature_->Compute(S_next_->time());
  bc_flux_->Compute(S_next_->time());
  UpdateBoundaryConditions_();

  // div K_e grad u
  S_next_->GetFieldEvaluator(conductivity_key_)
      ->HasFieldChanged(S_next_.ptr(), name_);
  Teuchos::RCP<const CompositeVector> conductivity =
      S_next_->GetFieldData(conductivity_key_);

  mfd_preconditioner_->CreateMFDstiffnessMatrices(conductivity.ptr());
  mfd_preconditioner_->CreateMFDrhsVectors();

  // skip accumulation terms, they're not needed
  // Assemble and precompute the Schur complement for inversion.
  mfd_preconditioner_->ApplyBoundaryConditions(bc_markers_, bc_values_);
  mfd_preconditioner_->AssembleGlobalMatrices();

  // derive the consistent faces, involves a solve
  mfd_preconditioner_->UpdateConsistentFaceConstraints(u.ptr());
}

} // namespace Energy
} // namespace Amanzi
