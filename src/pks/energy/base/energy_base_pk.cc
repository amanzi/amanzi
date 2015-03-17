/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon
------------------------------------------------------------------------- */

#include "energy_bc_factory.hh"
#include "advection_factory.hh"

#include "OperatorDiffusionMFD.hh"
#include "upwind_cell_centered.hh"
#include "upwind_arithmetic_mean.hh"
#include "upwind_total_flux.hh"
#include "upwind_gravity_flux.hh"

#include "composite_vector_function.hh"
#include "composite_vector_function_factory.hh"

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

  flux_tol_ = plist_->get<double>("flux tolerance", 1.);
};


void EnergyBase::SetupEnergy_(const Teuchos::Ptr<State>& S) {
  // Set up keys if they were not already set.
  if (energy_key_ == std::string()) {
    energy_key_ = plist_->get<std::string>("energy key",
            domain_prefix_+std::string("energy"));
  }
  if (cell_vol_key_ == std::string()) {
    cell_vol_key_ = plist_->get<std::string>("cell volume key",
            domain_prefix_+std::string("cell_volume"));
  }
  if (enthalpy_key_ == std::string()) {
    enthalpy_key_ = plist_->get<std::string>("enthalpy key",
            domain_prefix_+std::string("enthalpy"));
  }
  if (denthalpy_key_ == std::string()) {
    denthalpy_key_ = plist_->get<std::string>("enthalpy derivative key",
            std::string("d")+enthalpy_key_+std::string("_d")+key_);
  }
  if (flux_key_ == std::string()) {
    flux_key_ = plist_->get<std::string>("flux key",
            domain_prefix_+std::string("flux"));
  }
  if (energy_flux_key_ == std::string()) {
    energy_flux_key_ = plist_->get<std::string>("energy flux key",
            domain_prefix_+std::string("energy_flux"));
  }
  if (conductivity_key_ == std::string()) {
    conductivity_key_ = plist_->get<std::string>("conductivity key",
            domain_prefix_+std::string("thermal_conductivity"));
  }
  if (uw_conductivity_key_ == std::string()) {
    uw_conductivity_key_ = plist_->get<std::string>("upwind conductivity key",
            domain_prefix_+std::string("numerical_thermal_conductivity"));
  }
  if (de_dT_key_ == std::string()) {
    de_dT_key_ = plist_->get<std::string>("de/dT key",
            std::string("d")+energy_key_+std::string("_d")+key_);
  }
  if (source_key_ == std::string()) {
    source_key_ = plist_->get<std::string>("source key",
            domain_prefix_+std::string("total_energy_source"));
  }
  if (dsource_dT_key_ == std::string()) {
    dsource_dT_key_ = std::string("d")+source_key_+std::string("_d")+key_;
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
  S->RequireField(cell_vol_key_)->SetMesh(mesh_)
      ->AddComponent("cell", AmanziMesh::CELL, 1);
  S->RequireFieldEvaluator(cell_vol_key_);

  // Require a field for the mass flux for advection.
  S->RequireField(flux_key_)->SetMesh(mesh_)->SetGhosted()
      ->AddComponent("face", AmanziMesh::FACE, 1);

  // Require a field for the (conducted) energy flux.
  std::string updatestring = plist_->get<std::string>("update flux mode", "vis");
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
  S->RequireField(energy_flux_key_, name_)->SetMesh(mesh_)->SetGhosted()
      ->SetComponent("face", AmanziMesh::FACE, 1);

  // Require an upwinding strategy
  S->RequireField(uw_conductivity_key_, name_)->SetMesh(mesh_)->SetGhosted()
                    ->SetComponents(names2, locations2, num_dofs2);
  S->GetField(uw_conductivity_key_,name_)->set_io_vis(false);
  std::string method_name = plist_->get<std::string>("upwind conductivity method", "arithmetic mean");
  if (method_name == "cell centered") {
    upwinding_ = Teuchos::rcp(new Operators::UpwindCellCentered(name_,
            conductivity_key_, uw_conductivity_key_));
    Krel_method_ = Operators::UPWIND_METHOD_CENTERED;
  } else if (method_name == "arithmetic mean") {
    upwinding_ = Teuchos::rcp(new Operators::UpwindArithmeticMean(name_,
            conductivity_key_, uw_conductivity_key_));
    Krel_method_ = Operators::UPWIND_METHOD_ARITHMETIC_MEAN;
  } else {
    std::stringstream messagestream;
    messagestream << "Energy PK has no upwinding method named: " << method_name;
    Errors::Message message(messagestream.str());
    Exceptions::amanzi_throw(message);
  }


  // coupling terms
  // -- subsurface PK, coupled to the surface
  coupled_to_surface_via_flux_ = plist_->get<bool>("coupled to surface via flux", false);
  if (coupled_to_surface_via_flux_) {
    S->RequireField("surface_subsurface_energy_flux")
        ->SetMesh(S->GetMesh("surface"))
        ->AddComponent("cell", AmanziMesh::CELL, 1);
  }

  coupled_to_surface_via_temp_ =
      plist_->get<bool>("coupled to surface via temperature", false);
  if (coupled_to_surface_via_temp_) {
    // surface temperature used for BCs
    S->RequireField("surface_temperature");
    update_flux_ = UPDATE_FLUX_ITERATION;
  }

  // source terms
  is_source_term_ = plist_->get<bool>("source term");
  if (is_source_term_) {
    S->RequireField(source_key_)->SetMesh(mesh_)
        ->AddComponent("cell", AmanziMesh::CELL, 1);
    S->RequireFieldEvaluator(source_key_);
  }

  // boundary conditions
  Teuchos::ParameterList bc_plist = plist_->sublist("boundary conditions", true);
  EnergyBCFactory bc_factory(mesh_, bc_plist);
  bc_temperature_ = bc_factory.CreateTemperature();
  bc_flux_ = bc_factory.CreateEnthalpyFlux();

  int nfaces = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::USED);
  bc_markers_.resize(nfaces, Operators::OPERATOR_BC_NONE);
  bc_values_.resize(nfaces, 0.0);
  std::vector<double> mixed;
  bc_ = Teuchos::rcp(new Operators::BCs(Operators::OPERATOR_BC_TYPE_FACE, bc_markers_, bc_values_, mixed));

  bc_markers_adv_.resize(nfaces, Operators::OPERATOR_BC_NONE);
  bc_values_adv_.resize(nfaces, 0.0);
  bc_adv_ = Teuchos::rcp(new Operators::BCs(Operators::OPERATOR_BC_TYPE_FACE,
          bc_markers_adv_, bc_values_adv_, mixed));

  // Operator:
  //  -- operator for the diffusion terms
  Teuchos::ParameterList mfd_plist = plist_->sublist("Diffusion");
  matrix_diff_ = Teuchos::rcp(new Operators::OperatorDiffusionMFD(mfd_plist, mesh_));
  matrix_diff_->SetBCs(bc_);
  matrix_diff_->Setup(Teuchos::null);

  //  -- for advection terms
  implicit_advection_ = !plist_->get<bool>("explicit advection", false);
  if (implicit_advection_)
    implicit_advection_in_pc_ = !plist_->get<bool>("supress advective terms in preconditioner", false);
  Operators::AdvectionFactory advection_factory;
  Teuchos::ParameterList advect_plist = plist_->sublist("Advection");
  matrix_adv_ = Teuchos::rcp(new Operators::OperatorAdvection(advect_plist, mesh_));

  // Preconditioner
  // -- for the diffusive terms
  Teuchos::ParameterList mfd_pc_plist = plist_->sublist("Diffusion PC");
  preconditioner_diff_ = Teuchos::rcp(new Operators::OperatorDiffusionMFD(mfd_pc_plist, mesh_));
  preconditioner_diff_->SetBCs(bc_);
  preconditioner_diff_->Setup(Teuchos::null);
  preconditioner_ = preconditioner_diff_->global_operator();

  // -- for the accumulation terms
  preconditioner_acc_ = Teuchos::rcp(new Operators::OperatorAccumulation(AmanziMesh::CELL, preconditioner_));

  if (implicit_advection_ && implicit_advection_in_pc_) {
    preconditioner_adv_ = Teuchos::rcp(new Operators::OperatorAdvection(advect_plist, preconditioner_));
  }

  precon_used_ = mfd_pc_plist.isSublist("preconditioner");
  if (precon_used_) {
    preconditioner_->SymbolicAssembleMatrix();
  }

  // constraint on max delta T, which kicks us out of bad iterates faster?
  dT_max_ = plist_->get<double>("maximum temperature change", 10.);

  // ewc and other predictors can result in odd face values
  modify_predictor_with_consistent_faces_ = false;
  //    plist_->get<bool>("modify predictor with consistent faces", false);
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

  // Set extra fields as initialized -- these don't currently have evaluators,
  // and will be initialized in the call to commit_state()
  S->GetFieldData(uw_conductivity_key_,name_)->PutScalar(1.0);
  S->GetField(uw_conductivity_key_,name_)->set_initialized();

  // initialize flux
  S->GetFieldData(energy_flux_key_, name_)->PutScalar(0.0);
  S->GetField(energy_flux_key_, name_)->set_initialized();
};


// -----------------------------------------------------------------------------
// Update any secondary (dependent) variables given a solution.
//
//   After a timestep is evaluated (or at ICs), there is no way of knowing if
//   secondary variables have been updated to be consistent with the new
//   solution.
// -----------------------------------------------------------------------------
void EnergyBase::commit_state(double dt, const Teuchos::RCP<State>& S) {
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_EXTREME))
    *vo_->os() << "Commiting state." << std::endl;
  PKPhysicalBDFBase::commit_state(dt, S);

  niter_ = 0;
  bool update = UpdateConductivityData_(S.ptr());

  if (update_flux_ == UPDATE_FLUX_TIMESTEP ||
      (update_flux_ == UPDATE_FLUX_ITERATION && update)) {
    Teuchos::RCP<const CompositeVector> conductivity =
        S->GetFieldData(uw_conductivity_key_);
    matrix_diff_->global_operator()->Init();
    matrix_diff_->Setup(conductivity, Teuchos::null);
    matrix_diff_->UpdateMatrices(Teuchos::null, Teuchos::null);

    Teuchos::RCP<const CompositeVector> temp = S->GetFieldData(key_);
    Teuchos::RCP<CompositeVector> flux = S->GetFieldData(energy_flux_key_, name_);
    matrix_diff_->UpdateFlux(*temp, *flux);
  }
};


bool EnergyBase::UpdateConductivityData_(const Teuchos::Ptr<State>& S) {
  bool update = S->GetFieldEvaluator(conductivity_key_)->HasFieldChanged(S, name_);
  if (update) {
    upwinding_->Update(S);
  }
  return update;
}

// -----------------------------------------------------------------------------
// Evaluate boundary conditions at the current time.
// -----------------------------------------------------------------------------
void EnergyBase::UpdateBoundaryConditions_() {
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_EXTREME))
    *vo_->os() << "  Updating BCs." << std::endl;

  for (unsigned int n=0; n!=bc_markers_.size(); ++n) {
    bc_markers_[n] = Operators::OPERATOR_BC_NONE;
    bc_values_[n] = 0.0;
    bc_markers_adv_[n] = Operators::OPERATOR_BC_NONE;
    bc_values_adv_[n] = 0.0;
  }

  // Dirichlet temperature boundary conditions
  for (Functions::BoundaryFunction::Iterator bc=bc_temperature_->begin();
       bc!=bc_temperature_->end(); ++bc) {
    int f = bc->first;
    bc_markers_[f] = Operators::OPERATOR_BC_DIRICHLET;
    bc_values_[f] = bc->second;
    bc_markers_adv_[f] = Operators::OPERATOR_BC_DIRICHLET;
  }

  // Neumann flux boundary conditions
  for (Functions::BoundaryFunction::Iterator bc=bc_flux_->begin();
       bc!=bc_flux_->end(); ++bc) {
    int f = bc->first;
    bc_markers_[f] = Operators::OPERATOR_BC_NEUMANN;
    bc_values_[f] = bc->second;
    bc_markers_adv_[f] = Operators::OPERATOR_BC_NEUMANN;
    bc_values_adv_[f] = 0.;
    // push all onto diffusion
  }

  // Dirichlet temperature boundary conditions from a coupled surface.
  if (coupled_to_surface_via_temp_) {
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
      bc_markers_[f] = Operators::OPERATOR_BC_DIRICHLET;
      bc_values_[f] = temp[0][c];
      bc_markers_adv_[f] = Operators::OPERATOR_BC_DIRICHLET;

    }
  }

  // surface coupling
  if (coupled_to_surface_via_flux_) {
    // Diffusive fluxes are given by the residual of the surface equation.
    // Advective fluxes are given by the surface temperature and whatever flux we have.
    Teuchos::RCP<const AmanziMesh::Mesh> surface = S_next_->GetMesh("surface");
    const Epetra_MultiVector& flux =
        *S_next_->GetFieldData("surface_subsurface_energy_flux")
        ->ViewComponent("cell",false);

    int ncells_surface = flux.MyLength();
    for (int c=0; c!=ncells_surface; ++c) {
      // -- get the surface cell's equivalent subsurface face
      AmanziMesh::Entity_ID f =
        surface->entity_get_parent(AmanziMesh::CELL, c);

      // -- set that value to Neumann
      bc_markers_[f] = Operators::OPERATOR_BC_NEUMANN;
      // flux is in units of J / s, whereas Neumann BCs are J/s/A
      bc_values_[f] = flux[0][c] / mesh_->face_area(f);

      // -- mark advective BCs as Dirichlet: this ensures the surface
      //    temperature is picked up and advective fluxes are treated
      //    via advection operator, not diffusion operator.
      bc_markers_adv_[f] = Operators::OPERATOR_BC_DIRICHLET;
    }
  }
};


// -----------------------------------------------------------------------------
// Add a boundary marker to owned faces.
// -----------------------------------------------------------------------------
void EnergyBase::ApplyBoundaryConditions_(const Teuchos::RCP<CompositeVector>& temp) {
  Epetra_MultiVector& temp_f = *temp->ViewComponent("face",true);
  unsigned int nfaces = temp->size("face",true);
  for (unsigned int f=0; f!=nfaces; ++f) {
    if (bc_markers_[f] == Operators::OPERATOR_BC_DIRICHLET) {
      temp_f[0][f] = bc_values_[f];
    }
  }
};


// -----------------------------------------------------------------------------
// Check admissibility of the solution guess.
// -----------------------------------------------------------------------------
bool EnergyBase::IsAdmissible(Teuchos::RCP<const TreeVector> up) {
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_EXTREME))
    *vo_->os() << "  Checking admissibility..." << std::endl;

  // For some reason, wandering PKs break most frequently with an unreasonable
  // temperature.  This simply tries to catch that before it happens.
  Teuchos::RCP<const CompositeVector> temp = up->Data();
  double minT, maxT;
  
  const Epetra_MultiVector& temp_c = *temp->ViewComponent("cell",false);
  double minT_c(1.e6), maxT_c(-1.e6);
  int min_c(-1), max_c(-1);
  for (int c=0; c!=temp_c.MyLength(); ++c) {
    if (temp_c[0][c] < minT_c) {
      minT_c = temp_c[0][c];
      min_c = c;
    }
    if (temp_c[0][c] > maxT_c) {
      maxT_c = temp_c[0][c];
      max_c = c;
    }
  }

  double minT_f(1.e6), maxT_f(-1.e6);
  int min_f(-1), max_f(-1);
  if (temp->HasComponent("face")) {
    const Epetra_MultiVector& temp_f = *temp->ViewComponent("face",false);
    for (int f=0; f!=temp_f.MyLength(); ++f) {
      if (temp_f[0][f] < minT_f) {
        minT_f = temp_f[0][f];
        min_f = f;
      }
      if (temp_f[0][f] > maxT_f) {
        maxT_f = temp_f[0][f];
        max_f = f;
      }
    }
    minT = std::min(minT_c, minT_f);
    maxT = std::max(maxT_c, maxT_f);

  } else {
    minT = minT_c;
    maxT = maxT_c;
  }

  double minT_l = minT;
  double maxT_l = maxT;
  mesh_->get_comm()->MaxAll(&maxT_l, &maxT, 1);
  mesh_->get_comm()->MinAll(&minT_l, &minT, 1);
  
  if (vo_->os_OK(Teuchos::VERB_HIGH)) {
    *vo_->os() << "    Admissible T? (min/max): " << minT << ",  " << maxT << std::endl;
  }

  if (minT < 200.0 || maxT > 300.0) {
    if (vo_->os_OK(Teuchos::VERB_MEDIUM)) {
      *vo_->os() << " is not admissible, as it is not within bounds of constitutive models:" << std::endl;
      ENorm_t global_minT_c, local_minT_c;
      ENorm_t global_maxT_c, local_maxT_c;

      local_minT_c.value = minT_c;
      local_minT_c.gid = temp_c.Map().GID(min_c);
      local_maxT_c.value = maxT_c;
      local_maxT_c.gid = temp_c.Map().GID(max_c);

      MPI_Allreduce(&local_minT_c, &global_minT_c, 1, MPI_DOUBLE_INT, MPI_MINLOC, MPI_COMM_WORLD);
      MPI_Allreduce(&local_maxT_c, &global_maxT_c, 1, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD);
      *vo_->os() << "   cells (min/max): [" << global_minT_c.gid << "] " << global_minT_c.value
                 << ", [" << global_maxT_c.gid << "] " << global_maxT_c.value << std::endl;

      if (temp->HasComponent("face")) {
        const Epetra_MultiVector& temp_f = *temp->ViewComponent("face",false);
        ENorm_t global_minT_f, local_minT_f;
        ENorm_t global_maxT_f, local_maxT_f;

        local_minT_f.value = minT_f;
        local_minT_f.gid = temp_f.Map().GID(min_f);
        local_maxT_f.value = maxT_f;
        local_maxT_f.gid = temp_f.Map().GID(max_f);
        
        MPI_Allreduce(&local_minT_f, &global_minT_f, 1, MPI_DOUBLE_INT, MPI_MINLOC, MPI_COMM_WORLD);
        MPI_Allreduce(&local_maxT_f, &global_maxT_f, 1, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD);
        *vo_->os() << "   cells (min/max): [" << global_minT_f.gid << "] " << global_minT_f.value
                   << ", [" << global_maxT_f.gid << "] " << global_maxT_f.value << std::endl;
      }
    }
    return false;
  }
  return true;
}


// -----------------------------------------------------------------------------
// BDF takes a prediction step -- make sure it is physical and otherwise ok.
// -----------------------------------------------------------------------------
bool EnergyBase::ModifyPredictor(double h, Teuchos::RCP<const TreeVector> u0,
        Teuchos::RCP<TreeVector> u) {

  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_EXTREME))
    *vo_->os() << "Modifying predictor:" << std::endl;

  if (modify_predictor_with_consistent_faces_) {
    if (vo_->os_OK(Teuchos::VERB_EXTREME))
      *vo_->os() << "  modifications for consistent face temperatures." << std::endl;
    //    CalculateConsistentFaces(u->Data().ptr());
    return true;
  }

  return false;
}


// // -----------------------------------------------------------------------------
// // Given an arbitrary set of cell values, calculate consitent face constraints.
// //
// //  This is useful for prediction steps, hacky preconditioners, etc.
// // -----------------------------------------------------------------------------
// void EnergyBase::CalculateConsistentFaces(const Teuchos::Ptr<CompositeVector>& u) {
//   // update boundary conditions
//   bc_temperature_->Compute(S_next_->time());
//   bc_flux_->Compute(S_next_->time());
//   UpdateBoundaryConditions_();

//   // div K_e grad u
//   ChangedSolution();
//   bool update = UpdateConductivityData_(S_next_.ptr());
//   Teuchos::RCP<const CompositeVector> conductivity =
//       S_next_->GetFieldData(uw_conductivity_key_);

//   matrix_->CreateMFDstiffnessMatrices(conductivity.ptr());
//   matrix_->CreateMFDrhsVectors();

//   // skip accumulation terms, they're not needed
//   // Assemble and precompute the Schur complement for inversion.
//   matrix_->ApplyBoundaryConditions(bc_markers_, bc_values_);

//   // derive the consistent faces, involves a solve
//   matrix_->UpdateConsistentFaceConstraints(u.ptr());
// }

} // namespace Energy
} // namespace Amanzi
