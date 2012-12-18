/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
This is the flow component of the Amanzi code.
License: BSD
Authors: Neil Carlson (version 1)
         Konstantin Lipnikov (version 2) (lipnikov@lanl.gov)
         Ethan Coon (ATS version) (ecoon@lanl.gov)
------------------------------------------------------------------------- */
#include "boost/math/special_functions/fpclassify.hpp"

#include "bdf1_time_integrator.hh"
#include "flow_bc_factory.hh"

#include "upwinding.hh"
#include "Mesh_MSTK.hh"
#include "Point.hh"
#include "matrix_mfd.cc"

#include "upwind_cell_centered.hh"
#include "upwind_arithmetic_mean.hh"
#include "upwind_total_flux.hh"
#include "upwind_gravity_flux.hh"

#include "composite_vector_function.hh"
#include "composite_vector_function_factory.hh"

#include "wrm_richards_evaluator.hh"
#include "rel_perm_evaluator.hh"
#include "richards_water_content.hh"

#include "richards.hh"

#define DEBUG_FLAG 0


namespace Amanzi {
namespace Flow {

RegisteredPKFactory<Richards> Richards::reg_("richards flow");


// -------------------------------------------------------------
// Setup data
// -------------------------------------------------------------
void Richards::setup(const Teuchos::Ptr<State>& S) {
  PKPhysicalBDFBase::setup(S);
  SetupRichardsFlow_(S);
  SetupPhysicalEvaluators_(S);
};


// -------------------------------------------------------------
// Pieces of the construction process that are common to all
// Richards-like PKs.
// -------------------------------------------------------------
void Richards::SetupRichardsFlow_(const Teuchos::Ptr<State>& S) {

  // Require fields and evaluators for those fields.
  // -- primary variable: pressure on both cells and faces, ghosted, with 1 dof
  std::vector<AmanziMesh::Entity_kind> locations2(2);
  std::vector<std::string> names2(2);
  std::vector<int> num_dofs2(2,1);
  locations2[0] = AmanziMesh::CELL;
  locations2[1] = AmanziMesh::FACE;
  names2[0] = "cell";
  names2[1] = "face";

  S->RequireField(key_, name_)->SetMesh(mesh_)->SetGhosted()
                    ->SetComponents(names2, locations2, num_dofs2);

#if DEBUG_FLAG
  for (int i=1; i!=23; ++i) {
    std::stringstream namestream;
    namestream << "flow_residual_" << i;
    std::stringstream solnstream;
    solnstream << "flow_solution_" << i;
    S->RequireField(namestream.str(), name_)->SetMesh(mesh_)->SetGhosted()
                    ->SetComponents(names2, locations2, num_dofs2);
    S->RequireField(solnstream.str(), name_)->SetMesh(mesh_)->SetGhosted()
                    ->SetComponents(names2, locations2, num_dofs2);
  }
#endif

  // -- secondary variables, no evaluator used
  S->RequireField("darcy_flux_direction", name_)->SetMesh(mesh_)->SetGhosted()
      ->SetComponent("face", AmanziMesh::FACE, 1);
  S->RequireField("total_flux", name_)->SetMesh(mesh_)->SetGhosted()
                                ->SetComponent("face", AmanziMesh::FACE, 1);
  S->RequireField("darcy_flux", name_)->SetMesh(mesh_)->SetGhosted()
                                ->SetComponent("face", AmanziMesh::FACE, 1);
  S->RequireField("total_velocity", name_)->SetMesh(mesh_)->SetGhosted()
                                ->SetComponent("cell", AmanziMesh::CELL, 3);
  S->RequireField("darcy_velocity", name_)->SetMesh(mesh_)->SetGhosted()
                                ->SetComponent("cell", AmanziMesh::CELL, 3);

  // Get data for non-field quanitites.
  S->RequireFieldEvaluator("cell_volume");
  S->RequireGravity();
  S->RequireScalar("atmospheric_pressure");

  // Create the absolute permeability tensor.
  int c_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  K_ = Teuchos::rcp(new std::vector<WhetStone::Tensor>(c_owned));
  for (int c=0; c!=c_owned; ++c) {
    (*K_)[c].init(mesh_->space_dimension(),1);
  }

  // Create the boundary condition data structures.
  Teuchos::ParameterList bc_plist = plist_.sublist("boundary conditions", true);
  FlowBCFactory bc_factory(mesh_, bc_plist);
  bc_pressure_ = bc_factory.CreatePressure();
  bc_flux_ = bc_factory.CreateMassFlux();
  infiltrate_only_if_unfrozen_ = bc_plist.get<bool>("infiltrate only if unfrozen",false);
  bc_seepage_ = bc_factory.CreateSeepageFace();
  bc_seepage_->Compute(0.); // compute at t=0 to set up

  // coupling
  coupled_to_surface_via_residual_ = plist_.get<bool>("coupled to surface via residual", false);
  if (coupled_to_surface_via_residual_) {
    S->RequireField("ponded_depth");
    S->RequireField("overland_source_from_subsurface", name_)
        ->SetMesh(S->GetMesh("surface"))->SetComponent("cell", AmanziMesh::CELL, 1);

    surface_head_eps_ = plist_.get<double>("surface head epsilon", 0.);
  }

  // Create the upwinding method
  S->RequireField("numerical_rel_perm", name_)->SetMesh(mesh_)->SetGhosted()
                    ->SetComponents(names2, locations2, num_dofs2);
  S->GetField("numerical_rel_perm",name_)->set_io_vis(false);
  string method_name = plist_.get<string>("relative permeability method", "upwind with gravity");
  bool symmetric = false;
  if (method_name == "upwind with gravity") {
    upwinding_ = Teuchos::rcp(new Operators::UpwindGravityFlux(name_,
            "relative_permeability", "numerical_rel_perm", K_));
    Krel_method_ = FLOW_RELATIVE_PERM_UPWIND_GRAVITY;
  } else if (method_name == "cell centered") {
    upwinding_ = Teuchos::rcp(new Operators::UpwindCellCentered(name_,
            "relative_permeability", "numerical_rel_perm"));
    symmetric = true;
    Krel_method_ = FLOW_RELATIVE_PERM_CENTERED;
  } else if (method_name == "upwind with Darcy flux") {
    upwinding_ = Teuchos::rcp(new Operators::UpwindTotalFlux(name_,
            "relative_permeability", "numerical_rel_perm", "darcy_flux_direction"));
    Krel_method_ = FLOW_RELATIVE_PERM_UPWIND_DARCY_FLUX;
  } else if (method_name == "arithmetic mean") {
    upwinding_ = Teuchos::rcp(new Operators::UpwindArithmeticMean(name_,
            "relative_permeability", "numerical_rel_perm"));
    Krel_method_ = FLOW_RELATIVE_PERM_ARITHMETIC_MEAN;
  } else {
    std::stringstream messagestream;
    messagestream << "Richards FLow PK has no upwinding method named: " << method_name;
    Errors::Message message(messagestream.str());
    Exceptions::amanzi_throw(message);
  }

  // operator for the diffusion terms
  Teuchos::ParameterList mfd_plist = plist_.sublist("Diffusion");
  matrix_ = Teuchos::rcp(new Operators::MatrixMFD(mfd_plist, mesh_));
  matrix_->SetSymmetryProperty(symmetric);
  matrix_->SymbolicAssembleGlobalMatrices();

  // preconditioner for the NKA system
  Teuchos::ParameterList mfd_pc_plist = plist_.sublist("Diffusion PC");
  preconditioner_ = Teuchos::rcp(new Operators::MatrixMFD(mfd_pc_plist, mesh_));
  preconditioner_->SetSymmetryProperty(symmetric);
  preconditioner_->SymbolicAssembleGlobalMatrices();
  preconditioner_->InitPreconditioner(mfd_pc_plist);
  assemble_preconditioner_ = plist_.get<bool>("assemble preconditioner", true);
  modify_predictor_with_consistent_faces_ =
    plist_.get<bool>("modify predictor with consistent faces", false);
}

// -------------------------------------------------------------
// Create the physical evaluators for water content, water
// retention, rel perm, etc, that are specific to Richards.
// -------------------------------------------------------------
void Richards::SetupPhysicalEvaluators_(const Teuchos::Ptr<State>& S) {
  // -- Absolute permeability.
  //       For now, we assume scalar permeability.  This will change.
  S->RequireField("permeability")->SetMesh(mesh_)->SetGhosted()
      ->AddComponent("cell", AmanziMesh::CELL, 1);
  S->RequireFieldEvaluator("permeability");

  // -- water content, and evaluator
  S->RequireField("water_content")->SetMesh(mesh_)->SetGhosted()
      ->AddComponent("cell", AmanziMesh::CELL, 1);
  Teuchos::ParameterList wc_plist = plist_.sublist("water content evaluator");
  Teuchos::RCP<RichardsWaterContent> wc = Teuchos::rcp(new RichardsWaterContent(wc_plist));
  S->SetFieldEvaluator("water_content", wc);

  // -- Water retention evaluators, for saturation and rel perm.
  S->RequireField("relative_permeability")->SetMesh(mesh_)->SetGhosted()
      ->AddComponent("cell", AmanziMesh::CELL, 1);
  Teuchos::ParameterList wrm_plist = plist_.sublist("water retention evaluator");
  Teuchos::RCP<FlowRelations::WRMRichardsEvaluator> wrm =
      Teuchos::rcp(new FlowRelations::WRMRichardsEvaluator(wrm_plist));
  S->SetFieldEvaluator("saturation_liquid", wrm);
  S->SetFieldEvaluator("saturation_gas", wrm);

  Teuchos::RCP<FlowRelations::RelPermEvaluator> rel_perm_evaluator =
      Teuchos::rcp(new FlowRelations::RelPermEvaluator(wrm_plist, wrm->get_WRMs()));
  S->SetFieldEvaluator("relative_permeability", rel_perm_evaluator);

  // -- Liquid density and viscosity for the transmissivity.
  S->RequireField("molar_density_liquid")->SetMesh(mesh_)->SetGhosted()
      ->AddComponent("cell", AmanziMesh::CELL, 1);
  S->RequireFieldEvaluator("molar_density_liquid");

  S->RequireField("viscosity_liquid")->SetMesh(mesh_)->SetGhosted()
      ->AddComponent("cell", AmanziMesh::CELL, 1);
  S->RequireFieldEvaluator("viscosity_liquid");

  // -- liquid mass density for the gravity fluxes
  S->RequireField("mass_density_liquid")->SetMesh(mesh_)->SetGhosted()
      ->AddComponent("cell", AmanziMesh::CELL, 1);
  S->RequireFieldEvaluator("mass_density_liquid"); // simply picks up the molar density one.
}



// -------------------------------------------------------------
// Initialize PK
// -------------------------------------------------------------
void Richards::initialize(const Teuchos::Ptr<State>& S) {
  // Check for PK-specific initialization
  if (!plist_.isSublist("initial condition")) {
    std::stringstream messagestream;
    messagestream << name_ << " has no initial condition parameter list.";
    Errors::Message message(messagestream.str());
    Exceptions::amanzi_throw(message);
  }
  Teuchos::ParameterList ic_plist = plist_.sublist("initial condition");
  if (ic_plist.isParameter("initial water table z-coordinate")) {
    S->GetFieldData(key_,name_)->PutScalar(101325.);
    S->GetField(key_, name_)->set_initialized();
  }

  // initialize BDF time integrator (BDFBase) and primary variable (PhysicalBase)
  PKPhysicalBDFBase::initialize(S);

  S->GetFieldData(key_)->ScatterMasterToGhosted("face");

  // debugggin cruft
#if DEBUG_FLAG
  for (int i=1; i!=23; ++i) {
    std::stringstream namestream;
    namestream << "flow_residual_" << i;
    S->GetFieldData(namestream.str(),name_)->PutScalar(0.);
    S->GetField(namestream.str(),name_)->set_initialized();

    std::stringstream solnstream;
    solnstream << "flow_solution_" << i;
    S->GetFieldData(solnstream.str(),name_)->PutScalar(0.);
    S->GetField(solnstream.str(),name_)->set_initialized();
  }

#endif

  // initialize boundary conditions
  int nfaces = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::USED);
  bc_markers_.resize(nfaces, Operators::MFD_BC_NULL);
  bc_values_.resize(nfaces, 0.0);

  // Set extra fields as initialized -- these don't currently have evaluators,
  // and will be initialized in the call to commit_state()
  S->GetFieldData("numerical_rel_perm",name_)->PutScalar(1.0);
  S->GetField("numerical_rel_perm",name_)->set_initialized();

  S->GetFieldData("darcy_flux", name_)->PutScalar(0.0);
  S->GetField("darcy_flux", name_)->set_initialized();
  S->GetFieldData("total_flux", name_)->PutScalar(0.0);
  S->GetField("total_flux", name_)->set_initialized();
  S->GetFieldData("darcy_flux_direction", name_)->PutScalar(0.0);
  S->GetField("darcy_flux_direction", name_)->set_initialized();

  S->GetFieldData("darcy_velocity", name_)->PutScalar(0.0);
  S->GetField("darcy_velocity", name_)->set_initialized();
  S->GetFieldData("total_velocity", name_)->PutScalar(0.0);
  S->GetField("total_velocity", name_)->set_initialized();

  if (coupled_to_surface_via_residual_) {
    S->GetFieldData("overland_source_from_subsurface", name_)->PutScalar(0.);
    S->GetField("overland_source_from_subsurface", name_)->set_initialized();
  }

  // absolute perm
  SetAbsolutePermeabilityTensor_(S);

  // operators
  matrix_->CreateMFDmassMatrices(K_.ptr());
  preconditioner_->CreateMFDmassMatrices(K_.ptr());
};


// -----------------------------------------------------------------------------
// Update any secondary (dependent) variables given a solution.
//
//   After a timestep is evaluated (or at ICs), there is no way of knowing if
//   secondary variables have been updated to be consistent with the new
//   solution.
// -----------------------------------------------------------------------------
void Richards::commit_state(double dt, const Teuchos::RCP<State>& S) {
  niter_ = 0;

  // Update flux if rel perm, density, or pressure have changed.
  UpdateFlux_(S);

  // As a diagnostic, calculate the mass balance error
  if (S_next_ != Teuchos::null) {
    Teuchos::RCP<const CompositeVector> wc1 = S_next_->GetFieldData("water_content");
    Teuchos::RCP<const CompositeVector> wc0 = S_->GetFieldData("water_content");
    Teuchos::RCP<const CompositeVector> darcy_flux = S->GetFieldData("darcy_flux", name_);
    CompositeVector error(*wc1);

    for (int c=0; c!=error.size("cell"); ++c) {
      error("cell",c) = (*wc1)("cell",c) - (*wc0)("cell",c);

      AmanziMesh::Entity_ID_List faces;
      std::vector<int> dirs;
      mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
      for (int n=0; n!=faces.size(); ++n) {
        error("cell",c) += (*darcy_flux)("face",faces[n]) * dirs[n] * dt;
      }
    }

    double einf(0.0);
    error.NormInf(&einf);

    // VerboseObject stuff.
    Teuchos::OSTab tab = getOSTab();
    *out_ << "Final Mass Balance Error: " << einf << std::endl;
  }
};


// -----------------------------------------------------------------------------
// Update any diagnostic variables prior to vis (in this case velocity field).
// -----------------------------------------------------------------------------
void Richards::calculate_diagnostics(const Teuchos::RCP<State>& S) {
  // update the cell velocities
  Teuchos::RCP<CompositeVector> darcy_velocity = S->GetFieldData("darcy_velocity", name_);
  Teuchos::RCP<const CompositeVector> darcy_flux = S->GetFieldData("darcy_flux");
  matrix_->DeriveCellVelocity(*darcy_flux, darcy_velocity);

  Teuchos::RCP<CompositeVector> total_velocity = S->GetFieldData("total_velocity", name_);
  Teuchos::RCP<const CompositeVector> total_flux = S->GetFieldData("total_flux");
  matrix_->DeriveCellVelocity(*total_flux, total_velocity);
};


// -----------------------------------------------------------------------------
// Use the physical rel perm (on cells) to update a work vector for rel perm.
//
//   This deals with upwinding, etc.
// -----------------------------------------------------------------------------
bool Richards::UpdatePermeabilityData_(const Teuchos::Ptr<State>& S) {
  Teuchos::RCP<CompositeVector> uw_rel_perm = S->GetFieldData("numerical_rel_perm", name_);
  Teuchos::RCP<const CompositeVector> rel_perm = S->GetFieldData("relative_permeability");
  bool update_perm = S->GetFieldEvaluator("relative_permeability")
      ->HasFieldChanged(S, name_);

  update_perm |= S->GetFieldEvaluator("molar_density_liquid")->HasFieldChanged(S, name_);
  update_perm |= S->GetFieldEvaluator("viscosity_liquid")->HasFieldChanged(S, name_);

  // requirements due to the upwinding method
  if (Krel_method_ == FLOW_RELATIVE_PERM_UPWIND_DARCY_FLUX) {
    bool update_dir = S->GetFieldEvaluator("mass_density_liquid")
        ->HasFieldChanged(S, name_);
    update_dir |= S->GetFieldEvaluator(key_)->HasFieldChanged(S, name_);

    if (update_dir) {
      // update the direction of the flux -- note this is NOT the flux
      Teuchos::RCP<CompositeVector> flux_dir =
          S->GetFieldData("darcy_flux_direction", name_);

      // Create the stiffness matrix without a rel perm (rel perm = 1)
      matrix_->CreateMFDstiffnessMatrices(Teuchos::null);

      // Derive the pressure fluxes
      Teuchos::RCP<const CompositeVector> pres = S->GetFieldData(key_);
      matrix_->DeriveFlux(*pres, flux_dir);

      // Add in the gravity fluxes
      Teuchos::RCP<const Epetra_Vector> gvec = S->GetConstantVectorData("gravity");
      Teuchos::RCP<const CompositeVector> rho = S->GetFieldData("mass_density_liquid");
      AddGravityFluxesToVector_(gvec, Teuchos::null, rho, flux_dir);
    }


    update_perm |= update_dir;
  }

  if (update_perm) {
    // upwind
    upwinding_->Update(S);

    // patch up the BCs -- FIX ME --etc
    for (int f=0; f!=uw_rel_perm->size("face",false); ++f) {
      AmanziMesh::Entity_ID_List cells;
      uw_rel_perm->mesh()->face_get_cells(f, AmanziMesh::USED, &cells);
      if (cells.size() == 1) {
        // just grab the cell inside's perm... this will need to be fixed eventually.
        (*uw_rel_perm)("face",f) = (*rel_perm)("cell",cells[0]);
      }
    }

    if (coupled_to_surface_via_residual_) {
      // patch up the rel perm on surface as 1 -- FIX ME --etc
      Teuchos::RCP<const AmanziMesh::Mesh_MSTK> surface =
          Teuchos::rcp_static_cast<const AmanziMesh::Mesh_MSTK>(S->GetMesh("surface"));
      int ncells_surface = surface->num_entities(AmanziMesh::CELL,AmanziMesh::OWNED);
      for (int c=0; c!=ncells_surface; ++c) {
        // -- get the surface cell's equivalent subsurface face and neighboring cell
        AmanziMesh::Entity_ID f =
            surface->entity_get_parent(AmanziMesh::CELL, c);
        (*uw_rel_perm)("face",f) = 1.0;
      }
    }
  }

  // Scale cells by n/visc if needed.
  if (update_perm) {
    Teuchos::RCP<const CompositeVector> n_liq = S->GetFieldData("molar_density_liquid");
    Teuchos::RCP<const CompositeVector> visc = S->GetFieldData("viscosity_liquid");
    for (int c=0; c!=uw_rel_perm->size("cell", false); ++c) {
      (*uw_rel_perm)("cell",c) *= (*n_liq)("cell",c) / (*visc)("cell",c);
    }

    // communicate
    uw_rel_perm->ScatterMasterToGhosted();
  }

  return update_perm;
};


// -----------------------------------------------------------------------------
// Evaluate boundary conditions at the current time.
// -----------------------------------------------------------------------------
void Richards::UpdateBoundaryConditions_() {
  for (int n=0; n!=bc_markers_.size(); ++n) {
    bc_markers_[n] = Operators::MFD_BC_NULL;
    bc_values_[n] = 0.0;
  }

  Functions::BoundaryFunction::Iterator bc;
  for (bc=bc_pressure_->begin(); bc!=bc_pressure_->end(); ++bc) {
    int f = bc->first;
    bc_markers_[f] = Operators::MFD_BC_DIRICHLET;
    bc_values_[f] = bc->second;
  }

  if (!infiltrate_only_if_unfrozen_) {
    for (bc=bc_flux_->begin(); bc!=bc_flux_->end(); ++bc) {
      int f = bc->first;
      bc_markers_[f] = Operators::MFD_BC_FLUX;
      bc_values_[f] = bc->second;
    }
  } else {
    const Epetra_MultiVector& temp = *S_next_->GetFieldData("temperature")->ViewComponent("face");
    for (bc=bc_flux_->begin(); bc!=bc_flux_->end(); ++bc) {
      int f = bc->first;
      bc_markers_[f] = Operators::MFD_BC_FLUX;
      if (temp[0][f] > 273.15) {
        bc_values_[f] = bc->second;
      } else {
        bc_values_[f] = 0.;
      }
    }
  }

  // seepage face -- pressure <= p_atm, outward mass flux >= 0
  const Epetra_MultiVector pressure = *S_next_->GetFieldData(key_)->ViewComponent("face");
  const double& p_atm = *S_next_->GetScalarData("atmospheric_pressure");
  for (bc=bc_seepage_->begin(); bc!=bc_seepage_->end(); ++bc) {
    int f = bc->first;
    if (pressure[0][f] < p_atm) {
      bc_markers_[f] = Operators::MFD_BC_FLUX;
      bc_values_[f] = bc->second;
    } else {
      bc_markers_[f] = Operators::MFD_BC_DIRICHLET;
      bc_values_[f] = p_atm;
    }
  }

  // surface coupling
  if (coupled_to_surface_via_residual_) {
    // given the surface head, calculate a new pressure with surface head on
    // the surface faces
    Teuchos::RCP<const AmanziMesh::Mesh_MSTK> surface =
      Teuchos::rcp_static_cast<const AmanziMesh::Mesh_MSTK>(S_next_->GetMesh("surface"));
    const CompositeVector& ponded_depth = *S_next_->GetFieldData("ponded_depth");
    Teuchos::RCP<const CompositeVector> rho = S_next_->GetFieldData("mass_density_liquid");
    const CompositeVector& dens = *S_next_->GetFieldData("molar_density_liquid");
    const CompositeVector& surface_cell_volume =
        *S_next_->GetFieldData("surface_cell_volume");
    const double& p_atm = *S_next_->GetScalarData("atmospheric_pressure");
    Teuchos::RCP<const Epetra_Vector> gvec = S_next_->GetConstantVectorData("gravity");
    Teuchos::RCP<CompositeVector> source =
        S_next_->GetFieldData("overland_source_from_subsurface", name_);
    double dt = S_next_->time() - S_inter_->time();

    // -- set up a pressure vector containing current subsurface pressures
    // -- plus the effective pressure from the surface head
    CompositeVector pres(*S_next_->GetFieldData("pressure"));
    pres = *S_next_->GetFieldData("pressure");

    int ncells_surface = ponded_depth.size("cell");
    for (int c=0; c!=ncells_surface; ++c) {
      // -- get the surface cell's equivalent subsurface face and neighboring cell
      AmanziMesh::Entity_ID f =
        surface->entity_get_parent(AmanziMesh::CELL, c);
      AmanziMesh::Entity_ID_List cells;
      mesh_->face_get_cells(f, AmanziMesh::OWNED, &cells);
      ASSERT(cells.size() == 1);

      // -- assign lambda value
      //      pres("face",f) = std::max(
      //          -ponded_depth("cell",c)*(*rho)("cell",cells[0])*(*gvec)[2] + p_atm,
      //          p_atm);
      pres("face",f) = -ponded_depth("cell",c)*(*rho)("cell",cells[0])*(*gvec)[2] + p_atm;
    }

    // Calculate the fluxes that would result with this new head.
    UpdatePermeabilityData_(S_next_.ptr());
    S_next_->GetFieldEvaluator(key_)->HasFieldChanged(S_next_.ptr(), name_);
    S_next_->GetFieldEvaluator("mass_density_liquid")->HasFieldChanged(S_next_.ptr(), name_);

    Teuchos::RCP<const CompositeVector> rel_perm = S_next_->GetFieldData("numerical_rel_perm");
    // -- update the stiffness matrix
    matrix_->CreateMFDstiffnessMatrices(rel_perm.ptr());

    // -- derive the Darcy fluxes
    Teuchos::RCP<CompositeVector> darcy_flux = S_next_->GetFieldData("darcy_flux", name_);
    matrix_->DeriveFlux(pres, darcy_flux);

    // -- add in gravitational fluxes
    AddGravityFluxesToVector_(gvec, rel_perm, rho, darcy_flux);

    // Determine the BC
    for (int c=0; c!=ncells_surface; ++c) {
      // -- get the surface cell's equivalent subsurface face and neighboring cell
      AmanziMesh::Entity_ID f =
        surface->entity_get_parent(AmanziMesh::CELL, c);
      AmanziMesh::Entity_ID_List cells;
      mesh_->face_get_cells(f, AmanziMesh::OWNED, &cells);
      ASSERT(cells.size() == 1);

      // -- get the direction of the flux
      AmanziMesh::Entity_ID_List faces;
      std::vector<int> dirs;
      mesh_->cell_get_faces_and_dirs(cells[0], &faces, &dirs);
      int my_n = std::find(faces.begin(), faces.end(), f) - faces.begin();

      double q_out = (*darcy_flux)("face",f) * dirs[my_n];
      double Q_ss = ((*source)("cell",c)
                     - ponded_depth("cell",c) * surface_cell_volume("cell",c) / dt)
          * dens("cell",cells[0]); // residual mismatch, plus water available on surface

      if (q_out >= 0. || Q_ss >= 0.) {
        // Flux is outward, use the Dirichlet condition and the calculated flux
        bc_markers_[f] = Operators::MFD_BC_DIRICHLET;
        bc_values_[f] = pres("face",f);
        (*source)("cell",c) = q_out / dens("cell",cells[0]);
      } else {
        if (q_out < Q_ss) {
          // Flux wants to be inward and the subsurface
          // can accomodate all water on the surface.
          bc_markers_[f] = Operators::MFD_BC_FLUX;

          // note these are NOT the same necessarily!
          bc_values_[f] = Q_ss / mesh_->face_area(f); //surface_cell_volume("cell",c);
          (*source)("cell",c) = Q_ss / dens("cell",cells[0]);
        } else {
          // Flux wants to be inward, but the subsurface cannot handle the
          // entirety of the surface's water.
          bc_markers_[f] = Operators::MFD_BC_DIRICHLET;
          bc_values_[f] = pres("face",f);
          (*source)("cell",c) = q_out / dens("cell",cells[0]);
        }
      }


      //      if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_HIGH, true)) {

      if (bc_markers_[f] == Operators::MFD_BC_DIRICHLET) {
          //          AmanziGeometry::Point c_surf_point = surface->cell_centroid(c);
          //          AmanziGeometry::Point f_point = mesh_->face_centroid(f);
          //          AmanziGeometry::Point c_sub_point = mesh_->cell_centroid(cells[0]);

          Teuchos::OSTab tab = getOSTab();
          std::cout<< "SURFACE BC: c_surf: " << c << " f: " << f << " c_sub: " << cells[0] << std::endl;
          std::cout<< "   sizes:  cell_area: " << surface_cell_volume("cell",c) << "  face area: " << mesh_->face_area(f) << std::endl;
          //          std::cout<< "    c_surf: " << c_surf_point << "    f: " << f_point << "    c_sub: " << c_sub_point << std::endl;
          std::cout<< "            p = " << pres("cell",cells[0]) << " p_eff = " << pres("face",f) << std::endl;
          std::cout<< "            h = " << ponded_depth("cell",c) << ", q_out = " << q_out << ", Q_ss = " << Q_ss << std::endl;
          if (bc_markers_[f] == Operators::MFD_BC_FLUX) {
            std::cout<< "  RESULT:  flux, " << bc_values_[f] << std::endl;
          } else {
            std::cout<< "  RESULT:  dirichlet, " << bc_values_[f] << std::endl;
          }
        }
        // }
    }
  }

};


// -----------------------------------------------------------------------------
// Add a boundary marker to owned faces.
// -----------------------------------------------------------------------------
void
Richards::ApplyBoundaryConditions_(const Teuchos::RCP<CompositeVector>& pres) {
  int nfaces = pres->size("face");
  for (int f=0; f!=nfaces; ++f) {
    if (bc_markers_[f] == Operators::MFD_BC_DIRICHLET) {
      (*pres)("face",f) = bc_values_[f];
    }
  }
};


// -----------------------------------------------------------------------------
// Experimental approach -- calling this indicates that the time
// integration scheme is changing the value of the solution in
// state.
// -----------------------------------------------------------------------------
void Richards::changed_solution() {
  solution_evaluator_->SetFieldAsChanged();
  S_next_->GetFieldData(key_)->ScatterMasterToGhosted("face");
};


bool Richards::modify_predictor(double h, const Teuchos::RCP<TreeVector>& u) {
  if (modify_predictor_with_consistent_faces_) {
    // derive consistent constraints

    Teuchos::RCP<CompositeVector> uw_rel_perm = S_next_->GetFieldData("numerical_rel_perm", name_);

    // VerboseObject stuff.
    Teuchos::OSTab tab = getOSTab();

    // update the rel perm according to the scheme of choice
    UpdatePermeabilityData_(S_next_.ptr());

    // update boundary conditions
    bc_pressure_->Compute(S_next_->time());
    bc_flux_->Compute(S_next_->time());
    UpdateBoundaryConditions_();

    Teuchos::RCP<const CompositeVector> rho =
      S_next_->GetFieldData("mass_density_liquid");
    Teuchos::RCP<const Epetra_Vector> gvec =
      S_next_->GetConstantVectorData("gravity");

    // Update the preconditioner with darcy and gravity fluxes
    preconditioner_->CreateMFDstiffnessMatrices(uw_rel_perm.ptr());
    preconditioner_->CreateMFDrhsVectors();
    AddGravityFluxes_(gvec, uw_rel_perm, rho, preconditioner_);

    // skip accumulation terms, they're not needed

    // Assemble and precompute the Schur complement for inversion.
    preconditioner_->ApplyBoundaryConditions(bc_markers_, bc_values_);
    preconditioner_->AssembleGlobalMatrices();

    // derive the consistent faces, involves a solve
    preconditioner_->UpdateConsistentFaceConstraints(u->data());
    return true;
  }

  return PKPhysicalBDFBase::modify_predictor(h, u);
}

} // namespace
} // namespace
