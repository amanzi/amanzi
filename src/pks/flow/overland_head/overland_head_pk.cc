/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -----------------------------------------------------------------------------
This is the overland flow component of ATS.
License: BSD
Author: Ethan Coon (ecoon@lanl.gov)
----------------------------------------------------------------------------- */

#include "EpetraExt_MultiVectorOut.h"
#include "Epetra_MultiVector.h"

#include "bdf1_time_integrator.hh"
#include "flow_bc_factory.hh"
#include "Mesh.hh"
#include "Point.hh"

#include "composite_vector_function.hh"
#include "composite_vector_function_factory.hh"
#include "independent_variable_field_evaluator.hh"

#include "matrix_mfd_tpfa.hh"
#include "upwinding.hh"
#include "upwind_potential_difference.hh"
#include "elevation_evaluator.hh"
#include "meshed_elevation_evaluator.hh"
#include "standalone_elevation_evaluator.hh"
#include "overland_conductivity_evaluator.hh"
#include "overland_conductivity_model.hh"
#include "overland_head_water_content_evaluator.hh"
#include "height_model.hh"
#include "height_evaluator.hh"
#include "overland_source_from_subsurface_flux_evaluator.hh"

#include "overland_head.hh"

namespace Amanzi {
namespace Flow {

#define DEBUG_FLAG 1

RegisteredPKFactory<OverlandHeadFlow> OverlandHeadFlow::reg_("overland flow, head basis");

// -------------------------------------------------------------
// Constructor
// -------------------------------------------------------------
void OverlandHeadFlow::setup(const Teuchos::Ptr<State>& S) {
  // set up the meshes
  if (!S->HasMesh("surface")) {
    Teuchos::RCP<const AmanziMesh::Mesh> domain = S->GetMesh();
    ASSERT(domain->space_dimension() == 2);
    standalone_mode_ = true;
    S->RegisterMesh("surface", domain);
  } else {
    standalone_mode_ = false;
  }

  PKPhysicalBDFBase::setup(S);
  setLinePrefix("overland");  // make the prefix fit in the available space.
  SetupOverlandFlow_(S);
  SetupPhysicalEvaluators_(S);
}


void OverlandHeadFlow::SetupOverlandFlow_(const Teuchos::Ptr<State>& S) {
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

  // -- owned secondary variables, no evaluator used
  S->RequireField("surface_flux", name_)->SetMesh(S->GetMesh("surface"))
                ->SetGhosted()->SetComponent("face", AmanziMesh::FACE, 1);
  S->RequireField("surface_velocity", name_)->SetMesh(S->GetMesh("surface"))
                ->SetGhosted()->SetComponent("cell", AmanziMesh::CELL, 3);

  // -- cell volume and evaluator
  S->RequireFieldEvaluator("surface_cell_volume");
  S->RequireGravity();
  S->RequireScalar("atmospheric_pressure");

  // -- coupling to subsurface
  coupled_to_subsurface_via_flux_ =
      plist_.get<bool>("coupled to subsurface via flux", false);
  coupled_to_subsurface_via_head_ =
      plist_.get<bool>("coupled to subsurface via head", false);
  coupled_to_subsurface_via_full_ =
      plist_.get<bool>("coupled to subsurface via full coupler", false);
  ASSERT(!(coupled_to_subsurface_via_flux_ && coupled_to_subsurface_via_head_));
  ASSERT(!(coupled_to_subsurface_via_full_ && coupled_to_subsurface_via_head_));
  ASSERT(!(coupled_to_subsurface_via_flux_ && coupled_to_subsurface_via_full_));

  if (coupled_to_subsurface_via_flux_) {
    // -- coupling term in PC, filled in by subsurface PK.
    S->RequireField("doverland_source_from_subsurface_dsurface_pressure")
        ->SetMesh(mesh_)->AddComponent("cell", AmanziMesh::CELL, 1);

    // -- source term from subsurface, owned by me but filled in by the coupler
    S->RequireField("overland_source_from_subsurface", name_)
        ->SetMesh(mesh_)->SetComponent("cell", AmanziMesh::CELL, 1);
  }

  if (coupled_to_subsurface_via_head_) {
    // -- coupling term in PC, filled in by subsurface PK.
    S->RequireField("doverland_source_from_subsurface_dsurface_pressure")
        ->SetMesh(mesh_)->AddComponent("cell", AmanziMesh::CELL, 1);
  }

  if (coupled_to_subsurface_via_head_ || coupled_to_subsurface_via_full_) {
    // -- source term from subsurface, filled in by evaluator,
    //    which picks the fluxes from "darcy_flux" field.
    S->RequireField("overland_source_from_subsurface")
        ->SetMesh(mesh_)->AddComponent("cell", AmanziMesh::CELL, 1);

    Teuchos::ParameterList source_plist =
        plist_.sublist("source from subsurface evaluator");
    Teuchos::RCP<Relations::OverlandSourceFromSubsurfaceFluxEvaluator>
        source_evaluator = Teuchos::rcp(
            new Relations::OverlandSourceFromSubsurfaceFluxEvaluator(source_plist));
    S->SetFieldEvaluator("overland_source_from_subsurface", source_evaluator);
  }

  // Create the boundary condition data structures.
  Teuchos::ParameterList bc_plist = plist_.sublist("boundary conditions", true);
  FlowBCFactory bc_factory(mesh_, bc_plist);
  bc_pressure_ = bc_factory.CreatePressure();
  bc_head_ = bc_factory.CreateHead();
  bc_zero_gradient_ = bc_factory.CreateZeroGradient();
  bc_flux_ = bc_factory.CreateMassFlux();

  // Create the upwinding method.
  S->RequireField("upwind_overland_conductivity", name_)->SetMesh(mesh_)
    ->SetGhosted()->SetComponents(names2, locations2, num_dofs2);
  S->GetField("upwind_overland_conductivity",name_)->set_io_vis(false);

  upwinding_ = Teuchos::rcp(new Operators::UpwindPotentialDifference(name_,
          "overland_conductivity", "upwind_overland_conductivity",
          "pres_elev", "ponded_depth"));

  // operator for the diffusion terms
  Teuchos::ParameterList mfd_plist = plist_.sublist("Diffusion");
  matrix_ = Teuchos::rcp(new Operators::MatrixMFD(mfd_plist, mesh_));
  symmetric_ = false;
  matrix_->SetSymmetryProperty(symmetric_);
  matrix_->SymbolicAssembleGlobalMatrices();
  matrix_->CreateMFDmassMatrices(Teuchos::null);


  Teuchos::ParameterList mfd_pc_plist = plist_.sublist("Diffusion PC");
  Teuchos::RCP<Operators::Matrix> precon;
  if (coupled_to_subsurface_via_full_) {
    precon = Teuchos::rcp(new Operators::MatrixMFD_TPFA(mfd_pc_plist, mesh_));
  } else {
    precon = Teuchos::rcp(new Operators::MatrixMFD(mfd_pc_plist, mesh_));
  }
  set_preconditioner(precon);
  assemble_preconditioner_ = plist_.get<bool>("assemble preconditioner", true);
  modify_predictor_with_consistent_faces_ =
    plist_.get<bool>("modify predictor with consistent faces", false);

  // how often to update the fluxes?
  std::string updatestring = plist_.get<std::string>("update flux mode", "iteration");
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
};


// -------------------------------------------------------------
// Create the physical evaluators for water content, water
// retention, rel perm, etc, that are specific to Richards.
// -------------------------------------------------------------
void OverlandHeadFlow::SetupPhysicalEvaluators_(const Teuchos::Ptr<State>& S) {
  std::vector<AmanziMesh::Entity_kind> locations2(2);
  std::vector<std::string> names2(2);
  std::vector<int> num_dofs2(2,1);
  locations2[0] = AmanziMesh::CELL;
  locations2[1] = AmanziMesh::FACE;
  names2[0] = "cell";
  names2[1] = "face";

  // -- evaluator for surface geometry.
  S->RequireField("elevation")->SetMesh(S->GetMesh("surface"))->SetGhosted()
                ->AddComponents(names2, locations2, num_dofs2);
  S->RequireField("slope_magnitude")->SetMesh(S->GetMesh("surface"))
                ->SetGhosted()->AddComponent("cell", AmanziMesh::CELL, 1);
  S->RequireField("pres_elev")->SetMesh(S->GetMesh("surface"))->SetGhosted()
                ->AddComponents(names2, locations2, num_dofs2);

  Teuchos::RCP<FlowRelations::ElevationEvaluator> elev_evaluator;
  if (standalone_mode_) {
    ASSERT(plist_.isSublist("elevation evaluator"));
    Teuchos::ParameterList elev_plist = plist_.sublist("elevation evaluator");
    elev_evaluator = Teuchos::rcp(new FlowRelations::StandaloneElevationEvaluator(elev_plist));
  } else {
    Teuchos::ParameterList elev_plist = plist_.sublist("elevation evaluator");
    elev_evaluator = Teuchos::rcp(new FlowRelations::MeshedElevationEvaluator(elev_plist));
  }
  S->SetFieldEvaluator("elevation", elev_evaluator);
  S->SetFieldEvaluator("slope_magnitude", elev_evaluator);
  S->SetFieldEvaluator("pres_elev", elev_evaluator);

  // -- conductivity evaluator
  S->RequireField("overland_conductivity")->SetMesh(mesh_)->SetGhosted()
    ->AddComponent("cell", AmanziMesh::CELL, 1);
  ASSERT(plist_.isSublist("overland conductivity evaluator"));
  Teuchos::ParameterList cond_plist = plist_.sublist("overland conductivity evaluator");
  Teuchos::RCP<FlowRelations::OverlandConductivityEvaluator> cond_evaluator =
      Teuchos::rcp(new FlowRelations::OverlandConductivityEvaluator(cond_plist));
  S->SetFieldEvaluator("overland_conductivity", cond_evaluator);
  cond_model_ = cond_evaluator->get_Model();

  // -- source term evaluator
  if (plist_.isSublist("source evaluator")) {
    Teuchos::ParameterList source_plist = plist_.sublist("source evaluator");
    source_plist.set("evaluator name", "overland_source");
    is_source_term_ = true;
    S->RequireField("overland_source")->SetMesh(mesh_)
        ->SetGhosted()->AddComponent("cell", AmanziMesh::CELL, 1);
    Teuchos::RCP<FieldEvaluator> source_evaluator =
        Teuchos::rcp(new IndependentVariableFieldEvaluator(source_plist));
    S->SetFieldEvaluator("overland_source", source_evaluator);
  }

  // -- water content
  S->RequireField("surface_water_content")->SetMesh(mesh_)->SetGhosted()
      ->AddComponent("cell", AmanziMesh::CELL, 1);
  Teuchos::ParameterList wc_plist = plist_.sublist("overland water content");
  Teuchos::RCP<FlowRelations::OverlandHeadWaterContentEvaluator> wc_evaluator =
      Teuchos::rcp(new FlowRelations::OverlandHeadWaterContentEvaluator(wc_plist));
  S->SetFieldEvaluator("surface_water_content", wc_evaluator);

  // -- ponded depth
  S->RequireField("ponded_depth")->SetMesh(mesh_)->SetGhosted()
                ->AddComponents(names2, locations2, num_dofs2);
  Teuchos::ParameterList height_plist = plist_.sublist("ponded depth");
  Teuchos::RCP<FlowRelations::HeightEvaluator> height_evaluator =
      Teuchos::rcp(new FlowRelations::HeightEvaluator(height_plist));
  S->SetFieldEvaluator("ponded_depth", height_evaluator);
}


// -------------------------------------------------------------
// Initialize PK
// -------------------------------------------------------------
void OverlandHeadFlow::initialize(const Teuchos::Ptr<State>& S) {
  // initial condition is tricky
  // -- set the cell initial condition if it is taken from the subsurface
  if (!plist_.isSublist("initial condition")) {
    std::stringstream messagestream;
    messagestream << name_ << " has no initial condition parameter list.";
    Errors::Message message(messagestream.str());
    Exceptions::amanzi_throw(message);
  }

  Teuchos::ParameterList ic_plist = plist_.sublist("initial condition");
  if (ic_plist.get<bool>("initialize surface head from subsurface",false)) {
    Epetra_MultiVector& head = *S->GetFieldData(key_, name_)
        ->ViewComponent("cell",false);
    const Epetra_MultiVector& pres = *S->GetFieldData("pressure")
        ->ViewComponent("face",false);

    int ncells_surface = mesh_->num_entities(AmanziMesh::CELL,AmanziMesh::OWNED);
    for (int c=0; c!=ncells_surface; ++c) {
      // -- get the surface cell's equivalent subsurface face and neighboring cell
      AmanziMesh::Entity_ID f =
          mesh_->entity_get_parent(AmanziMesh::CELL, c);
      head[0][c] = pres[0][f];
    }
    S->GetField(key_,name_)->set_initialized();
  }
  S->GetFieldData(key_,name_)->ViewComponent("face",false)->PutScalar(0.);

  // Initialize BDF stuff and physical domain stuff.
  PKPhysicalBDFBase::initialize(S);

  // Initialize boundary conditions.
  int nfaces = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::USED);
  bc_markers_.resize(nfaces, Operators::MATRIX_BC_NULL);
  bc_values_.resize(nfaces, 0.0);

  bc_pressure_->Compute(S->time());
  bc_head_->Compute(S->time());
  bc_zero_gradient_->Compute(S->time());
  bc_flux_->Compute(S->time());
  //  UpdateBoundaryConditions_(S);

  // Set extra fields as initialized -- these don't currently have evaluators.
  S->GetFieldData("upwind_overland_conductivity",name_)->PutScalar(1.0);
  S->GetField("upwind_overland_conductivity",name_)->set_initialized();
  S->GetField("surface_flux", name_)->set_initialized();
  S->GetField("surface_velocity", name_)->set_initialized();

  // initialize coupling terms
  if (coupled_to_subsurface_via_flux_) {
    S->GetFieldData("overland_source_from_subsurface", name_)
        ->PutScalar(0.);
    S->GetField("overland_source_from_subsurface", name_)
        ->set_initialized();
  }
};


// -----------------------------------------------------------------------------
// Update any secondary (dependent) variables given a solution.
//
//   After a timestep is evaluated (or at ICs), there is no way of knowing if
//   secondary variables have been updated to be consistent with the new
//   solution.
// -----------------------------------------------------------------------------
void OverlandHeadFlow::commit_state(double dt, const Teuchos::RCP<State>& S) {
  // update boundary conditions
  bc_pressure_->Compute(S->time());
  bc_head_->Compute(S->time());
  bc_flux_->Compute(S->time());
  UpdateBoundaryConditions_(S.ptr());

  // Update flux if rel perm or h + Z has changed.
  bool update = UpdatePermeabilityData_(S.ptr());
  update |= S->GetFieldEvaluator("pres_elev")->HasFieldChanged(S.ptr(), name_);

  if (update_flux_ == UPDATE_FLUX_TIMESTEP ||
      (update_flux_ == UPDATE_FLUX_ITERATION && update)) {
    // update the stiffness matrix with the new rel perm
    Teuchos::RCP<const CompositeVector> conductivity =
        S->GetFieldData("upwind_overland_conductivity");
    matrix_->CreateMFDstiffnessMatrices(conductivity.ptr());

    // derive the fluxes
    Teuchos::RCP<const CompositeVector> potential = S->GetFieldData("pres_elev");
    Teuchos::RCP<CompositeVector> flux = S->GetFieldData("surface_flux", name_);
    potential->ScatterMasterToGhosted();
    matrix_->DeriveFlux(*potential, flux.ptr());
  }
};


// -----------------------------------------------------------------------------
// Update diagnostics -- used prior to vis.
// -----------------------------------------------------------------------------
void OverlandHeadFlow::calculate_diagnostics(const Teuchos::RCP<State>& S) {
  // update the cell velocities
  if (update_flux_ == UPDATE_FLUX_VIS) {
    Teuchos::RCP<CompositeVector> flux = S->GetFieldData("surface_flux",name_);
    Teuchos::RCP<const CompositeVector> conductivity =
        S->GetFieldData("upwind_overland_conductivity");
    matrix_->CreateMFDstiffnessMatrices(conductivity.ptr());

    // derive the fluxes
    S->GetFieldEvaluator("pres_elev")->HasFieldChanged(S.ptr(), name_);
    Teuchos::RCP<const CompositeVector> potential = S->GetFieldData("pres_elev");
    potential->ScatterMasterToGhosted();
    matrix_->DeriveFlux(*potential, flux.ptr());
  }

  if (update_flux_ != UPDATE_FLUX_NEVER) {
    Teuchos::RCP<const CompositeVector> flux = S->GetFieldData("surface_flux");
    Teuchos::RCP<CompositeVector> velocity = S->GetFieldData("surface_velocity", name_);
    matrix_->DeriveCellVelocity(*flux, velocity.ptr());

    Teuchos::RCP<const CompositeVector> depth = S->GetFieldData("ponded_depth");
    const Epetra_MultiVector& depth_c = *depth->ViewComponent("cell",false);
    Epetra_MultiVector& vel_c = *velocity->ViewComponent("cell",false);

    int ncells = velocity->size("cell");
    for (int c=0; c!=ncells; ++c) {
      vel_c[0][c] /= std::max( depth_c[0][c] , 1e-7);
      vel_c[1][c] /= std::max( depth_c[0][c] , 1e-7);
    }
  }
};


// -----------------------------------------------------------------------------
// Use the physical rel perm (on cells) to update a work vector for rel perm.
//
//   This deals with upwinding, etc.
// -----------------------------------------------------------------------------
bool OverlandHeadFlow::UpdatePermeabilityData_(const Teuchos::Ptr<State>& S) {
  bool update_perm = S->GetFieldEvaluator("overland_conductivity")
      ->HasFieldChanged(S, name_);

  update_perm |= S->GetFieldEvaluator("ponded_depth")->HasFieldChanged(S, name_);
  update_perm |= S->GetFieldEvaluator("pres_elev")->HasFieldChanged(S, name_);
  update_perm |= perm_update_required_;

  if (update_perm) {
    // Update the perm only if needed.
    perm_update_required_ = false;

    // This needs fixed to use the model, not assume a model.
    // Then it needs to be fixed to use a smart evaluator which picks
    // vals from cells and faces.
    // Then it needs to be fixed to work on a CompositeVector whose
    // components are cells and boundary faces. --etc
    const Epetra_MultiVector& depth = *S->GetFieldData("ponded_depth")
      ->ViewComponent("face", false);
    const Epetra_MultiVector& slope = *S->GetFieldData("slope_magnitude")
      ->ViewComponent("cell", false);
    const Epetra_MultiVector& coef = *S->GetFieldData("manning_coefficient")
      ->ViewComponent("cell", false);

    Teuchos::RCP<CompositeVector> upwind_conductivity =
        S->GetFieldData("upwind_overland_conductivity", name_);
    Epetra_MultiVector& upwind_conductivity_f =
        *upwind_conductivity->ViewComponent("face",false);
    Epetra_MultiVector& upwind_conductivity_c =
        *upwind_conductivity->ViewComponent("cell",false);

    // initialize the face coefficients
    upwind_conductivity_f.PutScalar(0.0);
    upwind_conductivity_c.PutScalar(1.0);

    // First do the boundary.
    AmanziMesh::Entity_ID_List cells;
    int nfaces = upwind_conductivity_f.MyLength();
    for (int f=0; f!=nfaces; ++f) {
      if (bc_markers_[f] != Operators::MATRIX_BC_NULL) {
        mesh_->face_get_cells(f, AmanziMesh::USED, &cells);
        ASSERT(cells.size() == 1);
        int c = cells[0];
        upwind_conductivity_f[0][f] =
          cond_model_->Conductivity(depth[0][f], slope[0][c], coef[0][c]);
      }
    }

    // Patch up zero-gradient case, which should not upwind.
    const Epetra_MultiVector& conductivity_c =
        *S->GetFieldData("overland_conductivity")->ViewComponent("cell",false);
    for (Functions::BoundaryFunction::Iterator bc=bc_zero_gradient_->begin();
         bc!=bc_zero_gradient_->end(); ++bc) {
      int f = bc->first;
      mesh_->face_get_cells(f, AmanziMesh::USED, &cells);
      ASSERT(cells.size() == 1);
      int c = cells[0];
      upwind_conductivity_f[0][f] = conductivity_c[0][c];
    }

    // Then upwind.  This overwrites the boundary if upwinding says so.
    upwinding_->Update(S);

    // Scale cells by n
    const Epetra_MultiVector& n_liq = *S->GetFieldData("surface_molar_density_liquid")
        ->ViewComponent("cell",false);
    int ncells = upwind_conductivity_c.MyLength();
    for (int c=0; c!=ncells; ++c) {
      upwind_conductivity_c[0][c] *= n_liq[0][c];
    }

    // Communicate.  This could be done later, but i'm not exactly sure where, so
    // we'll do it here.
    upwind_conductivity->ScatterMasterToGhosted();
  }

  return update_perm;
}


// -----------------------------------------------------------------------------
// Evaluate boundary conditions at the current time.
// -----------------------------------------------------------------------------
void OverlandHeadFlow::UpdateBoundaryConditions_(const Teuchos::Ptr<State>& S) {
  AmanziMesh::Entity_ID_List cells;

  const Epetra_MultiVector& elevation = *S->GetFieldData("elevation")
      ->ViewComponent("face",false);

  // initialize all as null
  for (int n=0; n!=bc_markers_.size(); ++n) {
    bc_markers_[n] = Operators::MATRIX_BC_NULL;
    bc_values_[n] = 0.0;
  }

  // Head BCs are standard Dirichlet, plus elevation
  for (Functions::BoundaryFunction::Iterator bc=bc_head_->begin();
       bc!=bc_head_->end(); ++bc) {
    int f = bc->first;
    bc_markers_[f] = Operators::MATRIX_BC_DIRICHLET;
    bc_values_[f] = bc->second + elevation[0][f];
  }

  // pressure BCs require calculating head(pressure)
  const double& p_atm = *S->GetScalarData("atmospheric_pressure");
  const Epetra_Vector& gravity = *S->GetConstantVectorData("gravity");
  double gz = -gravity[2];
  S->GetFieldEvaluator("surface_mass_density_liquid")->HasFieldChanged(S.ptr(), name_);
  const Epetra_MultiVector& rho = *S->GetFieldData("surface_mass_density_liquid")
      ->ViewComponent("cell",false);

  for (Functions::BoundaryFunction::Iterator bc=bc_pressure_->begin();
       bc!=bc_pressure_->end(); ++bc) {
    int f = bc->first;
    mesh_->face_get_cells(f, AmanziMesh::USED, &cells);
    ASSERT( cells.size()==1 ) ;
    double height = height_model_->Height(bc->second, rho[0][cells[0]], p_atm, gz);

    bc_markers_[f] = Operators::MATRIX_BC_DIRICHLET;
    bc_values_[f] = height + elevation[0][f];
  }

  // Standard Neumann data for flux
  for (Functions::BoundaryFunction::Iterator bc=bc_flux_->begin();
       bc!=bc_flux_->end(); ++bc) {
    int f = bc->first;
    bc_markers_[f] = Operators::MATRIX_BC_FLUX;
    bc_values_[f] = bc->second;
  }

  // zero gradient: grad h = 0 implies that q = -k grad z
  for (Functions::BoundaryFunction::Iterator bc=bc_zero_gradient_->begin();
       bc!=bc_zero_gradient_->end(); ++bc) {
    int f = bc->first;
    bc_markers_[f] = Operators::MATRIX_BC_FLUX;
    bc_values_[f] = 0.;
  }

}


// -----------------------------------------------------------------------------
// Evaluate boundary conditions at the current time without elevation.
// -----------------------------------------------------------------------------
void OverlandHeadFlow::UpdateBoundaryConditionsMarkers_(const Teuchos::Ptr<State>& S) {
  AmanziMesh::Entity_ID_List cells;

  // initialize all as null
  for (int n=0; n!=bc_markers_.size(); ++n) {
    bc_markers_[n] = Operators::MATRIX_BC_NULL;
  }

  // Head BCs are standard Dirichlet, no elevation
  for (Functions::BoundaryFunction::Iterator bc=bc_head_->begin();
       bc!=bc_head_->end(); ++bc) {
    int f = bc->first;
    bc_markers_[f] = Operators::MATRIX_BC_DIRICHLET;
  }

  // pressure BCs require calculating head(pressure)
  for (Functions::BoundaryFunction::Iterator bc=bc_pressure_->begin();
       bc!=bc_pressure_->end(); ++bc) {
    int f = bc->first;
    bc_markers_[f] = Operators::MATRIX_BC_DIRICHLET;
  }

  // Standard Neumann data for flux
  for (Functions::BoundaryFunction::Iterator bc=bc_flux_->begin();
       bc!=bc_flux_->end(); ++bc) {
    int f = bc->first;
    bc_markers_[f] = Operators::MATRIX_BC_FLUX;
  }

  // zero gradient: grad h = 0 implies that q = -k grad z
  for (Functions::BoundaryFunction::Iterator bc=bc_zero_gradient_->begin();
       bc!=bc_zero_gradient_->end(); ++bc) {
    int f = bc->first;
    bc_markers_[f] = Operators::MATRIX_BC_FLUX;
  }
}


void OverlandHeadFlow::FixBCsForOperator_(const Teuchos::Ptr<State>& S) {
  // If the rel perm is 0, the face value drops out and is unconstrained.
  // Therefore we set it to Dirichlet to eliminate it from the system.
  double eps = 1.e-30;
  const Epetra_MultiVector& elevation = *S->GetFieldData("elevation")
      ->ViewComponent("face",false);
  Teuchos::RCP<CompositeVector> relperm =
      S->GetFieldData("upwind_overland_conductivity", name_);

  for (int f=0; f!=relperm->size("face"); ++f) {
    if ((*relperm)("face",f) < eps) {
      if (bc_markers_[f] == Operators::MATRIX_BC_FLUX) {
        bc_markers_[f] = Operators::MATRIX_BC_DIRICHLET;
        bc_values_[f] =  elevation[0][f];
      } else if (bc_markers_[f] == Operators::MATRIX_BC_NULL) {
        bc_markers_[f] = Operators::MATRIX_BC_DIRICHLET;
        bc_values_[f] =  elevation[0][f];
      }
    }
  }

  // Now we can safely calculate q = -k grad z for zero-gradient problems
  AmanziMesh::Entity_ID_List cells;
  AmanziMesh::Entity_ID_List faces;
  std::vector<double> dp;
  std::vector<int> dirs;

  const Epetra_MultiVector& elevation_c = *S->GetFieldData("elevation")
      ->ViewComponent("cell",false);

  std::vector<Teuchos::SerialDenseMatrix<int, double> >& Aff_cells =
      matrix_->Aff_cells();

  for (Functions::BoundaryFunction::Iterator bc=bc_zero_gradient_->begin();
       bc!=bc_zero_gradient_->end(); ++bc) {
    int f = bc->first;
    if (bc_markers_[f] == Operators::MATRIX_BC_FLUX) {
      mesh_->face_get_cells(f, AmanziMesh::USED, &cells);
      ASSERT(cells.size() == 1);
      AmanziMesh::Entity_ID c = cells[0];
      mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);

      dp.resize(faces.size());
      for (int n=0; n!=faces.size(); ++n) {
        dp[n] = elevation[0][faces[n]] - elevation_c[0][c];
      }
      int my_n = std::find(faces.begin(), faces.end(), f) - faces.begin();
      ASSERT(my_n !=faces.size());

      bc_values_[f] = 0.;
      for (int m=0; m!=faces.size(); ++m) {
        bc_values_[f] -= Aff_cells[c](my_n,m) * dp[m];
      }
    }
  }

};


void OverlandHeadFlow::FixBCsForPrecon_(const Teuchos::Ptr<State>& S) {
  // Attempt of a hack to deal with zero rel perm
  double eps = 1.e-30;
  Teuchos::RCP<CompositeVector> relperm =
      S->GetFieldData("upwind_overland_conductivity", name_);
  for (int f=0; f!=relperm->size("face"); ++f) {
    if ((*relperm)("face",f) < eps) {
      if (bc_markers_[f] == Operators::MATRIX_BC_FLUX) {
        bc_markers_[f] = Operators::MATRIX_BC_DIRICHLET;
      } else if (bc_markers_[f] == Operators::MATRIX_BC_NULL) {
        bc_markers_[f] = Operators::MATRIX_BC_DIRICHLET;
      }
    }
  }
};

void OverlandHeadFlow::FixBCsForConsistentFaces_(const Teuchos::Ptr<State>& S) {
  // If the rel perm is 0, the face value drops out and is unconstrained.
  // Therefore we set it to Dirichlet to eliminate it from the system.
  double eps = 1.e-30;
  const Epetra_MultiVector& elevation = *S->GetFieldData("elevation")
      ->ViewComponent("face",false);
  Teuchos::RCP<CompositeVector> relperm =
      S->GetFieldData("upwind_overland_conductivity", name_);

  for (int f=0; f!=relperm->size("face"); ++f) {
    if ((*relperm)("face",f) < eps) {
      if (bc_markers_[f] == Operators::MATRIX_BC_FLUX) {
        bc_markers_[f] = Operators::MATRIX_BC_DIRICHLET;
        bc_values_[f] =  elevation[0][f];
      } else if (bc_markers_[f] == Operators::MATRIX_BC_NULL) {
        bc_markers_[f] = Operators::MATRIX_BC_DIRICHLET;
        bc_values_[f] =  elevation[0][f];
      }
    }
  }

  // Now we can safely calculate q = -k grad z for zero-gradient problems
  AmanziMesh::Entity_ID_List cells;
  AmanziMesh::Entity_ID_List faces;
  std::vector<double> dp;
  std::vector<int> dirs;

  const Epetra_MultiVector& elevation_c = *S->GetFieldData("elevation")
      ->ViewComponent("cell",false);

  std::vector<Teuchos::SerialDenseMatrix<int, double> >& Aff_cells =
      mfd_preconditioner_->Aff_cells();

  for (Functions::BoundaryFunction::Iterator bc=bc_zero_gradient_->begin();
       bc!=bc_zero_gradient_->end(); ++bc) {
    int f = bc->first;
    if (bc_markers_[f] == Operators::MATRIX_BC_FLUX) {
      mesh_->face_get_cells(f, AmanziMesh::USED, &cells);
      ASSERT(cells.size() == 1);
      AmanziMesh::Entity_ID c = cells[0];
      mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);

      dp.resize(faces.size());
      for (int n=0; n!=faces.size(); ++n) {
        dp[n] = elevation[0][faces[n]] - elevation_c[0][c];
      }
      int my_n = std::find(faces.begin(), faces.end(), f) - faces.begin();
      ASSERT(my_n !=faces.size());

      bc_values_[f] = 0.;
      for (int m=0; m!=faces.size(); ++m) {
        bc_values_[f] -= Aff_cells[c](my_n,m) * dp[m];
      }
    }
  }
};


/* ******************************************************************
 * Add a boundary marker to owned faces.
 ****************************************************************** */
void OverlandHeadFlow::ApplyBoundaryConditions_(const Teuchos::RCP<State>& S,
        const Teuchos::RCP<CompositeVector>& pres) {
  int nfaces = pres->size("face",true);
  for (int f=0; f!=nfaces; ++f) {
    if (bc_markers_[f] == Operators::MATRIX_BC_DIRICHLET) {
      (*pres)("face",f) = bc_values_[f];
    }
  }
};


// -----------------------------------------------------------------------------
// Experimental approach -- calling this indicates that the time
// integration scheme is changing the value of the solution in
// state.
// -----------------------------------------------------------------------------
void OverlandHeadFlow::changed_solution() {
  solution_evaluator_->SetFieldAsChanged();
  // communicate both faces and cells -- faces for flux and cells for rel perm
  // upwinding overlap
  S_next_->GetFieldData(key_)->ScatterMasterToGhosted();
};


bool OverlandHeadFlow::modify_predictor(double h, Teuchos::RCP<TreeVector> u) {
  bool changed = false;
  // const double& patm = *S_next_->GetScalarData("atmospheric_pressure");
  // const Epetra_MultiVector& u_prev_c =
  //   *S_->GetFieldData(key_)->ViewComponent("cell",false);
  // Epetra_MultiVector& u_c = *u->data()->ViewComponent("cell",false);

  // // Damp the spurt of water
  // int ncells = u_c.MyLength();
  // for (int c=0; c!=ncells; ++c) {
  //   if ((u_prev_c[0][c] < patm) &&
  //       (u_c[0][c] > patm)) {
  //     u_c[0][c] = patm - 1.;
  //     changed = true;
  //   }
  // }

  if (modify_predictor_with_consistent_faces_) {
    CalculateConsistentFaces(u->data().ptr());
    return true;
  }
  if (changed) return true;

  return PKPhysicalBDFBase::modify_predictor(h, u);
};


void OverlandHeadFlow::CalculateConsistentFaces(const Teuchos::Ptr<CompositeVector>& u) {
  // VerboseObject stuff.
  Teuchos::OSTab tab = getOSTab();

  // update boundary conditions
  bc_pressure_->Compute(S_next_->time());
  bc_head_->Compute(S_next_->time());
  bc_flux_->Compute(S_next_->time());
  UpdateBoundaryConditions_(S_next_.ptr());

  // update the rel perm according to the scheme of choice
  UpdatePermeabilityData_(S_next_.ptr());

  // update the stiffness matrix
  Teuchos::RCP<const CompositeVector> cond =
    S_next_->GetFieldData("upwind_overland_conductivity", name_);
  mfd_preconditioner_->CreateMFDstiffnessMatrices(cond.ptr());

  // Patch up BCs in the case of zero conductivity
  FixBCsForConsistentFaces_(S_next_.ptr());

  // Grab needed data.
  S_next_->GetFieldEvaluator("pres_elev")->HasFieldChanged(S_next_.ptr(), name_);
  Teuchos::RCP<CompositeVector> pres_elev = S_next_->GetFieldData("pres_elev","pres_elev");

#if DEBUG_FLAG
  *out_ << "Consistent faces:" << std::endl;
  *out_ << "  u_c = " << (*u)("cell",3) << ",  pres_e_c = " << (*pres_elev)("cell",3) << std::endl;
  *out_ << "  u_c = " << (*u)("cell",4) << ",  pres_e_c = " << (*pres_elev)("cell",4) << std::endl;

  AmanziMesh::Entity_ID f = mesh_->entity_get_parent(AmanziMesh::CELL, 3);
  AmanziMesh::Entity_ID_List cells;
  S_next_->GetMesh()->face_get_cells(f, AmanziMesh::USED, &cells);
  ASSERT(cells.size() == 1);
  std::cout << "Surface cell 3 is subsurface face " << f
            << " with internal cell " << cells[0] << std::endl;
#endif

  // Update the preconditioner with darcy and gravity fluxes
  // skip accumulation terms, they're not needed
  mfd_preconditioner_->CreateMFDrhsVectors();

  // Assemble
  mfd_preconditioner_->ApplyBoundaryConditions(bc_markers_, bc_values_);
  mfd_preconditioner_->AssembleGlobalMatrices();

  // derive the consistent faces, involves a solve
  mfd_preconditioner_->UpdateConsistentFaceConstraints(pres_elev.ptr());

  // back out heights from pres_elev
  const Epetra_MultiVector& elevation = *S_next_->GetFieldData("elevation")
      ->ViewComponent("face",false);
  u->ViewComponent("face",false)->Update(1., *pres_elev->ViewComponent("face",false),
          -1., elevation, 0.);

#if DEBUG_FLAG
  *out_ << "  cond = " << (*cond)("face",11) << std::endl;
  *out_ << "  pres_e_f = " << (*pres_elev)("face",11) << std::endl;
  *out_ << "  u_f = " << (*u)("face",11) << std::endl;
#endif
}

} // namespace
} // namespace

