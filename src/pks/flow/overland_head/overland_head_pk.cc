/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -----------------------------------------------------------------------------
This is the overland flow component of ATS.
License: BSD
Author: Ethan Coon (ecoon@lanl.gov)
----------------------------------------------------------------------------- */

#include "EpetraExt_MultiVectorOut.h"
#include "Epetra_MultiVector.h"

#include "flow_bc_factory.hh"
#include "Mesh.hh"
#include "Point.hh"

#include "composite_vector_function.hh"
#include "composite_vector_function_factory.hh"
#include "independent_variable_field_evaluator.hh"

#include "MatrixMFD_TPFA_ScaledConstraint.hh"
#include "upwind_potential_difference.hh"
#include "upwind_total_flux.hh"
#include "upwind_cell_centered.hh"
#include "pres_elev_evaluator.hh"
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
#define DEBUG_RES_FLAG 0

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
    S->AliasMesh("domain", "surface");
  } else {
    standalone_mode_ = false;
  }

  PKPhysicalBDFBase::setup(S);
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

  S->RequireField(key_, name_)->SetMesh(mesh_)
    ->SetComponents(names2, locations2, num_dofs2);

#if DEBUG_RES_FLAG
  // -- residuals of various iterations for debugging
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

  // -- cell volume and evaluator
  S->RequireFieldEvaluator("surface_cell_volume");
  S->RequireGravity();
  S->RequireScalar("atmospheric_pressure");

  // -- coupling to subsurface
  coupled_to_subsurface_via_flux_ =
      plist_.get<bool>("coupled to subsurface via flux", false);
  coupled_to_subsurface_via_head_ =
      plist_.get<bool>("coupled to subsurface via head", false);
  ASSERT(!(coupled_to_subsurface_via_flux_ && coupled_to_subsurface_via_head_));

  if (coupled_to_subsurface_via_head_) {
    // -- source term from subsurface, filled in by evaluator,
    //    which picks the fluxes from "darcy_flux" field.
    S->RequireFieldEvaluator("surface_subsurface_flux");
    S->RequireField("surface_subsurface_flux")
        ->SetMesh(mesh_)->AddComponent("cell", AmanziMesh::CELL, 1);
  }

  // Create the upwinding method.
  S->RequireField("upwind_overland_conductivity", name_)->SetMesh(mesh_)
    ->SetGhosted()->SetComponents(names2, locations2, num_dofs2);
  S->GetField("upwind_overland_conductivity",name_)->set_io_vis(false);

  std::string method_name = plist_.get<string>("upwind conductivity method",
          "upwind by potential difference");
  if (method_name == "cell centered") {
    upwind_method_ = Operators::UPWIND_METHOD_CENTERED;
    upwinding_ = Teuchos::rcp(new Operators::UpwindCellCentered(name_,
            "overland_conductivity", "upwind_overland_conductivity"));
  } else if (method_name == "upwind with total flux") {
    upwind_method_ = Operators::UPWIND_METHOD_TOTAL_FLUX;
    upwinding_ = Teuchos::rcp(new Operators::UpwindTotalFlux(name_,
          "overland_conductivity", "upwind_overland_conductivity",
            "surface_flux_direction"));
  } else if (method_name == "upwind by potential difference") {
    upwind_method_ = Operators::UPWIND_METHOD_POTENTIAL_DIFFERENCE;
    upwinding_ = Teuchos::rcp(new Operators::UpwindPotentialDifference(name_,
            "overland_conductivity", "upwind_overland_conductivity",
            "pres_elev", "ponded_depth"));
  } else {
    std::stringstream messagestream;
    messagestream << "Overland FLow PK has no upwinding method named: " << method_name;
    Errors::Message message(messagestream.str());
    Exceptions::amanzi_throw(message);
  }

  // -- owned secondary variables, no evaluator used
  if (upwind_method_ == Operators::UPWIND_METHOD_TOTAL_FLUX) {
    S->RequireField("surface_flux_direction", name_)->SetMesh(mesh_)->SetGhosted()
        ->SetComponent("face", AmanziMesh::FACE, 1);
  }
  S->RequireField("surface_flux", name_)->SetMesh(S->GetMesh("surface"))
                ->SetGhosted()->SetComponent("face", AmanziMesh::FACE, 1);
  S->RequireField("surface_velocity", name_)->SetMesh(S->GetMesh("surface"))
                ->SetComponent("cell", AmanziMesh::CELL, 3);

  // Create the boundary condition data structures.
  Teuchos::ParameterList bc_plist = plist_.sublist("boundary conditions", true);
  FlowBCFactory bc_factory(mesh_, bc_plist);
  bc_pressure_ = bc_factory.CreatePressure();
  bc_head_ = bc_factory.CreateHead();
  bc_zero_gradient_ = bc_factory.CreateZeroGradient();
  bc_flux_ = bc_factory.CreateMassFlux();


  // operator for the diffusion terms: must use ScaledConstraint version
  Teuchos::ParameterList mfd_plist = plist_.sublist("Diffusion");
  tpfa_ = mfd_plist.get<bool>("TPFA", false);
  if (tpfa_) {
    matrix_ = Teuchos::rcp(new Operators::MatrixMFD_TPFA_ScaledConstraint(mfd_plist, mesh_));
  } else {
    matrix_ = Teuchos::rcp(new Operators::MatrixMFD_ScaledConstraint(mfd_plist, mesh_));
  }
  symmetric_ = false;
  matrix_->set_symmetric(symmetric_);
  matrix_->SymbolicAssembleGlobalMatrices();
  matrix_->CreateMFDmassMatrices(Teuchos::null);
  matrix_->InitPreconditioner();

  Teuchos::ParameterList mfd_pc_plist = plist_.sublist("Diffusion PC");
  Teuchos::RCP<Operators::Matrix> precon;
  if (coupled_to_subsurface_via_head_ || coupled_to_subsurface_via_flux_) {
    mfd_pc_plist.set("TPFA", true);
  }

  if (!tpfa_) tpfa_ = mfd_pc_plist.get<bool>("TPFA", false);
  full_jacobian_ = false;
  if (tpfa_) {
    full_jacobian_ = mfd_pc_plist.get<bool>("TPFA use full Jacobian", false);
  }

  if (tpfa_) {
    tpfa_preconditioner_ =
        Teuchos::rcp(new Operators::MatrixMFD_TPFA_ScaledConstraint(mfd_pc_plist, mesh_));
    precon = tpfa_preconditioner_;
  } else {
    precon = Teuchos::rcp(new Operators::MatrixMFD_ScaledConstraint(mfd_pc_plist, mesh_));
  }
  set_preconditioner(precon);
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

  // upwinding by total flux requires a flux each iteration
  if (upwind_method_ == Operators::UPWIND_METHOD_TOTAL_FLUX)
    update_flux_ = UPDATE_FLUX_ITERATION;
};


// -------------------------------------------------------------
// Create the physical evaluators for water content, water
// retention, rel perm, etc, that are specific to Richards.
// -------------------------------------------------------------
void OverlandHeadFlow::SetupPhysicalEvaluators_(const Teuchos::Ptr<State>& S) {
  std::vector<AmanziMesh::Entity_kind> locations2(2), locations_bf(2);
  std::vector<std::string> names2(2), names_bf(2);
  std::vector<int> num_dofs2(2,1);
  locations2[0] = AmanziMesh::CELL;
  locations2[1] = AmanziMesh::FACE;
  names2[0] = "cell";
  names2[1] = "face";
  locations_bf[0] = AmanziMesh::CELL;
  locations_bf[1] = AmanziMesh::BOUNDARY_FACE;
  names_bf[0] = "cell";
  names_bf[1] = "boundary_face";

  // -- evaluator for surface geometry.
  S->RequireField("elevation")->SetMesh(S->GetMesh("surface"))->SetGhosted()
      ->AddComponents(names2, locations2, num_dofs2);
  S->RequireField("slope_magnitude")->SetMesh(S->GetMesh("surface"))
      ->AddComponents(names_bf, locations_bf, num_dofs2);

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

  // -- evaluator for potential field, h + z
  S->RequireField("pres_elev")->SetMesh(S->GetMesh("surface"))->SetGhosted()
                ->AddComponents(names2, locations2, num_dofs2);
  Teuchos::ParameterList pres_elev_plist = plist_.sublist("potential evaluator");
  Teuchos::RCP<FlowRelations::PresElevEvaluator> pres_elev_eval =
      Teuchos::rcp(new FlowRelations::PresElevEvaluator(pres_elev_plist));
  S->SetFieldEvaluator("pres_elev", pres_elev_eval);

  // -- source term evaluator
  if (plist_.isSublist("source evaluator")) {
    is_source_term_ = true;

    Teuchos::ParameterList source_plist = plist_.sublist("source evaluator");
    source_only_if_unfrozen_ = source_plist.get<bool>("source only if unfrozen",false);
    if (source_only_if_unfrozen_) {
      S->RequireScalar("air_temperature");
    }

    source_plist.set("evaluator name", "overland_source");
    S->RequireField("overland_source")->SetMesh(mesh_)
        ->SetGhosted()->AddComponent("cell", AmanziMesh::CELL, 1);
    Teuchos::RCP<FieldEvaluator> source_evaluator =
        Teuchos::rcp(new IndependentVariableFieldEvaluator(source_plist));
    S->SetFieldEvaluator("overland_source", source_evaluator);
  }

  // -- water content
  S->RequireField("surface_water_content")->SetMesh(mesh_)->SetGhosted()
      ->AddComponent("cell", AmanziMesh::CELL, 1);
  Teuchos::ParameterList wc_plist = plist_.sublist("overland water content evaluator");
  Teuchos::RCP<FlowRelations::OverlandHeadWaterContentEvaluator> wc_evaluator =
      Teuchos::rcp(new FlowRelations::OverlandHeadWaterContentEvaluator(wc_plist));
  S->SetFieldEvaluator("surface_water_content", wc_evaluator);

  // -- ponded depth
  S->RequireField("ponded_depth")->SetMesh(mesh_)->SetGhosted()
                ->AddComponents(names2, locations2, num_dofs2);
  Teuchos::ParameterList height_plist = plist_.sublist("ponded depth evaluator");
  Teuchos::RCP<FlowRelations::HeightEvaluator> height_evaluator =
      Teuchos::rcp(new FlowRelations::HeightEvaluator(height_plist));
  S->SetFieldEvaluator("ponded_depth", height_evaluator);
  height_model_ = height_evaluator->get_Model();

  // -- conductivity evaluator
  S->RequireField("overland_conductivity")->SetMesh(mesh_)->SetGhosted()
      ->AddComponents(names_bf, locations_bf, num_dofs2);
  ASSERT(plist_.isSublist("overland conductivity evaluator"));
  Teuchos::ParameterList cond_plist = plist_.sublist("overland conductivity evaluator");
  Teuchos::RCP<FlowRelations::OverlandConductivityEvaluator> cond_evaluator =
      Teuchos::rcp(new FlowRelations::OverlandConductivityEvaluator(cond_plist));
  S->SetFieldEvaluator("overland_conductivity", cond_evaluator);
}


// -------------------------------------------------------------
// Initialize PK
// -------------------------------------------------------------
void OverlandHeadFlow::initialize(const Teuchos::Ptr<State>& S) {
#if DEBUG_RES_FLAG
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

  // initial condition is tricky
  if (!plist_.isSublist("initial condition")) {
    std::stringstream messagestream;
    messagestream << name_ << " has no initial condition parameter list.";
    Errors::Message message(messagestream.str());
    Exceptions::amanzi_throw(message);
  }
  Teuchos::RCP<CompositeVector> head_cv = S->GetFieldData(key_, name_);
  head_cv->PutScalar(0.);

  // Initialize BDF stuff and physical domain stuff.
  PKPhysicalBDFBase::initialize(S);

  // -- set the cell initial condition if it is taken from the subsurface
  Teuchos::ParameterList ic_plist = plist_.sublist("initial condition");
  if (ic_plist.get<bool>("initialize surface head from subsurface",false)) {
    Epetra_MultiVector& head = *head_cv->ViewComponent("cell",false);
    const Epetra_MultiVector& pres = *S->GetFieldData("pressure")
        ->ViewComponent("face",false);

    unsigned int ncells_surface = mesh_->num_entities(AmanziMesh::CELL,AmanziMesh::OWNED);
    for (unsigned int c=0; c!=ncells_surface; ++c) {
      // -- get the surface cell's equivalent subsurface face and neighboring cell
      AmanziMesh::Entity_ID f =
          mesh_->entity_get_parent(AmanziMesh::CELL, c);
      head[0][c] = pres[0][f];
    }

    // mark as initialized
    if (ic_plist.get<bool>("initialize surface head from subsurface",false)) {
      S->GetField(key_,name_)->set_initialized();
    }
  }

  // Initialize BC data structures
  unsigned int nfaces = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::USED);
  bc_markers_.resize(nfaces, Operators::Matrix::MATRIX_BC_NULL);
  bc_values_.resize(nfaces, 0.0);

  // // Initialize elevation, whose faces need to be updated so that h=0 is a
  // // valid solution.
  // S->GetFieldEvaluator("elevation")->HasFieldChanged(S.ptr(), name_);
  // Teuchos::RCP<CompositeVector> cond =
  //   S->GetFieldData("upwind_overland_conductivity", name_);
  // cond->ViewComponent("cell",true)->PutScalar(1.0);
  // cond->ViewComponent("face",true)->PutScalar(1.0);
  // matrix_->CreateMFDstiffnessMatrices(cond.ptr());
  // matrix_->CreateMFDrhsVectors();
  // matrix_->ApplyBoundaryConditions(bc_markers_, bc_values_);
  // matrix_->AssembleGlobalMatrices();

  // Teuchos::RCP<CompositeVector> elev = S->GetFieldData("elevation","elevation");
  // matrix_->UpdateConsistentFaceConstraints(elev.ptr());

  // Initialize BC values
  bc_pressure_->Compute(S->time());
  bc_head_->Compute(S->time());
  bc_zero_gradient_->Compute(S->time());
  bc_flux_->Compute(S->time());

  // Set extra fields as initialized -- these don't currently have evaluators.
  S->GetFieldData("upwind_overland_conductivity",name_)->PutScalar(1.0);
  S->GetField("upwind_overland_conductivity",name_)->set_initialized();
  S->GetField("surface_flux", name_)->set_initialized();
  if (upwind_method_ == Operators::UPWIND_METHOD_TOTAL_FLUX)
    S->GetField("surface_flux_direction", name_)->set_initialized();
  S->GetField("surface_velocity", name_)->set_initialized();
};


// -----------------------------------------------------------------------------
// Update any secondary (dependent) variables given a solution.
//
//   After a timestep is evaluated (or at ICs), there is no way of knowing if
//   secondary variables have been updated to be consistent with the new
//   solution.
// -----------------------------------------------------------------------------
void OverlandHeadFlow::commit_state(double dt, const Teuchos::RCP<State>& S) {
  niter_ = 0;

  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_EXTREME))
    *vo_->os() << "Commiting state." << std::endl;

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
    matrix_->DeriveFlux(*potential, flux.ptr());
  }

  // As a diagnostic, calculate the mass balance error
#if DEBUG_FLAG
  if (S_next_ != Teuchos::null) {
    const Epetra_MultiVector& wc1 =
        *S_next_->GetFieldData("surface_water_content")
            ->ViewComponent("cell",false);
    const Epetra_MultiVector& wc0 =
        *S_->GetFieldData("surface_water_content")
            ->ViewComponent("cell",false);

    Teuchos::RCP<const CompositeVector> flux =
        S->GetFieldData("surface_flux", name_);
    flux->ScatterMasterToGhosted("face");
    const Epetra_MultiVector& flux_f = *flux->ViewComponent("face",true);

    // create work space
    Teuchos::RCP<CompositeVector> error = Teuchos::rcp(
        new CompositeVector(*S_->GetFieldData("surface_water_content")));
    error->PutScalar(0.);

    // source terms
    AddSourceTerms_(error.ptr());
    error->Scale(dt);

    Epetra_MultiVector& error_c = *error->ViewComponent("cell",false);
    for (unsigned int c=0; c!=error->size("cell"); ++c) {
      // accumulation
      error_c[0][c] += wc1[0][c] - wc0[0][c];

      // flux
      AmanziMesh::Entity_ID_List faces;
      std::vector<int> dirs;
      mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
      for (unsigned int n=0; n!=faces.size(); ++n) {
        error_c[0][c] += flux_f[0][faces[n]] * dirs[n] * dt;
      }
    }

    double einf(0.0);
    error->NormInf(&einf);

    // VerboseObject stuff.
    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "Final Mass Balance Error: " << einf << std::endl;
  }
#endif

};


// -----------------------------------------------------------------------------
// Update diagnostics -- used prior to vis.
// -----------------------------------------------------------------------------
void OverlandHeadFlow::calculate_diagnostics(const Teuchos::RCP<State>& S) {
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_EXTREME))
    *vo_->os() << "Calculating diagnostic variables." << std::endl;

  // update the cell velocities
  if (update_flux_ == UPDATE_FLUX_VIS) {
    Teuchos::RCP<CompositeVector> flux = S->GetFieldData("surface_flux",name_);
    Teuchos::RCP<const CompositeVector> conductivity =
        S->GetFieldData("upwind_overland_conductivity");
    matrix_->CreateMFDstiffnessMatrices(conductivity.ptr());

    // derive the fluxes
    S->GetFieldEvaluator("pres_elev")->HasFieldChanged(S.ptr(), name_);
    Teuchos::RCP<const CompositeVector> potential = S->GetFieldData("pres_elev");
    matrix_->DeriveFlux(*potential, flux.ptr());
  }

  if (update_flux_ != UPDATE_FLUX_NEVER) {
    Teuchos::RCP<const CompositeVector> flux = S->GetFieldData("surface_flux");
    Teuchos::RCP<CompositeVector> velocity = S->GetFieldData("surface_velocity", name_);
    matrix_->DeriveCellVelocity(*flux, velocity.ptr());

    S->GetFieldEvaluator("ponded_depth")->HasFieldChanged(S.ptr(), name_);
    const Epetra_MultiVector& depth_c = *S->GetFieldData("ponded_depth")
        ->ViewComponent("cell",false);

    S->GetFieldEvaluator("surface_molar_density_liquid")->HasFieldChanged(S.ptr(), name_);
    const Epetra_MultiVector& nliq_c = *S->GetFieldData("surface_molar_density_liquid")
        ->ViewComponent("cell",false);

    Epetra_MultiVector& vel_c = *velocity->ViewComponent("cell",false);

    unsigned int ncells = velocity->size("cell");
    for (unsigned int c=0; c!=ncells; ++c) {
      vel_c[0][c] /= (std::max( depth_c[0][c] , 1e-7) * nliq_c[0][c]);
      vel_c[1][c] /= (std::max( depth_c[0][c] , 1e-7) * nliq_c[0][c]);
    }
  }
};


// -----------------------------------------------------------------------------
// Use the physical rel perm (on cells) to update a work vector for rel perm.
//
//   This deals with upwinding, etc.
// -----------------------------------------------------------------------------
bool OverlandHeadFlow::UpdatePermeabilityData_(const Teuchos::Ptr<State>& S) {
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_EXTREME))
    *vo_->os() << "  Updating permeability?";

  bool update_perm = S->GetFieldEvaluator("overland_conductivity")
      ->HasFieldChanged(S, name_);
  update_perm |= S->GetFieldEvaluator("ponded_depth")->HasFieldChanged(S, name_);
  update_perm |= S->GetFieldEvaluator("pres_elev")->HasFieldChanged(S, name_);
  update_perm |= perm_update_required_;

  if (update_perm) {
    // Update the perm only if needed.
    perm_update_required_ = false;

    // get upwind conductivity data
    Teuchos::RCP<CompositeVector> uw_cond =
        S->GetFieldData("upwind_overland_conductivity", name_);

    // update the direction
    if (upwind_method_ == Operators::UPWIND_METHOD_TOTAL_FLUX) {
      // update the direction of the flux -- note this is NOT the flux
      Teuchos::RCP<CompositeVector> flux_dir =
          S->GetFieldData("surface_flux_direction", name_);
      const Epetra_MultiVector& n_liq =
          *S->GetFieldData("surface_molar_density_liquid")
              ->ViewComponent("cell",false);

      // Create the stiffness matrix without a rel perm (just n/mu)
      {
        Epetra_MultiVector& uw_cond_c = *uw_cond->ViewComponent("cell",false);
        int ncells = uw_cond_c.MyLength();
        for (unsigned int c=0; c!=ncells; ++c) {
          uw_cond_c[0][c] = n_liq[0][c];
        }
      }
      uw_cond->ViewComponent("face",false)->PutScalar(1.);
      matrix_->CreateMFDstiffnessMatrices(uw_cond.ptr());

      // Derive the flux
      Teuchos::RCP<const CompositeVector> pres_elev = S->GetFieldData("pres_elev");
      matrix_->DeriveFlux(*pres_elev, flux_dir.ptr());
    }

    // get conductivity data
    Teuchos::RCP<const CompositeVector> cond = S->GetFieldData("overland_conductivity");
    const Epetra_MultiVector& cond_bf = *cond->ViewComponent("boundary_face",false);
    const Epetra_MultiVector& cond_c = *cond->ViewComponent("cell",false);


    { // place boundary_faces on faces
      Epetra_MultiVector& uw_cond_f = *uw_cond->ViewComponent("face",false);
      const Epetra_Import& vandelay = mesh_->exterior_face_importer();
      const Epetra_Map& vandelay_map = mesh_->exterior_face_epetra_map();
      uw_cond_f.Export(cond_bf, vandelay, Insert);

      // Patch up zero-gradient case, which should not upwind.
      AmanziMesh::Entity_ID_List cells;
      for (Functions::BoundaryFunction::Iterator bc=bc_zero_gradient_->begin();
           bc!=bc_zero_gradient_->end(); ++bc) {
        int f = bc->first;
        mesh_->face_get_cells(f, AmanziMesh::USED, &cells);
        ASSERT(cells.size() == 1);
        int c = cells[0];
        uw_cond_f[0][f] = cond_c[0][c];
      }
    }

    // Then upwind.  This overwrites the boundary if upwinding says so.
    upwinding_->Update(S);

    { // Scale cells by n_liq
      const Epetra_MultiVector& n_liq =
          *S->GetFieldData("surface_molar_density_liquid")
              ->ViewComponent("cell",false);

      Epetra_MultiVector& uw_cond_c = *uw_cond->ViewComponent("cell",false);
      int ncells = uw_cond_c.MyLength();
      for (unsigned int c=0; c!=ncells; ++c) {
        ASSERT(uw_cond_c[0][c] == 1.); // REMOVE ME!
        uw_cond_c[0][c] *= n_liq[0][c];
      }
    }
  }

  if (update_perm && out_.get() &&
      includesVerbLevel(verbosity_, Teuchos::VERB_EXTREME, true)) {
    *vo_->os() << " TRUE." << std::endl;
  }
  return update_perm;
}


// -----------------------------------------------------------------------------
// Evaluate boundary conditions at the current time.
// -----------------------------------------------------------------------------
void OverlandHeadFlow::UpdateBoundaryConditions_(const Teuchos::Ptr<State>& S) {
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_EXTREME))
    *vo_->os() << "  Updating BCs." << std::endl;

  AmanziMesh::Entity_ID_List cells;
  const Epetra_MultiVector& elevation = *S->GetFieldData("elevation")
      ->ViewComponent("face",false);

  // initialize all as null
  for (unsigned int n=0; n!=bc_markers_.size(); ++n) {
    bc_markers_[n] = Operators::Matrix::MATRIX_BC_NULL;
    bc_values_[n] = 0.0;
  }

  // Head BCs are standard Dirichlet, plus elevation
  for (Functions::BoundaryFunction::Iterator bc=bc_head_->begin();
       bc!=bc_head_->end(); ++bc) {
    int f = bc->first;
    bc_markers_[f] = Operators::Matrix::MATRIX_BC_DIRICHLET;
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

    bc_markers_[f] = Operators::Matrix::MATRIX_BC_DIRICHLET;
    bc_values_[f] = height + elevation[0][f];
  }

  // Standard Neumann data for flux
  for (Functions::BoundaryFunction::Iterator bc=bc_flux_->begin();
       bc!=bc_flux_->end(); ++bc) {
    int f = bc->first;
    bc_markers_[f] = Operators::Matrix::MATRIX_BC_FLUX;
    bc_values_[f] = bc->second;
  }

  // zero gradient: grad h = 0 implies that q = -k grad z
  for (Functions::BoundaryFunction::Iterator bc=bc_zero_gradient_->begin();
       bc!=bc_zero_gradient_->end(); ++bc) {
    int f = bc->first;
    //    bc_markers_[f] = Operators::Matrix::MATRIX_BC_FLUX;
    //    bc_values_[f] = 0.;
    //    bc_markers_[f] = Operators::Matrix::MATRIX_BC_DIRICHLET;

    // mesh_->face_get_cells(f, AmanziMesh::USED, &cells);
    // ASSERT( cells.size()==1 ) ;
    // bc_values_[f] = height[0][cells[0]] + elevation[0][f];
  }
}


// -----------------------------------------------------------------------------
// Evaluate boundary conditions at the current time without elevation.
// -----------------------------------------------------------------------------
void OverlandHeadFlow::UpdateBoundaryConditionsMarkers_(const Teuchos::Ptr<State>& S) {
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_EXTREME))
    *vo_->os() << "  Updating BCs Markers only." << std::endl;

  AmanziMesh::Entity_ID_List cells;
  // initialize all as null
  for (unsigned int n=0; n!=bc_markers_.size(); ++n) {
    bc_markers_[n] = Operators::Matrix::MATRIX_BC_NULL;
    bc_values_[n] = 0.;
  }

  // Head BCs are standard Dirichlet, no elevation
  for (Functions::BoundaryFunction::Iterator bc=bc_head_->begin();
       bc!=bc_head_->end(); ++bc) {
    int f = bc->first;
    bc_markers_[f] = Operators::Matrix::MATRIX_BC_DIRICHLET;
  }

  // pressure BCs require calculating head(pressure)
  for (Functions::BoundaryFunction::Iterator bc=bc_pressure_->begin();
       bc!=bc_pressure_->end(); ++bc) {
    int f = bc->first;
    bc_markers_[f] = Operators::Matrix::MATRIX_BC_DIRICHLET;
  }

  // Standard Neumann data for flux
  for (Functions::BoundaryFunction::Iterator bc=bc_flux_->begin();
       bc!=bc_flux_->end(); ++bc) {
    int f = bc->first;
    bc_markers_[f] = Operators::Matrix::MATRIX_BC_FLUX;
  }

  // zero gradient: grad h = 0 implies that q = -k grad z
  for (Functions::BoundaryFunction::Iterator bc=bc_zero_gradient_->begin();
       bc!=bc_zero_gradient_->end(); ++bc) {
    //    int f = bc->first;
    //    bc_markers_[f] = Operators::Matrix::MATRIX_BC_FLUX;
    //    bc_markers_[f] = Operators::Matrix::MATRIX_BC_DIRICHLET;
  }
}


void OverlandHeadFlow::FixBCsForOperator_(const Teuchos::Ptr<State>& S) {
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_EXTREME))
    *vo_->os() << "    Tweaking BCs for the Operator." << std::endl;

  // Now we can safely calculate q = -k grad z for zero-gradient problems
  Teuchos::RCP<const CompositeVector> elev = S->GetFieldData("elevation");
  elev->ScatterMasterToGhosted();
  const Epetra_MultiVector& elevation_f = *elev->ViewComponent("face",false);
  const Epetra_MultiVector& elevation_c = *elev->ViewComponent("cell",false);

  std::vector<Teuchos::SerialDenseMatrix<int, double> >& Aff_cells =
      matrix_->Aff_cells();
  std::vector<Epetra_SerialDenseVector>& Ff_cells =
      matrix_->Ff_cells();

  int ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  for (Functions::BoundaryFunction::Iterator bc=bc_zero_gradient_->begin();
       bc!=bc_zero_gradient_->end(); ++bc) {

    int f = bc->first;

    AmanziMesh::Entity_ID_List cells;
    mesh_->face_get_cells(f, AmanziMesh::USED, &cells);
    ASSERT(cells.size() == 1);
    AmanziMesh::Entity_ID c = cells[0];

    if (c < ncells_owned) {
      AmanziMesh::Entity_ID_List faces;
      std::vector<int> dirs;
      mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);

      std::vector<double> dp(faces.size());
      for (unsigned int n=0; n!=faces.size(); ++n) {
        dp[n] = elevation_f[0][faces[n]] - elevation_c[0][c];
      }
      unsigned int my_n = std::find(faces.begin(), faces.end(), f) - faces.begin();
      ASSERT(my_n !=faces.size());

      double bc_val = 0.;
      for (unsigned int m=0; m!=faces.size(); ++m) {
        bc_val -= Aff_cells[c](my_n,m) * dp[m];
      }

      // Apply the BC to the matrix
      Ff_cells[c][my_n] -= bc_val;
    }
  }
};


void OverlandHeadFlow::FixBCsForPrecon_(const Teuchos::Ptr<State>& S) {
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_EXTREME))
    *vo_->os() << "    Tweaking BCs for the PC." << std::endl;

  // // Attempt of a hack to deal with zero rel perm
  // double eps = 1.e-30;
  // Teuchos::RCP<CompositeVector> relperm =
  //     S->GetFieldData("upwind_overland_conductivity", name_);
  // for (unsigned int f=0; f!=relperm->size("face"); ++f) {
  //   if ((*relperm)("face",f) < eps) {
  //     if (bc_markers_[f] == Operators::Matrix::MATRIX_BC_FLUX) {
  //       bc_markers_[f] = Operators::Matrix::MATRIX_BC_DIRICHLET;
  //     } else if (bc_markers_[f] == Operators::Matrix::MATRIX_BC_NULL) {
  //       bc_markers_[f] = Operators::Matrix::MATRIX_BC_DIRICHLET;
  //     }
  //   }
  // }
};

void OverlandHeadFlow::FixBCsForConsistentFaces_(const Teuchos::Ptr<State>& S) {
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_EXTREME))
    *vo_->os() << "    Tweaking BCs for calculation of consistent faces." << std::endl;

  // // If the rel perm is 0, the face value drops out and is unconstrained.
  // // Therefore we set it to Dirichlet to eliminate it from the system.
  // double eps = 1.e-30;
  // const Epetra_MultiVector& elevation = *S->GetFieldData("elevation")
  //     ->ViewComponent("face",false);
  // Teuchos::RCP<CompositeVector> relperm =
  //     S->GetFieldData("upwind_overland_conductivity", name_);

  // for (unsigned int f=0; f!=relperm->size("face"); ++f) {
  //   if ((*relperm)("face",f) < eps) {
  //     if (bc_markers_[f] == Operators::Matrix::MATRIX_BC_FLUX) {
  //       bc_markers_[f] = Operators::Matrix::MATRIX_BC_DIRICHLET;
  //       bc_values_[f] =  elevation[0][f];
  //     } else if (bc_markers_[f] == Operators::Matrix::MATRIX_BC_NULL) {
  //       bc_markers_[f] = Operators::Matrix::MATRIX_BC_DIRICHLET;
  //       bc_values_[f] =  elevation[0][f];
  //     }
  //   }
  // }

  // Now we can safely calculate q = -k grad z for zero-gradient problems
  Teuchos::RCP<const CompositeVector> elev = S->GetFieldData("elevation");
  elev->ScatterMasterToGhosted();
  const Epetra_MultiVector& elevation_f = *elev->ViewComponent("face",false);
  const Epetra_MultiVector& elevation_c = *elev->ViewComponent("cell",false);

  std::vector<Teuchos::SerialDenseMatrix<int, double> >& Aff_cells =
      matrix_->Aff_cells();
  std::vector<Epetra_SerialDenseVector>& Ff_cells =
      matrix_->Ff_cells();

  int ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  for (Functions::BoundaryFunction::Iterator bc=bc_zero_gradient_->begin();
       bc!=bc_zero_gradient_->end(); ++bc) {
    int f = bc->first;

    AmanziMesh::Entity_ID_List cells;
    mesh_->face_get_cells(f, AmanziMesh::USED, &cells);
    ASSERT(cells.size() == 1);
    AmanziMesh::Entity_ID c = cells[0];

    if (c < ncells_owned) {
      AmanziMesh::Entity_ID_List faces;
      std::vector<int> dirs;
      mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);

      std::vector<double> dp(faces.size());
      for (unsigned int n=0; n!=faces.size(); ++n) {
        dp[n] = elevation_f[0][faces[n]] - elevation_c[0][c];
      }
      unsigned int my_n = std::find(faces.begin(), faces.end(), f) - faces.begin();
      ASSERT(my_n !=faces.size());

      double bc_val = 0.;
      for (unsigned int m=0; m!=faces.size(); ++m) {
        bc_val -= Aff_cells[c](my_n,m) * dp[m];
      }

      // Apply the BC to the matrix
      Ff_cells[c][my_n] -= bc_val;
    }
  }
};


/* ******************************************************************
 * Add a boundary marker to owned faces.
 ****************************************************************** */
void OverlandHeadFlow::ApplyBoundaryConditions_(const Teuchos::RCP<State>& S,
        const Teuchos::RCP<CompositeVector>& pres) {
  Epetra_MultiVector& pres_f = *pres->ViewComponent("face",false);
  unsigned int nfaces = pres_f.MyLength();
  for (unsigned int f=0; f!=nfaces; ++f) {
    if (bc_markers_[f] == Operators::Matrix::MATRIX_BC_DIRICHLET) {
      pres_f[0][f] = bc_values_[f];
    }
  }
};


bool OverlandHeadFlow::modify_predictor(double h, Teuchos::RCP<TreeVector> u) {
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_EXTREME))
    *vo_->os() << "Modifying predictor:" << std::endl;

  if (modify_predictor_with_consistent_faces_) {
    if (vo_->os_OK(Teuchos::VERB_EXTREME))
      *vo_->os() << "  modifications for consistent face pressures." << std::endl;
    CalculateConsistentFaces(u->Data().ptr());
    return true;
  }

  return false;
};


void OverlandHeadFlow::CalculateConsistentFaces(const Teuchos::Ptr<CompositeVector>& u) {
  // VerboseObject stuff.
  Teuchos::OSTab tab = vo_->getOSTab();

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
  matrix_->CreateMFDstiffnessMatrices(cond.ptr());
  matrix_->CreateMFDrhsVectors();

  // Patch up BCs in the case of zero conductivity
  FixBCsForConsistentFaces_(S_next_.ptr());

  // Grab needed data.
  S_next_->GetFieldEvaluator("pres_elev")->HasFieldChanged(S_next_.ptr(), name_);
  Teuchos::RCP<CompositeVector> pres_elev = S_next_->GetFieldData("pres_elev","pres_elev");

  // Update the preconditioner with darcy and gravity fluxes
  // skip accumulation terms, they're not needed
  // Assemble
  matrix_->ApplyBoundaryConditions(bc_markers_, bc_values_);
  matrix_->AssembleGlobalMatrices();

  // derive the consistent faces, involves a solve
  matrix_->UpdateConsistentFaceConstraints(pres_elev.ptr());

  // back out heights from pres_elev
  const Epetra_MultiVector& elevation = *S_next_->GetFieldData("elevation")
      ->ViewComponent("face",false);
  u->ViewComponent("face",false)->Update(1., *pres_elev->ViewComponent("face",false),
          -1., elevation, 0.);
}

} // namespace
} // namespace

