/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -----------------------------------------------------------------------------
This is the overland flow component of ATS.
License: BSD
Authors: Gianmarco Manzini
         Ethan Coon (ecoon@lanl.gov)
----------------------------------------------------------------------------- */

#include "EpetraExt_MultiVectorOut.h"
#include "Epetra_MultiVector.h"

#include "bdf1_time_integrator.hh"
#include "flow_bc_factory.hh"
#include "Mesh.hh"
#include "Mesh_MSTK.hh"
#include "Point.hh"

#include "composite_vector_function.hh"
#include "composite_vector_function_factory.hh"
#include "independent_variable_field_evaluator.hh"

#include "upwinding.hh"
#include "upwind_potential_difference.hh"
#include "elevation_evaluator.hh"
#include "meshed_elevation_evaluator.hh"
#include "standalone_elevation_evaluator.hh"
#include "overland_conductivity_evaluator.hh"
#include "overland_conductivity_model.hh"
#include "overland_head_water_content_evaluator.hh"
#include "height_evaluator.hh"

#include "overland_head.hh"

namespace Amanzi {
namespace Flow {

RegisteredPKFactory<OverlandHeadFlow> OverlandHeadFlow::reg_("overland flow, head basis");

// -------------------------------------------------------------
// Constructor
// -------------------------------------------------------------
void OverlandHeadFlow::setup(const Teuchos::Ptr<State>& S) {
  PKPhysicalBDFBase::setup(S);
  CreateMesh_(S);
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
  S->RequireField("overland_flux", name_)->SetMesh(S->GetMesh("surface"))
                ->SetGhosted()->SetComponent("face", AmanziMesh::FACE, 1);
  S->RequireField("overland_velocity", name_)->SetMesh(S->GetMesh("surface"))
                ->SetGhosted()->SetComponent("cell", AmanziMesh::CELL, 3);

  // -- cell volume and evaluator
  S->RequireFieldEvaluator("surface_cell_volume");
  S->RequireGravity();
  S->RequireScalar("atmospheric_pressure");

  // Create the boundary condition data structures.
  Teuchos::ParameterList bc_plist = plist_.sublist("boundary conditions", true);
  FlowBCFactory bc_factory(mesh_, bc_plist);
  bc_pressure_ = bc_factory.CreatePressure();
  bc_zero_gradient_ = bc_factory.CreateZeroGradient();
  bc_flux_ = bc_factory.CreateMassFlux();

  // coupling to subsurface
  is_coupling_term_ = plist_.get<bool>("coupled to subsurface", false);

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
  bool symmetric = false;
  matrix_->SetSymmetryProperty(symmetric);
  matrix_->SymbolicAssembleGlobalMatrices();

  // preconditioner for the NKA system
  Teuchos::ParameterList mfd_pc_plist = plist_.sublist("Diffusion PC");
  preconditioner_ = Teuchos::rcp(new Operators::MatrixMFD(mfd_pc_plist, mesh_));
  preconditioner_->SetSymmetryProperty(symmetric);
  preconditioner_->SymbolicAssembleGlobalMatrices();
  preconditioner_->InitPreconditioner(mfd_pc_plist);
  assemble_preconditioner_ = plist_.get<bool>("assemble preconditioner", true);

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
                ->SetComponents(names2, locations2, num_dofs2);
  S->RequireField("slope_magnitude")->SetMesh(S->GetMesh("surface"))
                ->SetGhosted()->SetComponent("cell", AmanziMesh::CELL, 1);
  S->RequireField("pres_elev")->SetMesh(S->GetMesh("surface"))->SetGhosted()
                ->SetComponents(names2, locations2, num_dofs2);

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
    ->SetComponent("cell", AmanziMesh::CELL, 1);
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
        ->SetGhosted()->SetComponent("cell", AmanziMesh::CELL, 1);
    Teuchos::RCP<FieldEvaluator> source_evaluator =
        Teuchos::rcp(new IndependentVariableFieldEvaluator(source_plist));
    S->SetFieldEvaluator("overland_source", source_evaluator);
  }

  // -- coupling term, filled in by subsurface PK
  if (is_coupling_term_) {
    S->RequireField("overland_source_from_subsurface")
        ->SetMesh(mesh_)->SetComponent("cell", AmanziMesh::CELL, 1);
    S->RequireField("doverland_source_from_subsurface_dsurface_pressure")
        ->SetMesh(mesh_)->SetComponent("cell", AmanziMesh::CELL, 1);
  }

  // -- water content
  S->RequireField("surface_water_content")->SetMesh(mesh_)->SetGhosted()
      ->SetComponent("cell", AmanziMesh::CELL, 1);
  Teuchos::ParameterList wc_plist = plist_.sublist("overland water content");
  Teuchos::RCP<FlowRelations::OverlandHeadWaterContentEvaluator> wc_evaluator =
      Teuchos::rcp(new FlowRelations::OverlandHeadWaterContentEvaluator(wc_plist));
  S->SetFieldEvaluator("surface_water_content", wc_evaluator);

  // -- ponded depth
  S->RequireField("ponded_depth")->SetMesh(mesh_)->SetGhosted()
                ->SetComponents(names2, locations2, num_dofs2);
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

    Teuchos::RCP<const AmanziMesh::Mesh_MSTK> domain_mesh =
        Teuchos::rcp_dynamic_cast<const AmanziMesh::Mesh_MSTK>(S->GetMesh());
    Teuchos::RCP<const AmanziMesh::Mesh_MSTK> surface_mesh =
        Teuchos::rcp_dynamic_cast<const AmanziMesh::Mesh_MSTK>(mesh_);

    int ncells_surface = surface_mesh->num_entities(AmanziMesh::CELL,AmanziMesh::OWNED);
    for (int c=0; c!=ncells_surface; ++c) {
      // -- get the surface cell's equivalent subsurface face and neighboring cell
      AmanziMesh::Entity_ID f =
          surface_mesh->entity_get_parent(AmanziMesh::CELL, c);
      head[0][c] = pres[0][f];
    }
    S->GetField(key_,name_)->set_initialized();
  }
  S->GetFieldData(key_,name_)->ViewComponent("face",false)->PutScalar(0.);

  // Initialize BDF stuff and physical domain stuff.
  PKPhysicalBDFBase::initialize(S);

  // Initialize boundary conditions.
  int nfaces = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::USED);
  bc_markers_.resize(nfaces, Operators::MFD_BC_NULL);
  bc_values_.resize(nfaces, 0.0);

  bc_pressure_->Compute(S->time());
  bc_zero_gradient_->Compute(S->time());
  bc_flux_->Compute(S->time());
  //  UpdateBoundaryConditions_(S);

  // Set extra fields as initialized -- these don't currently have evaluators.
  S->GetFieldData("upwind_overland_conductivity",name_)->PutScalar(1.0);
  S->GetField("upwind_overland_conductivity",name_)->set_initialized();
  S->GetField("overland_flux", name_)->set_initialized();
  S->GetField("overland_velocity", name_)->set_initialized();

  // Initialize operators.
  matrix_->CreateMFDmassMatrices(Teuchos::null);
  preconditioner_->CreateMFDmassMatrices(Teuchos::null);
};


void OverlandHeadFlow::CreateMesh_(const Teuchos::Ptr<State>& S) {
  bool running_timers(false);

  // Create mesh
  if (S->GetMesh()->space_dimension() == 3) {
    // The domain mesh is a 3D volume mesh or a 3D surface mesh -- construct
    // the surface mesh.
    // -- Ensure the domain mesh is MSTK.
    Teuchos::RCP<const AmanziMesh::Mesh_MSTK> mesh =
      Teuchos::rcp_dynamic_cast<const AmanziMesh::Mesh_MSTK>(S->GetMesh());

    if (mesh == Teuchos::null) {
      Errors::Message message("Overland Flow PK requires surface mesh, which is currently only supported by MSTK.  Make the domain mesh an MSTK mesh.");
      Exceptions::amanzi_throw(message);
    }

    // -- Start the clock
    running_timers = true;
    //    Teuchos::RCP<Teuchos::Time> meshtime = Teuchos::TimeMonitor::getNewCounter("surface mesh creation");
    //    Teuchos::TimeMonitor timer(*meshtime);

    // -- Check that the surface mesh has a subset
    std::vector<std::string> setnames;
    if (plist_.isParameter("surface sideset name")) {
      setnames.push_back(plist_.get<std::string>("surface sideset name"));
    } else {
      setnames = plist_.get<Teuchos::Array<std::string> >("surface sideset names").toVector();
    }

    // -- Call the MSTK constructor to rip off the surface of the MSTK domain
    // -- mesh.
    Teuchos::RCP<AmanziMesh::Mesh> surface_mesh;

    if (mesh->cell_dimension() == 3) {
      Teuchos::RCP<AmanziMesh::Mesh> surface_mesh_3d = Teuchos::rcp(
          new AmanziMesh::Mesh_MSTK(*mesh,setnames,AmanziMesh::FACE,false,false));
      S->RegisterMesh("surface_3d", surface_mesh_3d);

      surface_mesh = Teuchos::rcp(
          new AmanziMesh::Mesh_MSTK(*mesh,setnames,AmanziMesh::FACE,true,false));
    } else {
      S->RegisterMesh("surface_3d", mesh);

      surface_mesh = Teuchos::rcp(
          new AmanziMesh::Mesh_MSTK(*mesh,setnames,AmanziMesh::CELL,true,false));
    }

    // -- push the mesh into state
    S->RegisterMesh("surface", surface_mesh);
    mesh_ = surface_mesh;
    standalone_mode_ = false;

  } else if (S->GetMesh()->space_dimension() == 2) {
    // The domain mesh is already a 2D mesh, so simply register it as the surface
    // mesh as well.
    S->RegisterMesh("surface", S->GetMesh());
    mesh_ = S->GetMesh();
    standalone_mode_ = true;
  } else {
    Errors::Message message("Invalid mesh dimension for overland flow.");
    Exceptions::amanzi_throw(message);
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
    Teuchos::RCP<CompositeVector> flux = S->GetFieldData("overland_flux", name_);
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
    Teuchos::RCP<CompositeVector> flux = S->GetFieldData("overland_flux",name_);
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
    Teuchos::RCP<const CompositeVector> flux = S->GetFieldData("overland_flux");
    Teuchos::RCP<CompositeVector> velocity = S->GetFieldData("overland_velocity", name_);
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

  if (update_perm) {
    // Update the perm only if needed.

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

    // First do the boundary
    AmanziMesh::Entity_ID_List cells;
    int nfaces = upwind_conductivity_f.MyLength();
    for (int f=0; f!=nfaces; ++f) {
      if (bc_markers_[f] != Operators::MFD_BC_NULL) {
        upwind_conductivity->mesh()->face_get_cells(f, AmanziMesh::USED, &cells);
        ASSERT(cells.size() == 1);
        int c = cells[0];
        upwind_conductivity_f[0][f] =
          cond_model_->Conductivity(depth[0][f], slope[0][c], coef[0][c]);
      }
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

  const Epetra_MultiVector& elevation = *S->GetFieldData("elevation")
      ->ViewComponent("face",false);
  S->GetFieldEvaluator("ponded_depth")->HasFieldChanged(S.ptr(), name_);
  const Epetra_MultiVector& height = *S->GetFieldData("ponded_depth")
      ->ViewComponent("cell",false);

  for (int n=0; n!=bc_markers_.size(); ++n) {
    bc_markers_[n] = Operators::MFD_BC_NULL;
    bc_values_[n] = 0.0;
  }

  Functions::BoundaryFunction::Iterator bc;
  int count = 0;
  for (bc=bc_pressure_->begin(); bc!=bc_pressure_->end(); ++bc) {
    int f = bc->first;
    bc_markers_[f] = Operators::MFD_BC_DIRICHLET;
    bc_values_[f] = bc->second + elevation[0][f];
    count ++;
  }

  AmanziMesh::Entity_ID_List cells;
  for (bc=bc_zero_gradient_->begin(); bc!=bc_zero_gradient_->end(); ++bc) {
    int f = bc->first;

    cells.clear();
    S->GetMesh("surface")->face_get_cells(f, AmanziMesh::USED, &cells);
    int ncells = cells.size();
    ASSERT( ncells==1 ) ;

    bc_markers_[f] = Operators::MFD_BC_DIRICHLET;
    bc_values_[f] = height[0][cells[0]] + elevation[0][f];
  }

  for (bc=bc_flux_->begin(); bc!=bc_flux_->end(); ++bc) {
    int f = bc->first;
    bc_markers_[f] = Operators::MFD_BC_FLUX;
    bc_values_[f] = bc->second;
  }
};

void OverlandHeadFlow::UpdateBoundaryConditionsNoElev_(const Teuchos::Ptr<State>& S) {
  S->GetFieldEvaluator("ponded_depth")->HasFieldChanged(S.ptr(), name_);
  const Epetra_MultiVector& height = *S->GetFieldData("ponded_depth")
      ->ViewComponent("cell",false);

  for (int n=0; n!=bc_markers_.size(); ++n) {
    bc_markers_[n] = Operators::MFD_BC_NULL;
    bc_values_[n] = 0.0;
  }

  Functions::BoundaryFunction::Iterator bc;
  int count = 0;
  for (bc=bc_pressure_->begin(); bc!=bc_pressure_->end(); ++bc) {
    int f = bc->first;
    bc_markers_[f] = Operators::MFD_BC_DIRICHLET;
    bc_values_[f] = bc->second;
    count ++;
  }

  AmanziMesh::Entity_ID_List cells;
  for (bc=bc_zero_gradient_->begin(); bc!=bc_zero_gradient_->end(); ++bc) {
    int f = bc->first;

    cells.clear();
    S->GetMesh("surface")->face_get_cells(f, AmanziMesh::USED, &cells);
    int ncells = cells.size();
    ASSERT( ncells==1 ) ;

    bc_markers_[f] = Operators::MFD_BC_DIRICHLET;
    bc_values_[f] = height[0][cells[0]];
  }

  for (bc=bc_flux_->begin(); bc!=bc_flux_->end(); ++bc) {
    int f = bc->first;
    bc_markers_[f] = Operators::MFD_BC_FLUX;
    bc_values_[f] = bc->second;
  }

  // Attempt of a hack to deal with zero rel perm
  double eps = 1.e-12;
  UpdatePermeabilityData_(S);
  Teuchos::RCP<CompositeVector> relperm =
      S->GetFieldData("upwind_overland_conductivity", name_);
  for (int f=0; f!=relperm->size("face"); ++f) {
    if ((*relperm)("face",f) < eps && bc_markers_[f] == Operators::MFD_BC_NULL) {
      bc_markers_[f] = Operators::MFD_BC_DIRICHLET;
      bc_values_[f] = 0.0;
      //      (*relperm)("face",f) = 1.0;
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
void OverlandHeadFlow::changed_solution() {
  solution_evaluator_->SetFieldAsChanged();
  // communicate both faces and cells -- faces for flux and cells for rel perm
  // upwinding overlap
  S_next_->GetFieldData(key_)->ScatterMasterToGhosted();
};

} // namespace
} // namespace

