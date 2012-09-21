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

#include "upwind_potential_difference.hh"
#include "elevation_evaluator.hh"
#include "meshed_elevation_evaluator.hh"
#include "standalone_elevation_evaluator.hh"
#include "manning_conductivity_evaluator.hh"

#include "overland.hh"

namespace Amanzi {
namespace Flow {

RegisteredPKFactory<OverlandFlow> OverlandFlow::reg_("overland flow");

// -------------------------------------------------------------
// Constructor
// -------------------------------------------------------------
void OverlandFlow::setup(const Teuchos::Ptr<State>& S) {
  PKPhysicalBDFBase::setup(S);
  CreateMesh_(S);
  SetupOverlandFlow_(S);
  SetupPhysicalEvaluators_(S);
}


void OverlandFlow::SetupOverlandFlow_(const Teuchos::Ptr<State>& S) {
  // Require fields and evaluators for those fields.
  // -- primary variable: pressure on both cells and faces, ghosted, with 1 dof
  std::vector<AmanziMesh::Entity_kind> locations2(2);
  std::vector<std::string> names2(2);
  std::vector<int> num_dofs2(2,1);
  locations2[0] = AmanziMesh::CELL;
  locations2[1] = AmanziMesh::FACE;
  names2[0] = "cell";
  names2[1] = "face";

  S->RequireField(key_, name_)->SetMesh(S->GetMesh("surface"))
                ->SetGhosted()->SetComponents(names2, locations2, num_dofs2);

  // -- owned secondary variables, no evaluator used
  S->RequireField("overland_flux", name_)->SetMesh(S->GetMesh("surface"))
                ->SetGhosted()->SetComponent("face", AmanziMesh::FACE, 1);
  S->RequireField("overland_velocity", name_)->SetMesh(S->GetMesh("surface"))
                ->SetGhosted()->SetComponent("cell", AmanziMesh::CELL, 3);

  // boundary conditions
  Teuchos::ParameterList bc_plist = plist_.sublist("boundary conditions", true);
  FlowBCFactory bc_factory(S->GetMesh("surface"), bc_plist);
  bc_pressure_ = bc_factory.CreatePressure();
  bc_zero_gradient_ = bc_factory.CreateZeroGradient();
  bc_flux_ = bc_factory.CreateMassFlux();

  // rel perm upwinding
  S->RequireField("upwind_overland_conductivity", name_)
                ->SetMesh(S->GetMesh("surface"))->SetGhosted()
                ->SetComponents(names2, locations2, num_dofs2);
  S->GetField("upwind_overland_conductivity",name_)->set_io_vis(false);

  upwinding_ = Teuchos::rcp(new Operators::UpwindPotentialDifference(name_,
          "overland_conductivity", "upwind_overland_conductivity",
          "pres_elev", key_));

  // operator for the diffusion terms
  Teuchos::ParameterList mfd_plist = plist_.sublist("Diffusion");
  matrix_ = Teuchos::rcp(new Operators::MatrixMFD(mfd_plist, S->GetMesh("surface")));
  bool symmetric = false;
  matrix_->SetSymmetryProperty(symmetric);
  matrix_->SymbolicAssembleGlobalMatrices();

  // preconditioner for the NKA system
  Teuchos::ParameterList mfd_pc_plist = plist_.sublist("Diffusion PC");
  preconditioner_ = Teuchos::rcp(new Operators::MatrixMFD(mfd_pc_plist, S->GetMesh("surface")));
  preconditioner_->SetSymmetryProperty(symmetric);
  preconditioner_->SymbolicAssembleGlobalMatrices();
  preconditioner_->InitPreconditioner(mfd_pc_plist);

  steptime_ = Teuchos::TimeMonitor::getNewCounter("overland flow advance time");

};


void OverlandFlow::SetupPhysicalEvaluators_(const Teuchos::Ptr<State>& S) {
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

  // -- "rel perm" evaluator
  S->RequireField("overland_conductivity")->SetMesh(S->GetMesh("surface"))
                ->SetGhosted()->SetComponent("cell", AmanziMesh::CELL, 1);
  ASSERT(plist_.isSublist("overland conductivity evaluator"));
  Teuchos::ParameterList cond_plist = plist_.sublist("overland conductivity evaluator");
  Teuchos::RCP<FieldEvaluator> cond_evaluator =
      Teuchos::rcp(new FlowRelations::ManningConductivityEvaluator(cond_plist));
  S->SetFieldEvaluator("overland_conductivity", cond_evaluator);
  // -- -- hack to deal with BCs of the rel perm upwinding... this will need
  // -- -- some future thought --etc
  manning_exp_ = cond_plist.get<double>("Manning exponent", 0.6666666666666667);
  slope_regularization_ = cond_plist.get<double>("slope regularization epsilon", 1.e-8);

  // -- source term evaluator
  if (plist_.isSublist("source evaluator")) {
    Teuchos::ParameterList source_plist = plist_.sublist("source evaluator");
    source_plist.set("evaluator name", "overland_source");
    is_source_term_ = true;
    S->RequireField("overland_source")->SetMesh(S->GetMesh("surface"))
        ->SetGhosted()->SetComponent("cell", AmanziMesh::CELL, 1);
    Teuchos::RCP<FieldEvaluator> source_evaluator =
        Teuchos::rcp(new IndependentVariableFieldEvaluator(source_plist));
    S->SetFieldEvaluator("overland_source", source_evaluator);
  }

  // -- coupling term evaluator
  if (plist_.get<bool>("coupled to subsurface", false)) {
    is_coupling_term_ = true;
    S->RequireField("overland_source_from_subsurface", name_)
        ->SetMesh(S->GetMesh("surface"))->SetGhosted()
        ->SetComponent("cell", AmanziMesh::CELL, 1);
    S->RequireFieldEvaluator("overland_coupling");
  }


  // -- cell volume and evaluator
  S->RequireField("surface_cell_volume")->SetMesh(S->GetMesh("surface"))->SetGhosted()
                                ->AddComponent("cell", AmanziMesh::CELL, 1);
  S->RequireFieldEvaluator("surface_cell_volume");
}

// -------------------------------------------------------------
// Initialize PK
// -------------------------------------------------------------
void OverlandFlow::initialize(const Teuchos::Ptr<State>& S) {
  // initialize BDF stuff and physical domain stuff
  PKPhysicalBDFBase::initialize(S);

  // initialize boundary conditions
  int nfaces = S->GetMesh("surface")->num_entities(AmanziMesh::FACE, AmanziMesh::USED);
  bc_markers_.resize(nfaces, Operators::MFD_BC_NULL);
  bc_values_.resize(nfaces, 0.0);

  bc_pressure_->Compute(S->time());
  bc_zero_gradient_->Compute(S->time());
  bc_flux_->Compute(S->time());
  UpdateBoundaryConditions_(S);

  // Set extra fields as initialized -- these don't currently have evaluators.
  S->GetFieldData("upwind_overland_conductivity",name_)->PutScalar(1.0);
  S->GetField("upwind_overland_conductivity",name_)->set_initialized();
  S->GetField("overland_flux", name_)->set_initialized();
  S->GetField("overland_velocity", name_)->set_initialized();

  // Initialize operators.
  matrix_->CreateMFDmassMatrices(Teuchos::null);
  preconditioner_->CreateMFDmassMatrices(Teuchos::null);
};


void OverlandFlow::CreateMesh_(const Teuchos::Ptr<State>& S) {
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
    Teuchos::RCP<Teuchos::Time> meshtime = Teuchos::TimeMonitor::getNewCounter("surface mesh creation");
    Teuchos::TimeMonitor timer(*meshtime);

    // Check that the surface mesh has a subset
    std::string surface_sideset_name =
      plist_.get<std::string>("Surface sideset name");

    // -- Call the MSTK constructor to rip off the surface of the MSTK domain
    // -- mesh.
    std::vector<std::string> setnames(1,surface_sideset_name);
    Teuchos::RCP<AmanziMesh::Mesh> surface_mesh_3d =
      Teuchos::rcp(new AmanziMesh::Mesh_MSTK(*mesh,setnames,AmanziMesh::FACE,false,false));
    Teuchos::RCP<AmanziMesh::Mesh> surface_mesh =
      Teuchos::rcp(new AmanziMesh::Mesh_MSTK(*mesh,setnames,AmanziMesh::FACE,true,false));

    // push the mesh into state
    S->RegisterMesh("surface", surface_mesh);
    S->RegisterMesh("surface_3d", surface_mesh_3d);
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

  Teuchos::TimeMonitor::summarize();
  Teuchos::TimeMonitor::zeroOutTimers();
};


// -----------------------------------------------------------------------------
// Update any secondary (dependent) variables given a solution.
//
//   After a timestep is evaluated (or at ICs), there is no way of knowing if
//   secondary variables have been updated to be consistent with the new
//   solution.
// -----------------------------------------------------------------------------
void OverlandFlow::commit_state(double dt, const Teuchos::RCP<State>& S) {
  // Update flux if rel perm or h + Z has changed.
  bool update = UpdatePermeabilityData_(S.ptr());
  update |= S->GetFieldEvaluator("pres_elev")->HasFieldChanged(S.ptr(), name_);

  if (update) {
    // update the stiffness matrix with the new rel perm
    Teuchos::RCP<const CompositeVector> conductivity =
        S->GetFieldData("upwind_overland_conductivity");
    matrix_->CreateMFDstiffnessMatrices(conductivity.ptr());

    // derive the fluxes
    Teuchos::RCP<const CompositeVector> potential = S->GetFieldData("pres_elev");
    Teuchos::RCP<CompositeVector> flux = S->GetFieldData("overland_flux", name_);
    matrix_->DeriveFlux(*potential, flux);
  }
};


// -----------------------------------------------------------------------------
// Update diagnostics -- used prior to vis.
// -----------------------------------------------------------------------------
void OverlandFlow::calculate_diagnostics(const Teuchos::RCP<State>& S) {
  // update the cell velocities
  Teuchos::RCP<CompositeVector> velocity = S->GetFieldData("overland_velocity", name_);
  Teuchos::RCP<const CompositeVector> flux = S->GetFieldData("overland_flux");
  Teuchos::RCP<const CompositeVector> pressure = S->GetFieldData(key_);
  matrix_->DeriveCellVelocity(*flux, velocity);

  int ncells = velocity->size("cell");
  for (int c=0; c!=ncells; ++c) {
    (*velocity)("cell",0,c) /= std::max( (*pressure)("cell",c) , 1e-7);
    (*velocity)("cell",1,c) /= std::max( (*pressure)("cell",c) , 1e-7);
  }
};


// -----------------------------------------------------------------------------
// Use the physical rel perm (on cells) to update a work vector for rel perm.
//
//   This deals with upwinding, etc.
// -----------------------------------------------------------------------------
bool OverlandFlow::UpdatePermeabilityData_(const Teuchos::Ptr<State>& S) {
  bool update_perm = S->GetFieldEvaluator("overland_conductivity")
      ->HasFieldChanged(S, name_);
  update_perm |= S->GetFieldEvaluator("pres_elev")->HasFieldChanged(S, name_);

  if (update_perm) {
    // Update the perm only if needed.

    // This needs fixed to use the model, not assume a model.
    // Then it needs to be fixed to use a smart evaluator which picks
    // vals from cells and faces.
    // Then it needs to be fixed to work on a CompositeVector whose
    // components are cells and boundary faces. --etc
    Teuchos::RCP<const CompositeVector> pressure =
        S->GetFieldData(key_);
    Teuchos::RCP<const CompositeVector> slope =
        S->GetFieldData("slope_magnitude");
    Teuchos::RCP<const CompositeVector> manning =
        S->GetFieldData("manning_coefficient");
    Teuchos::RCP<CompositeVector> upwind_conductivity =
        S->GetFieldData("upwind_overland_conductivity", name_);

    // initialize the face coefficients
    upwind_conductivity->ViewComponent("face",true)->PutScalar(0.0);
    if (upwind_conductivity->has_component("cell")) {
      upwind_conductivity->ViewComponent("cell",true)->PutScalar(1.0);
    }

    // First do the boundary
    AmanziMesh::Entity_ID_List cells;

    int nfaces = upwind_conductivity->size("face", true);
    for (int f=0; f!=nfaces; ++f) {
      if (bc_markers_[f] != Operators::MFD_BC_NULL) {
        upwind_conductivity->mesh()->face_get_cells(f, AmanziMesh::USED, &cells);
        int c = cells[0];
        double scaling = (*manning)("cell",c) *
            std::sqrt(std::max((*slope)("cell",c), slope_regularization_));
        (*upwind_conductivity)("face",f) = std::pow(std::abs((*pressure)("face",f)), manning_exp_ + 1.0) / scaling ;
      }
    }

    // Then upwind.  This overwrites the boundary if upwinding says so.
    upwinding_->Update(S);

    // Communicate.  This could be done later, but i'm not exactly sure where, so
    // we'll do it here.
    upwind_conductivity->ScatterMasterToGhosted("face");
  }

  return update_perm;
}


// -----------------------------------------------------------------------------
// Evaluate boundary conditions at the current time.
// -----------------------------------------------------------------------------
void OverlandFlow::UpdateBoundaryConditions_(const Teuchos::Ptr<State>& S) {

  Teuchos::RCP<const CompositeVector> elevation = S->GetFieldData("elevation");
  Teuchos::RCP<const CompositeVector> pres = S->GetFieldData(key_);

  for (int n=0; n!=bc_markers_.size(); ++n) {
    bc_markers_[n] = Operators::MFD_BC_NULL;
    bc_values_[n] = 0.0;
  }

  Functions::BoundaryFunction::Iterator bc;
  int count = 0;
  for (bc=bc_pressure_->begin(); bc!=bc_pressure_->end(); ++bc) {
    int f = bc->first;
    bc_markers_[f] = Operators::MFD_BC_DIRICHLET;
    bc_values_[f] = bc->second + (*elevation)("face",f);
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
    bc_values_[f] = (*pres)("cell",cells[0]) + (*elevation)("face",f);
  }

  for (bc=bc_flux_->begin(); bc!=bc_flux_->end(); ++bc) {
    int f = bc->first;
    bc_markers_[f] = Operators::MFD_BC_FLUX;
    bc_values_[f] = bc->second;
  }
};

void OverlandFlow::UpdateBoundaryConditionsNoElev_(const Teuchos::Ptr<State>& S) {

  Teuchos::RCP<const CompositeVector> elevation = S->GetFieldData("elevation");
  Teuchos::RCP<const CompositeVector> pres = S->GetFieldData(key_);

  for (int n=0; n!=bc_markers_.size(); ++n) {
    bc_markers_[n] = Operators::MFD_BC_NULL;
    bc_values_[n] = 0.0;
  }

  Functions::BoundaryFunction::Iterator bc;
  int count = 0;
  for (bc=bc_pressure_->begin(); bc!=bc_pressure_->end(); ++bc) {
    int f = bc->first;
    bc_markers_[f] = Operators::MFD_BC_DIRICHLET;
    bc_values_[f] = bc->second;// + (*elevation)("face",f);
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
    bc_values_[f] = (*pres)("cell",cells[0]); // + (*elevation)("face",f);
  }

  for (bc=bc_flux_->begin(); bc!=bc_flux_->end(); ++bc) {
    int f = bc->first;
    bc_markers_[f] = Operators::MFD_BC_FLUX;
    bc_values_[f] = bc->second;
  }

  // Attempt of a hack to deal with zero rel perm
  double eps = 1.e-16;
  Teuchos::RCP<const CompositeVector> relperm =
      S->GetFieldData("upwind_overland_conductivity");
  for (int f=0; f!=relperm->size("face"); ++f) {
    if ((*relperm)("face",f) < eps) {
      bc_markers_[f] = Operators::MFD_BC_DIRICHLET;
      bc_values_[f] = 0.0;
    }
  }
};


/* ******************************************************************
 * Add a boundary marker to owned faces.
 ****************************************************************** */
void OverlandFlow::ApplyBoundaryConditions_(const Teuchos::RCP<State>& S,
        const Teuchos::RCP<CompositeVector>& pres) {
  int nfaces = pres->size("face",true);
  for (int f=0; f!=nfaces; ++f) {
    if (bc_markers_[f] == Operators::MFD_BC_DIRICHLET) {
      (*pres)("face",f) = bc_values_[f];
    }
  }
};


} // namespace
} // namespace
