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
#include "independent_field_model.hh"

#include "upwind_potential_difference.hh"
#include "elevation_model.hh"
#include "meshed_elevation_model.hh"
#include "standalone_elevation_model.hh"
#include "manning_conductivity_model.hh"

#include "overland.hh"

namespace Amanzi {
namespace Flow {

#if 0 // a trick for xemacs indent algorithm
}}
#endif

RegisteredPKFactory<OverlandFlow> OverlandFlow::reg_("overland flow");

// -------------------------------------------------------------
// Constructor
// -------------------------------------------------------------
OverlandFlow::OverlandFlow(Teuchos::ParameterList& flow_plist,
                           const Teuchos::RCP<State>& S,
                           const Teuchos::RCP<TreeVector>& solution) :
  flow_plist_(flow_plist), flow_time_(0.), nsteps_(0) {

  solution_ = solution;
  CreateMesh_(S);

  // Require fields and models for those fields.
  // -- primary variable: pressure on both cells and faces, ghosted, with 1 dof
  std::vector<AmanziMesh::Entity_kind> locations2(2);
  std::vector<std::string> names2(2);
  std::vector<int> num_dofs2(2,1);
  locations2[0] = AmanziMesh::CELL;
  locations2[1] = AmanziMesh::FACE;
  names2[0] = "cell";
  names2[1] = "face";

  S->RequireField("overland_pressure", "overland_flow")->SetMesh(S->Mesh("surface"))
                ->SetGhosted()->SetComponents(names2, locations2, num_dofs2);
  Teuchos::RCP<PrimaryVariableFieldModel> pressure_model =
      Teuchos::rcp(new PrimaryVariableFieldModel("overland_pressure"));
  S->SetFieldModel("overland_pressure", pressure_model);

  // -- secondary variables model for surface geometry.
  S->RequireField("elevation", "overland_flow")->SetMesh(S->Mesh("surface"))->SetGhosted()
                ->SetComponents(names2, locations2, num_dofs2);
  S->RequireField("slope_magnitude", "overland_flow")->SetMesh(S->Mesh("surface"))
                ->SetGhosted()->SetComponent("cell", AmanziMesh::CELL, 1);
  S->RequireField("pres_elev", "overland_flow")->SetMesh(S->Mesh("surface"))->SetGhosted()
                ->SetComponents(names2, locations2, num_dofs2);

  Teuchos::RCP<FlowRelations::ElevationModel> elev_model;
  if (standalone_mode_) {
    Teuchos::ParameterList elev_plist = flow_plist_.sublist("elevation model");
    elev_model = Teuchos::rcp(new FlowRelations::StandaloneElevationModel(elev_plist));
  } else {
    elev_model = Teuchos::rcp(new FlowRelations::MeshedElevationModel());
  }
  S->SetFieldModel("elevation", elev_model);
  S->SetFieldModel("slope_magnitude", elev_model);
  S->SetFieldModel("pres_elev", elev_model);

  // -- secondary variables, no model used
  S->RequireField("overland_flux", "overland_flow")->SetMesh(S->Mesh("surface"))
                ->SetGhosted()->SetComponent("face", AmanziMesh::FACE, 1);
  S->RequireField("overland_velocity", "overland_flow")->SetMesh(S->Mesh("surface"))
                ->SetGhosted()->SetComponent("cell", AmanziMesh::CELL, 3);

  // -- "rel perm" model
  S->RequireField("overland_conductivity", "overland_flow")->SetMesh(S->Mesh("surface"))
                ->SetGhosted()->SetComponent("cell", AmanziMesh::CELL, 1);
  Teuchos::ParameterList cond_plist = flow_plist_.sublist("overland conductivity model");
  Teuchos::RCP<FieldModel> cond_model =
      Teuchos::rcp(new FlowRelations::ManningConductivityModel(cond_plist));
  S->SetFieldModel("overland_conductivity", cond_model);

  // -- source term
  if (flow_plist_.isSublist("source model")) {
    Teuchos::ParameterList source_plist = flow_plist_.sublist("source model");
    source_plist.set("model name", "source model");
    is_source_term_ = true;
    S->RequireField("overland_source", "overland_flow")->SetMesh(S->Mesh("surface"))
        ->SetGhosted()->SetComponent("cell", AmanziMesh::CELL, 1);
    Teuchos::RCP<FieldModel> source_model =
        Teuchos::rcp(new IndependentFieldModel(source_plist));
    S->SetFieldModel("overland_source", source_model);
  }

  // -- coupling term
  if (flow_plist_.get<bool>("coupled to subsurface", false)) {
    is_coupling_term_ = true;
    S->RequireField("overland_source_from_subsurface", "overland_flow")
        ->SetMesh(S->Mesh("surface"))->SetGhosted()
        ->SetComponent("cell", AmanziMesh::CELL, 1);
    S->RequireFieldModel("overland_coupling");
  }


  // -- cell volume and model
  S->RequireField("surface_cell_volume")->SetMesh(S->Mesh("surface"))->SetGhosted()
                                ->AddComponent("cell", AmanziMesh::CELL, 1);
  S->RequireFieldModel("surface_cell_volume");


  // boundary conditions
  Teuchos::ParameterList bc_plist = flow_plist_.sublist("boundary conditions", true);
  FlowBCFactory bc_factory(S->Mesh("surface"), bc_plist);
  bc_pressure_ = bc_factory.CreatePressure();
  bc_zero_gradient_ = bc_factory.CreateZeroGradient();
  bc_flux_ = bc_factory.CreateMassFlux();

  // rel perm upwinding
  S->RequireField("upwind_overland_conductivity", "overland_flow")
                ->SetMesh(S->Mesh("surface"))->SetGhosted()
                ->SetComponents(names2, locations2, num_dofs2);
  S->GetField("upwind_overland_conductivity","overland_flow")->set_io_vis(false);

  upwinding_ = Teuchos::rcp(new Operators::UpwindPotentialDifference("overland_flow",
          "overland_conductivity", "upwind_overland_conductivity",
          "pres_elev", "overland_pressure"));

  // operator for the diffusion terms
  Teuchos::ParameterList mfd_plist = flow_plist_.sublist("Diffusion");
  matrix_ = Teuchos::rcp(new Operators::MatrixMFD(mfd_plist, S->Mesh("surface")));
  bool symmetric = false;
  matrix_->SetSymmetryProperty(symmetric);
  matrix_->SymbolicAssembleGlobalMatrices();

  // preconditioner for the NKA system
  Teuchos::ParameterList mfd_pc_plist = flow_plist_.sublist("Diffusion PC");
  preconditioner_ = Teuchos::rcp(new Operators::MatrixMFD(mfd_pc_plist, S->Mesh("surface")));
  preconditioner_->SetSymmetryProperty(symmetric);
  preconditioner_->SymbolicAssembleGlobalMatrices();
  preconditioner_->InitPreconditioner(mfd_pc_plist);

  steptime_ = Teuchos::TimeMonitor::getNewCounter("overland flow advance time");

};


// -------------------------------------------------------------
// Initialize PK
// -------------------------------------------------------------
void OverlandFlow::initialize(const Teuchos::RCP<State>& S) {
  // initial timestep size
  dt_ = flow_plist_.get<double>("Initial time step", 1.);

  // initialize boundary conditions
  int nfaces = S->Mesh("surface")->num_entities(AmanziMesh::FACE, AmanziMesh::USED);
  bc_markers_.resize(nfaces, Operators::MFD_BC_NULL);
  bc_values_.resize(nfaces, 0.0);

  bc_pressure_->Compute(S->time());
  bc_zero_gradient_->Compute(S->time());
  bc_flux_->Compute(S->time());
  UpdateBoundaryConditions_(S);

  // initial conditions
  // -- Get the IC function plist.
  if (!flow_plist_.isSublist("initial condition")) {
    std::stringstream messagestream;
    messagestream << "Overland Flow PK has no initial condition parameter list.";
    Errors::Message message(messagestream.str());
    Exceptions::amanzi_throw(message);
  }

  // -- Calculate the IC.
  Teuchos::ParameterList ic_plist = flow_plist_.sublist("initial condition");
  Teuchos::RCP<CompositeVector> pres =
      S->GetFieldData("overland_pressure", "overland_flow");
  Teuchos::RCP<Functions::CompositeVectorFunction> ic =
      Functions::CreateCompositeVectorFunction(ic_plist, *pres);
  ic->Compute(S->time(), pres.ptr());

  // -- Initialize face values as the mean of neighboring cell values.
  DeriveFaceValuesFromCellValues_(S, pres);

  // Set extra fields as initialized -- these don't currently have models.
  S->GetFieldData("upwind_overland_conductivity","flow")->PutScalar(1.0);
  S->GetField("upwind_overland_conductivity","flow")->set_initialized();
  S->GetField("overland_flux", "flow")->set_initialized();
  S->GetField("overland_velocity", "flow")->set_initialized();

  // Initialize operators.
  matrix_->CreateMFDmassMatrices(Teuchos::null);
  preconditioner_->CreateMFDmassMatrices(Teuchos::null);

  // Initialize the timestepper.
  solution_->set_data(pres);
  atol_ = flow_plist_.get<double>("Absolute error tolerance",1.0);
  rtol_ = flow_plist_.get<double>("Relative error tolerance",1.0);
  precon_lag_ = flow_plist_.get<int>("preconditioner lag", 0);

  if (!flow_plist_.get<bool>("Strongly Coupled PK", false)) {
    // -- Instantiate the time stepper.
    Teuchos::RCP<Teuchos::ParameterList> bdf1_plist_p =
      Teuchos::rcp(new Teuchos::ParameterList(flow_plist_.sublist("Time integrator")));
    time_stepper_ = Teuchos::rcp(new BDF1TimeIntegrator(this, bdf1_plist_p, solution_));
    time_step_reduction_factor_ = bdf1_plist_p->get<double>("time step reduction factor");

    // -- Initialize the time derivative.
    Teuchos::RCP<TreeVector> solution_dot = Teuchos::rcp( new TreeVector(*solution_));
    solution_dot->PutScalar(0.0);

    // -- Set the initial state.
    time_stepper_->set_initial_state(S->time(), solution_, solution_dot);
  }
};


void OverlandFlow::CreateMesh_(const Teuchos::RCP<State>& S) {
  // Create mesh
  if (S->Mesh()->space_dimension() == 3) {

    // The domain mesh is a 3D volume mesh or a 3D surface mesh -- construct
    // the surface mesh.

    // -- Ensure the domain mesh is MSTK.
    Teuchos::RCP<const AmanziMesh::Mesh_MSTK> mesh =
      Teuchos::rcp_dynamic_cast<const AmanziMesh::Mesh_MSTK>(S->Mesh());

    if (mesh == Teuchos::null) {
      Errors::Message message("Overland Flow PK requires surface mesh, which is currently only supported by MSTK.  Make the domain mesh an MSTK mesh.");
      Exceptions::amanzi_throw(message);
    }
    Teuchos::RCP<Teuchos::Time> meshtime = Teuchos::TimeMonitor::getNewCounter("surface mesh creation");
    Teuchos::TimeMonitor timer(*meshtime);

    // Check that the surface mesh has a subset
    std::string surface_sideset_name =
      flow_plist_.get<std::string>("Surface sideset name");

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
    standalone_mode_ = false;

  } else if (S->Mesh()->space_dimension() == 2) {
    // The domain mesh is already a 2D mesh, so simply register it as the surface
    // mesh as well.
    S->RegisterMesh("surface", S->Mesh());
    standalone_mode_ = true;
  } else {
    Errors::Message message("Invalid mesh dimension for overland flow.");
    Exceptions::amanzi_throw(message);
  }

  Teuchos::TimeMonitor::summarize();
  Teuchos::TimeMonitor::zeroOutTimers();
};

// -----------------------------------------------------------------------------
// Transfer operators -- ONLY COPIES POINTERS
// -----------------------------------------------------------------------------
void OverlandFlow::state_to_solution(const Teuchos::RCP<State>& S,
                                     const Teuchos::RCP<TreeVector>& solution) {
  solution->set_data(S->GetFieldData("overland_pressure", "overland_flow"));
};


// -----------------------------------------------------------------------------
// Transfer operators -- ONLY COPIES POINTERS
// -----------------------------------------------------------------------------
void OverlandFlow::solution_to_state(const Teuchos::RCP<TreeVector>& solution,
                                     const Teuchos::RCP<State>& S) {
  S->SetData("overland_pressure", "overland_flow", solution->data());
};


// -----------------------------------------------------------------------------
// Advance from state S to state S_next at time S.time + dt.
// -----------------------------------------------------------------------------
bool OverlandFlow::advance(double dt) {
  state_to_solution(S_next_, solution_);

  // save flow_rate
  // 2D test case
  // if ( nsteps_%6==0 || true ) {
  //   output_flow_rate() ;
  // }
  ++nsteps_ ;

  // take a bdf timestep
  double h = dt;

  niter_ = 0;
  ntries_++;

  double dt_solver;
  try {
    Teuchos::TimeMonitor timer(*steptime_);
    //time_stepper_->time_step(h, solution_);
    dt_solver = time_stepper_->time_step(h, solution_);
    flow_time_ += h;
  } catch (Exceptions::Amanzi_exception &error) {
    if (S_next_->Mesh("surface")->get_comm()->MyPID() == 0) {
      std::cout << "Timestepper called error: " << error.what() << std::endl;
    }
    if (error.what() == std::string("BDF time step failed") ||
        error.what() == std::string("Cut time step")) {
      // try cutting the timestep
      dt_ = h*time_step_reduction_factor_;
      return true;
    } else {
      throw error;
    }
  }

  ntries_ = 0;
  Teuchos::TimeMonitor::summarize();
  Teuchos::TimeMonitor::zeroOutTimers();

  // commit the step as successful
  time_stepper_->commit_solution(h, solution_);
  solution_to_state(solution_, S_next_);
  commit_state(h, S_next_);

  // update the timestep size
  if (dt_solver < dt_ && dt_solver >= h) {
    // We took a smaller step than we recommended, and it worked fine (not
    // suprisingly).  Likely this was due to constraints from other PKs or
    // vis.  Do not reduce our recommendation.
  } else {
    dt_ = dt_solver;
  }

  return false;
};


// -----------------------------------------------------------------------------
// Update any secondary (dependent) variables given a solution.
//
//   After a timestep is evaluated (or at ICs), there is no way of knowing if
//   secondary variables have been updated to be consistent with the new
//   solution.
// -----------------------------------------------------------------------------
void OverlandFlow::commit_state(double dt, const Teuchos::RCP<State>& S) {
  // update secondary variables for IC

  // update the flux
  UpdatePermeabilityData_(S);
  Teuchos::RCP<const CompositeVector> upwind_conductivity =
    S->GetFieldData("upwind_overland_conductivity");

  Teuchos::RCP<const CompositeVector> pres =
    S->GetFieldData("overland_pressure", "overland_flow");
  Teuchos::RCP<const CompositeVector> pres_elev =
    S->GetFieldData("pres_elev");
  Teuchos::RCP<const CompositeVector> elevation =
    S->GetFieldData("elevation");

  Teuchos::RCP<CompositeVector> darcy_flux =
    S->GetFieldData("overland_flux", "overland_flow");

  // communicate face perms
  matrix_->CreateMFDstiffnessMatrices(*upwind_conductivity);
  matrix_->DeriveFlux(*pres_elev, darcy_flux);
};


// -----------------------------------------------------------------------------
// Update diagnostics -- used prior to vis.
// -----------------------------------------------------------------------------
void OverlandFlow::calculate_diagnostics(const Teuchos::RCP<State>& S) {
  // update the cell velocities
  Teuchos::RCP<CompositeVector> velocity = S->GetFieldData("overland_velocity", "overland_flow");
  Teuchos::RCP<const CompositeVector> flux = S->GetFieldData("overland_flux");
  Teuchos::RCP<const CompositeVector> pressure = S->GetFieldData("overland_pressure");
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
void OverlandFlow::UpdatePermeabilityData_(const Teuchos::RCP<State>& S) {
  // This needs fixed to use the model, not assume a model. -- etc
  Teuchos::RCP<const CompositeVector> pressure =
    S->GetFieldData("overland_pressure");
  Teuchos::RCP<const CompositeVector> slope =
    S->GetFieldData("slope_magnitude");
  Teuchos::RCP<const CompositeVector> manning =
    S->GetFieldData("manning_coef");
  Teuchos::RCP<CompositeVector> upwind_conductivity =
    S->GetFieldData("upwind_overland_conductivity", "overland_flow");

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
      //std::sqrt((*slope)("cell",c) + slope_regularization_);
      (*upwind_conductivity)("face",f) = std::pow(std::abs((*pressure)("face",f)), manning_exp_ + 1.0) / scaling ;
      //(*upwind_conductivity)("face",f) = std::pow((*pressure)("face",f), manning_exp_ + 1.0) / scaling ;
    }
  }

  // Then upwind.  This overwrites the boundary if upwinding says so.
  upwinding_->Update(S.ptr());

  // Communicate.  This could be done later, but i'm not exactly sure where, so
  // we'll do it here.
  upwind_conductivity->ScatterMasterToGhosted("face");

}


// -----------------------------------------------------------------------------
// Interpolate pressure ICs on cells to ICs for lambda (faces).
// -----------------------------------------------------------------------------
void OverlandFlow::DeriveFaceValuesFromCellValues_(const Teuchos::RCP<State>& S,
        const Teuchos::RCP<CompositeVector> & pres) {
  AmanziMesh::Entity_ID_List cells;
  pres->ScatterMasterToGhosted("cell");

  int f_owned = pres->size("face");
  for (int f=0; f!=f_owned; ++f) {
    cells.clear();
    S->Mesh("surface")->face_get_cells(f, AmanziMesh::USED, &cells);

    int ncells = cells.size();
    double face_value = 0.0;
    for (int n=0; n!=ncells; ++n) {
      face_value += (*pres)("cell",0,cells[n]);
    }
    (*pres)("face",f) = face_value / ncells;
  }
};


// -----------------------------------------------------------------------------
// Evaluate boundary conditions at the current time.
// -----------------------------------------------------------------------------
void OverlandFlow::UpdateBoundaryConditions_(const Teuchos::RCP<State>& S) {

  Teuchos::RCP<const CompositeVector> elevation = S->GetFieldData("elevation");
  Teuchos::RCP<const CompositeVector> pres = S->GetFieldData("overland_pressure");

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
    S->Mesh("surface")->face_get_cells(f, AmanziMesh::USED, &cells);
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

void OverlandFlow::UpdateBoundaryConditionsNoElev_(const Teuchos::RCP<State>& S) {

  Teuchos::RCP<const CompositeVector> elevation = S->GetFieldData("elevation");
  Teuchos::RCP<const CompositeVector> pres = S->GetFieldData("overland_pressure");

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
    S->Mesh("surface")->face_get_cells(f, AmanziMesh::USED, &cells);
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
