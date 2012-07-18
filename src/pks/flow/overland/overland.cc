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
#include "composite_vector_function.hh"
#include "composite_vector_function_factory.hh"
#include "upwind_total_flux.hh"
#include "Point.hh"

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
  flow_plist_(flow_plist), flow_time_(0.), niter_(0) {

  solution_ = solution;

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


  // Create data structures
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

  // -- secondary variable: elevation on both cells and faces, ghosted, with 1 dof
  S->RequireField("elevation", "overland_flow")->SetMesh(S->Mesh("surface"))->SetGhosted()
                ->SetComponents(names2, locations2, num_dofs2);
  S->RequireField("slope_magnitude", "overland_flow")->SetMesh(S->Mesh("surface"))
                ->SetGhosted()->SetComponent("cell", AmanziMesh::CELL, 1);

  // -- secondary variable: pres_elev on both cells and faces, ghosted, with 1 dof
  S->RequireField("pres_elev", "overland_flow")->SetMesh(S->Mesh("surface"))->SetGhosted()
                ->SetComponents(names2, locations2, num_dofs2);

  // -- other secondary variables
  S->RequireField("overland_flux", "overland_flow")->SetMesh(S->Mesh("surface"))
                ->SetGhosted()->SetComponent("face", AmanziMesh::FACE, 1);
  S->RequireField("overland_velocity", "overland_flow")->SetMesh(S->Mesh("surface"))
                ->SetGhosted()->SetComponent("cell", AmanziMesh::CELL, 3);
                                  // NOTE this is 3 because VisIt is dumb.
  S->RequireField("overland_conductivity", "overland_flow")->SetMesh(S->Mesh("surface"))
                ->SetGhosted()->SetComponent("cell", AmanziMesh::CELL, 1);
  S->RequireField("rainfall_rate", "overland_flow")->SetMesh(S->Mesh("surface"))
                ->SetGhosted()->SetComponent("cell", AmanziMesh::CELL, 1);
  S->RequireField("manning_coef", "overland_flow")->SetMesh(S->Mesh("surface"))
                ->SetGhosted()->SetComponent("cell", AmanziMesh::CELL, 1);

  // -- independent variables not owned by this PK
  S->RequireField("surface_cell_volume")->SetMesh(S->Mesh("surface"))->SetGhosted()
                                ->AddComponent("cell", AmanziMesh::CELL, 1);

  // -- work vectors
  S->RequireField("upwind_overland_conductivity", "overland_flow")
                ->SetMesh(S->Mesh("surface"))->SetGhosted()
                ->SetComponents(names2, locations2, num_dofs2);
  S->GetField("upwind_overland_conductivity","overland_flow")->set_io_vis(false);

  // boundary conditions
  Teuchos::ParameterList bc_plist = flow_plist_.sublist("boundary conditions", true);
  FlowBCFactory bc_factory(S->Mesh("surface"), bc_plist);
  bc_pressure_ = bc_factory.CreatePressure();
  bc_zero_gradient_ = bc_factory.CreateZeroGradient();
  bc_flux_ = bc_factory.CreateMassFlux();

  // rel perm upwinding
  upwinding_ = Teuchos::rcp(new Operators::UpwindTotalFlux("overland_flow",
          "overland_conductivity", "upwind_overland_conductivity",
          "overland_flux"));

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
  Teuchos::ParameterList mfd_pc_ml_plist = mfd_pc_plist.sublist("ML Parameters");
  preconditioner_->InitMLPreconditioner(mfd_pc_ml_plist);
};


// -------------------------------------------------------------
// Initialize PK
// -------------------------------------------------------------
void OverlandFlow::initialize(const Teuchos::RCP<State>& S) {
  // initial timestep size
  dt_ = flow_plist_.get<double>("Initial time step", 1.);

  // initialize boundary conditions
  int nfaces = S->Mesh("surface")->num_entities(AmanziMesh::FACE, AmanziMesh::USED);
  int nfaces_owned = S->Mesh("surface")->num_entities(AmanziMesh::FACE, AmanziMesh::OWNED);

  int ncells = S->Mesh("surface")->num_entities(AmanziMesh::CELL, AmanziMesh::USED);
  int ncells_owned = S->Mesh("surface")->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);

  int nnodes = S->Mesh("surface")->num_entities(AmanziMesh::NODE, AmanziMesh::USED);
  int nnodes_owned = S->Mesh("surface")->num_entities(AmanziMesh::NODE, AmanziMesh::OWNED);

  std::cout << "Overland Initialize: " << S->Mesh("surface")->get_comm()->MyPID() << std::endl;
  std::cout << "  nnodes (owned, used):" << nnodes_owned << " " << nnodes << std::endl;
  std::cout << "  nfaces (owned, used):" << nfaces_owned << " " << nfaces << std::endl;
  std::cout << "  ncells (owned, used):" << ncells_owned << " " << ncells << std::endl;

  int dnfaces = S->Mesh("domain")->num_entities(AmanziMesh::FACE, AmanziMesh::USED);
  int dnfaces_owned = S->Mesh("domain")->num_entities(AmanziMesh::FACE, AmanziMesh::OWNED);

  int dncells = S->Mesh("domain")->num_entities(AmanziMesh::CELL, AmanziMesh::USED);
  int dncells_owned = S->Mesh("domain")->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);

  int dnnodes = S->Mesh("domain")->num_entities(AmanziMesh::NODE, AmanziMesh::USED);
  int dnnodes_owned = S->Mesh("domain")->num_entities(AmanziMesh::NODE, AmanziMesh::OWNED);

  std::cout << "Domain Initialize: " << S->Mesh("domain")->get_comm()->MyPID() << std::endl;
  std::cout << "  nnodes (owned, used):" << dnnodes_owned << " " << dnnodes << std::endl;
  std::cout << "  nfaces (owned, used):" << dnfaces_owned << " " << dnfaces << std::endl;
  std::cout << "  ncells (owned, used):" << dncells_owned << " " << dncells << std::endl;

  bc_markers_.resize(nfaces, Operators::MFD_BC_NULL);
  bc_values_.resize(nfaces, 0.0);

  // update face pressures as a hint?
  Teuchos::RCP<CompositeVector> pres =
    S->GetFieldData("overland_pressure", "overland_flow");
  DeriveFaceValuesFromCellValues_(S, pres ) ;

  // Initialize the elevation field and its slope.
  if (standalone_mode_) {
    // -- Elevation field is not already known -- get a function to specify it.
    Teuchos::RCP<CompositeVector> elev = S->GetFieldData("elevation", "overland_flow");
    Teuchos::ParameterList elev_plist = flow_plist_.sublist("Elevation model");
    elevation_function_ =
      Functions::CreateCompositeVectorFunction(elev_plist, elev.ptr());

    // -- Slope field is not already known -- get a function to specify it.
    Teuchos::RCP<CompositeVector> slope = S->GetFieldData("slope_magnitude", "overland_flow");
    Teuchos::ParameterList slope_plist = flow_plist_.sublist("Slope model");
    slope_function_ =
      Functions::CreateCompositeVectorFunction(slope_plist, slope.ptr());
  }

  // -- Evaluate the functions or use the mesh to get the values.
  UpdateElevationAndSlope_(S);
  S->GetField("elevation","overland_flow")->set_initialized();
  S->GetField("slope_magnitude","overland_flow")->set_initialized();

  // Initialize the rainfall model.
  Teuchos::RCP<const CompositeVector> rain = S->GetFieldData("rainfall_rate");
  Teuchos::ParameterList rain_plist = flow_plist_.sublist("Rainfall model");
  rain_rate_function_ =
    Functions::CreateCompositeVectorFunction(rain_plist, rain.ptr());

  // Initialize Manning Coefficient.
  Teuchos::RCP<CompositeVector> mann = S->GetFieldData("manning_coef", "overland_flow");
  Teuchos::ParameterList mann_plist = flow_plist_.sublist("Manning coefficient");
  manning_coef_function_ =
    Functions::CreateCompositeVectorFunction(mann_plist, mann.ptr());
  manning_coef_function_->Compute(S->time(), mann.ptr());
  S->GetField("manning_coef","overland_flow")->set_initialized();

  // Get the Manning exponent.
  manning_exp_ = flow_plist_.get<double>("Manning exponent");

  // Declare secondary variables initialized, as they will be done by
  // the commit_state call.
  S->GetField("overland_conductivity", "overland_flow")->set_initialized();
  S->GetField("overland_flux", "overland_flow")->set_initialized();
  S->GetField("pres_elev", "overland_flow")->set_initialized();
  S->GetField("rainfall_rate", "overland_flow")->set_initialized();

  // velocity needs component names to get it into a vector
  Teuchos::RCP<FieldCompositeVector> velocity =
    Teuchos::rcp_static_cast<FieldCompositeVector>(S->GetField("overland_velocity", "overland_flow"));
  std::vector<std::vector<std::string> > names(1);
  names[0].resize(3);
  names[0][0] = "X-Component";
  names[0][1] = "Y-Component";
  names[0][2] = "Z-Component";
  velocity->set_subfield_names(names);
  velocity->set_initialized();
  S->GetFieldData("overland_velocity","overland_flow")->PutScalar(1.0);

  // Rel perm is special -- if the mode is symmetric, it needs to be
  // initialized to 1.
  S->GetFieldData("upwind_overland_conductivity","overland_flow")->PutScalar(1.0);
  S->GetField("upwind_overland_conductivity", "overland_flow")->set_initialized();

  // Initialize BCs.
  bc_pressure_->Compute(S->time());
  bc_zero_gradient_->Compute(S->time());
  bc_flux_->Compute(S->time());
  UpdateBoundaryConditions_(S);

  // Initialize operators.
  matrix_->CreateMFDmassMatrices(Teuchos::null);
  preconditioner_->CreateMFDmassMatrices(Teuchos::null);

  // Initialize the timestepper.
  solution_->set_data(pres);
  atol_ = flow_plist_.get<double>("Absolute error tolerance",1e-6);
  rtol_ = flow_plist_.get<double>("Relative error tolerance",1e0);

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
  if ( niter_%6==0 || true ) {
    output_flow_rate() ;
  }
  ++niter_ ;

  // take a bdf timestep
  double h = dt;

  //  try {
    time_stepper_->time_step(h, solution_);
    //dt_ = time_stepper_->time_step(h, solution_);
    flow_time_ += h;
  // } catch (Exceptions::Amanzi_exception &error) {
  //   std::cout << "Timestepper called error: " << error.what() << std::endl;
  //   if (error.what() == std::string("BDF time step failed")) {
  //     // try cutting the timestep
  //     dt_ = dt_*time_step_reduction_factor_;
  //     return true;
  //   } else {
  //     throw error;
  //   }
  // }

  // commit the step as successful
  time_stepper_->commit_solution(h, solution_);
  solution_to_state(solution_, S_next_);
  commit_state(h, S_next_);

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
  UpdateSecondaryVariables_(S);

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
    (*velocity)("cell",0,c) /= (*pressure)("cell",c);
    (*velocity)("cell",1,c) /= (*pressure)("cell",c);
  }
};


// -----------------------------------------------------------------------------
// Use the physical rel perm (on cells) to update a work vector for rel perm.
//
//   This deals with upwinding, etc.
// -----------------------------------------------------------------------------
void OverlandFlow::UpdatePermeabilityData_(const Teuchos::RCP<State>& S) {
  upwinding_->Update(S.ptr());

  // patch up at the boundary condition
  Teuchos::RCP<const CompositeVector> pressure =
    S->GetFieldData("overland_pressure");
  Teuchos::RCP<const CompositeVector> slope =
    S->GetFieldData("slope_magnitude");
  Teuchos::RCP<const CompositeVector> manning =
    S->GetFieldData("manning_coef");
  Teuchos::RCP<CompositeVector> upwind_conductivity =
    S->GetFieldData("upwind_overland_conductivity", "overland_flow");

  AmanziMesh::Entity_ID_List cells;
  double eps = 1.e-14;

  int nfaces = upwind_conductivity->size("face", true);
  for (int f=0; f!=nfaces; ++f) {
    if (bc_markers_[f] != Operators::MFD_BC_NULL) {
      upwind_conductivity->mesh()->face_get_cells(f, AmanziMesh::USED, &cells);
      int c = cells[0];
      double scaling = (*manning)("cell",c) * std::sqrt((*slope)("cell",c) + eps);
      (*upwind_conductivity)("face",f) = std::pow( (*pressure)("face",f), manning_exp_ + 1.0) / scaling ;
    }
  }
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
  for (bc=bc_pressure_->begin(); bc!=bc_pressure_->end(); ++bc) {
    int f = bc->first;
    bc_markers_[f] = Operators::MFD_BC_DIRICHLET;
    bc_values_[f] = bc->second + (*elevation)("face",f);
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

/* ******************************************************************
 * ELEVATION STUFF
 ****************************************************************** */

// Test 1
// single zone model
void OverlandFlow::TestOneSetUpElevationPars_() {
  zone_x.resize(0) ;
  // ---
  slope_x.resize(1) ;
  slope_x[0] = 0.0016 ;
  // ---
  slope_y.resize(1) ;
  slope_y[0] = 0.   ;
  // ---
  manning.resize(1) ;
  manning[0] = 0.025 ;
}
double OverlandFlow::TestOneElevation_( double x, double y ) {
  double Lx = 600./3.28084 ;
  return slope_x[0]*( Lx - x ) ;
}

// Test 2
// three zones model
int OverlandFlow::TestTwoZoneFlag_( double x, double y ) { 
  int retval = -1 ; // zone unset by default
  if ( zone_x[0]<=x && x<  zone_x[1] ) { retval = 0 ; } 
  if ( zone_x[1]<=x && x<= zone_x[2] ) { retval = 1 ; } 
  if ( zone_x[2]<x  && x<= zone_x[3] ) { retval = 2 ; }
  assert( retval==0 || retval==1 || retval==2 ) ; // check proper assignement
  return retval ; 
}

// set elevation pars
void OverlandFlow::TestTwoSetUpElevationPars_() {
  double dx =  2.e-6 ;
  // ---
  zone_x.resize(4) ;
  zone_x[0] = 0.-dx ;
  zone_x[1] = 800.  ;
  zone_x[2] = 820.  ;
  zone_x[3] = 1620+dx ;
  // ---
  slope_x.resize(3) ;
  slope_x[0] = 0.05 ;
  slope_x[1] = 0.0  ;
  slope_x[2] = 0.05  ;
  // ---
  slope_y.resize(3) ;
  slope_y[0] = 0.02 ;
  slope_y[1] = 0.02 ;
  slope_y[2] = 0.02 ;
  // ---
  manning.resize(3) ;
  manning[0] = 0.015 ;
  manning[1] = 0.15  ;
  manning[2] = 0.015 ; 
}


void OverlandFlow::SetUpElevation_( const Teuchos::RCP<State>& S ) {

#if 0

  // TEST 1
  TestOneSetUpElevationPars_();

#else

  // TEST 2
  TestTwoSetUpElevationPars_();

#endif

}


// OTHER DEBUGGING STUFF

void OverlandFlow::print_pressure( const Teuchos::RCP<State>& S, string prt_str ) const {

  MSGF("begin  ---> "<<prt_str) ;

  int c_owned = S->Mesh("surface")->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  int f_owned = S->Mesh("surface")->num_entities(AmanziMesh::FACE, AmanziMesh::OWNED);

  const Teuchos::RCP<CompositeVector> pressure = 
    S->GetFieldData("overland_pressure", "overland_flow");

#if 0
  Teuchos::RCP<CompositeVector> elevation = 
    S->GetFieldData("elevation", "overland_flow");
  solution_->set_data(elevation);
  for (int c=0; c!=c_owned; ++c) {
    printf("cell elevation(%5i)=%14.7e\n",c,(*elevation)("cell",0,c));
  }
  LINE(--) ;
  for (int f=0; f!=f_owned; ++f) {
    printf("face elevation(%5i)=%14.7e\n",f,(*elevation)("face",0,f));
  }
#endif

#if 0
  Teuchos::RCP<CompositeVector> pres_elev = 
    S->GetFieldData("pres_elev", "overland_flow");
  solution_->set_data(pres_elev);
  for (int c=0; c!=c_owned; ++c) {
    printf("cell pres_elev(%5i)=%14.7e\n",c,(*pres_elev)("cell",0,c));
  }
  LINE(--) ;
  for (int f=0; f!=f_owned; ++f) {
    printf("face pres_elev(%5i)=%14.7e\n",f,(*pres_elev)("face",0,f));
  }
#endif

  MSGF("end of ---> "<<prt_str) ;
}


void OverlandFlow::print_vector( const Teuchos::RCP<State>& S, 
                                 const Teuchos::RCP<const CompositeVector>& pressure, 
                                 string prt_str ) const {

  MSGF("begin printing ---> "<<prt_str) ;

  int c_owned = S->Mesh("surface")->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  int f_owned = S->Mesh("surface")->num_entities(AmanziMesh::FACE, AmanziMesh::OWNED);

#if 1
  for (int c=0; c!=c_owned; ++c) {
    printf("cell value(%5i)=%14.7e\n",c,(*pressure)("cell",0,c));
  }
  LINE(--) ;
  for (int f=0; f!=f_owned; ++f) {
    printf("face value(%5i)=%14.7e\n",f,(*pressure)("face",0,f));
  }
#endif
  MSGF("end of printing ---> "<<prt_str) ;
}

void OverlandFlow::print_vector2( const Teuchos::RCP<State>& S, 
                                 const Teuchos::RCP<const CompositeVector>& pressure, 
                                 const Teuchos::RCP<const CompositeVector>& elevation, 
                                 string prt_str ) const {

  MSGF("begin printing ---> "<<prt_str) ;

  int c_owned = S->Mesh("surface")->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  int f_owned = S->Mesh("surface")->num_entities(AmanziMesh::FACE, AmanziMesh::OWNED);

#if 1
  for (int c=0; c!=c_owned; ++c) {
    printf("cell: pres(%5i)=%14.7e   elev(%5i)=%14.7e\n",c,(*pressure)("cell",0,c),c,(*elevation)("cell",0,c));
  }
  LINE(--) ;
  for (int f=0; f!=f_owned; ++f) {
    printf("face: pres(%5i)=%14.7e   elev(%5i)=%14.7e\n",f,(*pressure)("face",0,f),f,(*elevation)("face",0,f));
  }
#endif
  MSGF("end of printing ---> "<<prt_str) ;
}

void OverlandFlow::print_faceval( const Teuchos::RCP<State>& S, 
                                  const Teuchos::RCP<const CompositeVector>& vec, 
                                  string prt_str ) const {

  MSGF("begin printing ---> "<<prt_str) ;

  int f_owned = S->Mesh("surface")->num_entities(AmanziMesh::FACE, AmanziMesh::OWNED);

#if 1
  for (int f=0; f!=f_owned; ++f) {
    printf("face value(%5i)=%20.13e\n",f,(*vec)("face",0,f));
  }
#endif
  MSGF("end of printing ---> "<<prt_str) ;
}

void OverlandFlow::output_flow_rate() {
  MSGF(endl << "begin  output_flow_rate") ;

  // get flux
  Teuchos::RCP<CompositeVector> flux = 
    S_next_->GetFieldData("overland_flux", "overland_flow");

  Teuchos::RCP<CompositeVector> pres = 
    S_next_->GetFieldData("overland_pressure", "overland_flow");

  Teuchos::RCP<CompositeVector> elev = 
    S_next_->GetFieldData("elevation", "overland_flow");

  Teuchos::RCP<CompositeVector> krel = 
    S_next_->GetFieldData("upwind_overland_conductivity", "overland_flow");

  AmanziMesh::Entity_ID_List cells;
  AmanziMesh::Entity_ID_List nodes;
  AmanziGeometry::Point coord0_surf(2), coord1_surf(2);
  AmanziGeometry::Point coord0(3), coord1(3);
  Teuchos::RCP<const AmanziMesh::Mesh_MSTK> surface_mesh =
    Teuchos::rcp_static_cast<const AmanziMesh::Mesh_MSTK>(S_next_->Mesh("surface"));

  //print_vector(S_next_,pres,"output_flow_rate: pressure") ;
  if (S_next_->Mesh("surface")->get_comm()->MyPID() > 0) MPI_Barrier(MPI_COMM_WORLD);


  // compute total flow_rate by integrating the face 
  // with zero gradient conditions
  Functions::BoundaryFunction::Iterator bc;
  double total_flow_rate = 0.;
  int nfaces_owned = flux->size("face",false);
  std::cout << "DATA ON DOWNGRADIENT EDGE " << S_next_->Mesh("surface")->get_comm()->MyPID() << std::endl;
  for (bc=bc_zero_gradient_->begin(); bc!=bc_zero_gradient_->end(); ++bc) {
    int f = bc->first;
    if (f < nfaces_owned) {
      double face_area = S_next_->Mesh("surface")->face_area(f) ;
      AmanziGeometry::Point cent = S_next_->Mesh("surface")->face_centroid(f);
      total_flow_rate += (*flux)("face",0,f);
      printf("f=%5i area=%14.7e flux(%5i)=%14.7e pres(%5i)=%14.7e elev(%5i)=%14.7e krel=%14.7e x=(%14.7e,%14.7e)\n",
             f,face_area,f,(*flux)("face",0,f)/face_area,f,(*pres)("face",0,f),f,(*elev)("face",0,f),(*krel)("face",0,f),cent[0],cent[1]) ;


#if 0
      cells.clear();
      S_next_->Mesh("surface")->face_get_cells(f, AmanziMesh::USED, &cells);
      ASSERT(cells.size() == 1);
      cent = S_next_->Mesh("surface")->cell_centroid(cells[0]);
      printf("  c=%5i x=(%14.7e,%14.7e) elev=%14.7e\n",cells[0],cent[0],cent[1],(*elev)("cell",0,cells[0]));

      nodes.clear();
      S_next_->Mesh("surface")->face_get_nodes(f, &nodes);
      ASSERT(nodes.size() == 2);
      S_next_->Mesh("surface")->node_get_coordinates(nodes[0], &coord0_surf);
      S_next_->Mesh("surface")->node_get_coordinates(nodes[1], &coord1_surf);

      AmanziMesh::Entity_ID node0 = surface_mesh->entity_get_parent(AmanziMesh::NODE, nodes[0]);
      AmanziMesh::Entity_ID node1 = surface_mesh->entity_get_parent(AmanziMesh::NODE, nodes[1]);
      S_next_->Mesh("domain")->node_get_coordinates(node0, &coord0);
      S_next_->Mesh("domain")->node_get_coordinates(node1, &coord1);

      printf("  nodes: surface: (%5i) (%14.7e, %14.7e),                 (%5i) (%14.7e, %14.7e)\n         domain:  (%5i) (%14.7e, %14.7e, %14.7e), (%5i) (%14.7e, %14.7e, %14.7e)\n",
             nodes[0],coord0_surf[0],coord0_surf[1],
             nodes[1],coord1_surf[0],coord1_surf[1],
             node0,coord0[0],coord0[1],coord0[2],
             node1,coord1[0],coord1[1],coord1[2]);
#endif

    }
  }

  fflush(stdout);
  if (S_next_->Mesh("surface")->get_comm()->MyPID() == 0) MPI_Barrier(MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);

  LINE(--);

  // std::cout << "DATA ON OUTER EDGE" << std::endl;
  // for (bc=bc_pressure_->begin(); bc!=bc_pressure_->end(); ++bc) {
  //   int f = bc->first;
  //   if (f < nfaces_owned) {
  //     double face_area = S_next_->Mesh("surface")->face_area(f) ;
  //     printf("f=%5i area=%14.7e flux(%5i)=%14.7e pres(%5i)=%14.7e elev(%5i)=%14.7e krel=%14.7e\n",
  //            f,face_area,f,(*flux)("face",0,f)/face_area,f,(*pres)("face",0,f),f,(*elev)("face",0,f),(*krel)("face",0,f)) ;
  //   }
  // }


#ifdef HAVE_MPI
  double buf = total_flow_rate;
  MPI_Allreduce(&buf, &total_flow_rate, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif

  if (S_next_->Mesh("surface")->get_comm()->MyPID() == 0) {
    PRT(total_flow_rate) ;

    // output file
    string fname("flow_rate.dat");

    // check if the dataset file exists
    bool ok_file(false);
    std::ifstream inpf( fname.c_str() );
    if ( inpf.good() ) { ok_file = bool(true); }
    inpf.close();

    // open the file
    std::ofstream outf;
    if ( !ok_file ) { // if false, open it for output
      outf.open( fname.c_str(), std::ios_base::out );
      outf << "## flow_rate " << endl;
    } else {          // if true,  open it output and append a blank line
      outf.open( fname.c_str(), std::ios_base::app );
    }

    outf << flow_time_/60. << "  " << total_flow_rate << endl;
    outf.close();
  }
}



} // namespace
} // namespace
