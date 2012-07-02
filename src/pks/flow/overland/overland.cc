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
    std::cout << "Region Surface mesh size: " << mesh->get_set_size("ss7",
            AmanziMesh::FACE, AmanziMesh::OWNED) << std::endl;

    // -- Call the MSTK constructor to rip off the surface of the MSTK domain
    // -- mesh.
    std::vector<std::string> setnames(1,"ss7");
    Teuchos::RCP<AmanziMesh::Mesh> surface_mesh =
      Teuchos::rcp(new AmanziMesh::Mesh_MSTK(*mesh,setnames,AmanziMesh::FACE,true,false));

    // push the mesh into state
    S->RegisterMesh("surface", surface_mesh);

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

  // abs perm tensor
  variable_abs_perm_ = false;
  int c_owned = S->Mesh("surface")->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  K_.resize(c_owned);
  for (int c=0; c!=c_owned; ++c) {
    K_[c].init(S->Mesh("surface")->space_dimension(), 1);
    K_[c](0,0) = 1.0;
  }

  // boundary conditions
  Teuchos::ParameterList bc_plist = flow_plist_.sublist("boundary conditions", true);
  FlowBCFactory bc_factory(S->Mesh("surface"), bc_plist);
  bc_pressure_ = bc_factory.CreatePressure();
  bc_zero_gradient_ = bc_factory.CreateZeroGradient();
  bc_flux_ = bc_factory.CreateMassFlux();

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
  S->GetField("overland_velocity", "overland_flow")->set_initialized();
  S->GetField("pres_elev", "overland_flow")->set_initialized();
  S->GetField("rainfall_rate", "overland_flow")->set_initialized();

  // Rel perm is special -- if the mode is symmetric, it needs to be
  // initialized to 1.
  S->GetFieldData("upwind_overland_conductivity","overland_flow")->PutScalar(1.0);
  S->GetField("upwind_overland_conductivity", "overland_flow")->set_initialized();

  // Initialize BCs.
  bc_pressure_->Compute(S->time());
  bc_zero_gradient_->Compute(S->time());
  bc_flux_->Compute(S->time());

  // Initialize operators.
  matrix_->CreateMFDmassMatrices(K_);
  preconditioner_->CreateMFDmassMatrices(K_);

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
    //    dt_ = time_stepper_->time_step(h, solution_);
    time_stepper_->time_step(h, solution_);
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

  if (false) {
    string str_label = string("advance, after time_step----") ;
    PRT(flow_time_) ;
    print_pressure(S_next_,str_label);
    LINE(---) ;
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

  std::cout << "pres cell: " << (*pres)("cell",0) << std::endl;
  std::cout << "pres cell: " << (*pres)("face",1) << ", "
            << (*pres)("face",2) << ", " << (*pres)("face",3) << ", "
            << (*pres)("face",4) << std::endl;
  std::cout << "krel cell: " << (*upwind_conductivity)("cell",0) << std::endl;
  std::cout << "krel cell: " << (*upwind_conductivity)("face",1) << ", "
            << (*upwind_conductivity)("face",2) << ", "
            << (*upwind_conductivity)("face",3) << ", "
            << (*upwind_conductivity)("face",4) << std::endl;


  matrix_->CreateMFDstiffnessMatrices(*upwind_conductivity);
  matrix_->DeriveFlux(*pres_elev, darcy_flux);

#if 0
  print_pressure(S,"commit_state: pressure") ;
  LINE(--) ;

  print_vector2 (S,pressure,elevation,"commit_state: pres & elev") ;
  LINE(--) ;

  print_vector (S,pres_elev,"commit_state: pres_elev") ;
  LINE(--) ;
  print_faceval(S,darcy_flux,"commit_state: darcy_flux") ;
#endif

#if 0
  // cruft to dump overland flow solution
  std::stringstream filename;
  filename << "overland_pressure_" << S->time();
  std::string fname = filename.str();

  Teuchos::RCP<const CompositeVector> pres = S->GetFieldData("overland_pressure");
  Teuchos::RCP<const Epetra_MultiVector> pres_mv = pres->ViewComponent("cell", false);
  EpetraExt::MultiVectorToMatrixMarketFile(fname.c_str(), *pres_mv,"abc","abc",false);
#endif
};

// -- update diagnostics -- used prior to vis
void OverlandFlow::calculate_diagnostics(const Teuchos::RCP<State>& S) {
  // update the cell velocities
  Teuchos::RCP<CompositeVector> velocity = S->GetFieldData("overland_velocity", "overland_flow");
  Teuchos::RCP<const CompositeVector> flux = S->GetFieldData("overland_flux");
  matrix_->DeriveCellVelocity(*flux, velocity);
};


// -----------------------------------------------------------------------------
// Use the physical rel perm (on cells) to update a work vector for rel perm.
//
//   This deals with upwinding, etc.
// -----------------------------------------------------------------------------
void OverlandFlow::UpdatePermeabilityData_(const Teuchos::RCP<State>& S) {
  // get reference to relative permeability
  Teuchos::RCP<const CompositeVector> rel_perm =
    S->GetFieldData("overland_conductivity");

  Teuchos::RCP<CompositeVector> upwind_conductivity =
    S->GetFieldData("upwind_overland_conductivity", "overland_flow");

  rel_perm->ScatterMasterToGhosted();

  Teuchos::RCP<const CompositeVector> flux = S->GetFieldData("overland_flux");
  CalculateRelativePermeabilityUpwindFlux_(S, *flux, *rel_perm,
          upwind_conductivity);
};


// -----------------------------------------------------------------------------
// Defines upwinded relative permeabilities for faces using a given flux.
// -----------------------------------------------------------------------------
void OverlandFlow::CalculateRelativePermeabilityUpwindFlux_(
        const Teuchos::RCP<State>& S,
        const CompositeVector& flux,
        const CompositeVector& conductivity,
        const Teuchos::RCP<CompositeVector>& upwind_conductivity ) {
  AmanziMesh::Entity_ID_List faces;
  std::vector<int> dirs;

  upwind_conductivity->ViewComponent("face",true)->PutScalar(0.0);
  int c_owned = upwind_conductivity->size("cell");
  for (int c=0; c!=c_owned; ++c) {
    S->Mesh("surface")->cell_get_faces_and_dirs(c, &faces, &dirs);

    for (int n=0; n!=faces.size(); ++n) {
      int f = faces[n];
      if ((flux("face",f) * dirs[n] >= 0.0) ||
          (bc_markers_[f] != Operators::MFD_BC_NULL)) {
        (*upwind_conductivity)("face",f) = conductivity("cell",c);
      }
    }
  }

  //print_vector(S, upwind_conductivity, "face conductivity");
}


// -----------------------------------------------------------------------------
// Interpolate pressure ICs on cells to ICs for lambda (faces).
// -----------------------------------------------------------------------------
void OverlandFlow::DeriveFaceValuesFromCellValues_(const Teuchos::RCP<State>& S,
                                                   const Teuchos::RCP<CompositeVector> & pres) {
  AmanziMesh::Entity_ID_List cells;

  int f_owned = pres->size("face");
  for (int f=0; f!=f_owned; ++f) {
    cells.clear();
    S->Mesh("surface")->face_get_cells(f, AmanziMesh::USED, &cells);
    int ncells = cells.size();
    // ---
    double face_value = 0.0;
    for (int n=0; n!=ncells; ++n) {
      face_value += (*pres)("cell",0,cells[n]);
    }
    (*pres)("face",0,f) = face_value / ncells;
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
  int nfaces = pres->size("face");
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

#if 1
  for (int c=0; c!=c_owned; ++c) {
    printf("cell pressure(%5i)=%14.7e\n",c,(*pressure)("cell",0,c));
  }
  LINE(--) ;
  for (int f=0; f!=f_owned; ++f) {
    printf("face pressure(%5i)=%14.7e\n",f,(*pressure)("face",0,f));
  }
#endif

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


  //print_vector(S_next_,pres,"output_flow_rate: pressure") ;

  // compute total flow_rate by integrating the face 
  // with zero gradient conditions
  Functions::BoundaryFunction::Iterator bc;
  double total_flow_rate = 0.;
  std::cout << "DATA ON DOWNGRADIENT EDGE" << std::endl;
  for (bc=bc_zero_gradient_->begin(); bc!=bc_zero_gradient_->end(); ++bc) {
    int f = bc->first;
    //PRT(f) ;
    double face_area = S_next_->Mesh("surface")->face_area(f) ;
    total_flow_rate += (*flux)("face",0,f);
    printf("f=%5i area=%14.7e flux(%5i)=%14.7e pres(%5i)=%14.7e elev(%5i)=%14.7e krel=%14.7e\n",
           f,face_area,f,(*flux)("face",0,f)/face_area,f,(*pres)("face",0,f),f,(*elev)("face",0,f),(*krel)("face",0,f)) ;
  }

  LINE(--);

  std::cout << "DATA ON OUTER EDGE" << std::endl;
  for (bc=bc_pressure_->begin(); bc!=bc_pressure_->end(); ++bc) {
    int f = bc->first;
    //PRT(f) ;
    double face_area = S_next_->Mesh("surface")->face_area(f) ;
    //    total_flow_rate += (*flux)("face",0,f);
    printf("f=%5i area=%14.7e flux(%5i)=%14.7e pres(%5i)=%14.7e elev(%5i)=%14.7e krel=%14.7e\n",
           f,face_area,f,(*flux)("face",0,f)/face_area,f,(*pres)("face",0,f),f,(*elev)("face",0,f),(*krel)("face",0,f)) ;
  }


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



} // namespace
} // namespace
