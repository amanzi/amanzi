/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/*
This is the flow component of the Amanzi code. 
License: BSD
Authors: Neil Carlson (version 1) 
         Konstantin Lipnikov (version 2) (lipnikov@lanl.gov)
         Ethan Coon (ATS version) (ecoon@lanl.gov)
*/

#include "EpetraExt_MultiVectorOut.h"
#include "Epetra_MultiVector.h"

#include "bdf1_time_integrator.hh"
#include "flow_bc_factory.hh"

#include "overland.hh"

namespace Amanzi {
namespace Flow {

#if 0 // a trick for xemacs indent algorithm
}}
#endif

RegisteredPKFactory<OverlandFlow> OverlandFlow::reg_("overland flow");

// constructor
OverlandFlow::OverlandFlow(Teuchos::ParameterList& flow_plist,
                           const Teuchos::RCP<State>& S,
                           const Teuchos::RCP<TreeVector>& solution) :
    flow_plist_(flow_plist) {

  // just the extras...
  // data layouts for fields
  std::vector<AmanziMesh::Entity_kind> locations2(2);
  std::vector<std::string> names2(2);
  locations2[0] = AmanziMesh::CELL;
  locations2[1] = AmanziMesh::FACE;
  names2[0] = "cell";
  names2[1] = "face";

  // -- primary variable: pressure on both cells and faces, ghosted, with 1 dof
  std::vector< std::vector<std::string> > subfield_names(2);
  S->RequireField("overland_pressure", "overland_flow", names2, locations2, 1, true);
  Teuchos::RCP<CompositeVector> pressure = S->GetFieldData("overland_pressure", "overland_flow");

  // update the tree-vector solution with flow's primary variable
  solution->set_data(pressure);
  solution_ = solution;

  // --------------------------------------------------------------
  // SET THE TWO MAIN SECONDARY VARS
  // --------------------------------------------------------------
  // -- secondary variable: elevation on both cells and faces, ghosted, with 1 dof
  S->RequireField("elevation", "overland_flow", names2, locations2, 1, true);
  //Teuchos::RCP<CompositeVector> elevation = S->GetFieldData("elevation", "overland_flow");

  // -- secondary variable: pres_elev on both cells and faces, ghosted, with 1 dof
  S->RequireField("pres_elev", "overland_flow", names2, locations2, 1, true);
  //Teuchos::RCP<CompositeVector> pres_elev = S->GetFieldData("pres_elev", "overland_flow");

  // --------------------------------------------------------------

  // -- other secondary variables
  S->RequireField("overland_flux", "overland_flow", AmanziMesh::FACE, 1, true);
  S->RequireField("overland_velocity", "overland_flow", AmanziMesh::CELL, 3, false);
  S->RequireField("overland_conductivity", "overland_flow", AmanziMesh::CELL, 1, true);

  // -- independent variables not owned by this PK
  S->RequireField("cell_volume", AmanziMesh::CELL, 1, true);

  // -- work vectors
  //  S->RequireField("upwind_overland_conductivity", "overland_flow", AmanziMesh::FACE, 1, true);
  S->RequireField("upwind_overland_conductivity", "overland_flow", names2, locations2, 1, true);
  S->GetRecord("upwind_overland_conductivity","overland_flow")->set_io_vis(false);

  // abs perm tensor
  variable_abs_perm_ = false;
  int c_owned = S->mesh()->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  K_.resize(c_owned);
  for (int c=0; c!=c_owned; ++c) {
    K_[c].init( S->mesh()->space_dimension(), 1 );
  }

  // boundary conditions
  Teuchos::ParameterList bc_plist = flow_plist_.sublist("boundary conditions", true);
  FlowBCFactory bc_factory(S->mesh(), bc_plist);
  bc_pressure_ = bc_factory.CreatePressure();
  bc_zero_gradient_ = bc_factory.CreateZeroGradient();
  bc_flux_ = bc_factory.CreateMassFlux();

  // operator for the diffusion terms
  Teuchos::ParameterList mfd_plist = flow_plist_.sublist("Diffusion");
  matrix_ = Teuchos::rcp(new Operators::MatrixMFD(mfd_plist, S->mesh()));

  bool symmetric = false;
  matrix_->SetSymmetryProperty(symmetric);
  matrix_->SymbolicAssembleGlobalMatrices();

  // preconditioner for the NKA system
  Teuchos::ParameterList mfd_pc_plist = flow_plist_.sublist("Diffusion PC");
  preconditioner_ = Teuchos::rcp(new Operators::MatrixMFD(mfd_pc_plist, S->mesh()));
  preconditioner_->SetSymmetryProperty(symmetric);
  preconditioner_->SymbolicAssembleGlobalMatrices();
  Teuchos::ParameterList mfd_pc_ml_plist = mfd_pc_plist.sublist("ML Parameters");
  preconditioner_->InitMLPreconditioner(mfd_pc_ml_plist);

#if 0
  { // -------->>>>> remove
    DBGF(CTOR) ; 
    Teuchos::RCP<CompositeVector> my_pres =  S->GetFieldData("overland_pressure", "overland_flow");

    int c_owned = S_next_->mesh()->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
    int f_owned = S_next_->mesh()->num_entities(AmanziMesh::FACE, AmanziMesh::OWNED);
    VAL(c_owned) ; PRT(f_owned) ;
    for (int c=0; c!=c_owned; ++c) {
      printf("my_pres(%5i)=%14.7e\n",c,(*my_pres)("cell",0,c));
    }
    LINE(--) ;
    for (int f=0; f!=f_owned; ++f) {
      printf("my_pres(%5i)=%14.7e\n",f,(*my_pres)("face",0,f));
    }
  }
  exit(0) ;
#endif



};

// -- Initialize the PK
void OverlandFlow::initialize(const Teuchos::RCP<State>& S) {
  // initial timestep size
  dt_ = flow_plist_.get<double>("Initial time step", 1.);

  // initialize boundary conditions
  int nfaces = S->mesh()->num_entities(AmanziMesh::FACE, AmanziMesh::USED);
  bc_markers_.resize(nfaces, Operators::MFD_BC_NULL);
  bc_values_.resize(nfaces, 0.0);

  // update face pressures as a hint?
  Teuchos::RCP<CompositeVector> pres =
    S->GetFieldData("overland_pressure", "overland_flow");
  DeriveFaceValuesFromCellValues_(S, pres ) ;

  //print_pressure( S, "INITIALIZE" ) ;

  // initialize the elevation field
  SetUpElevation_(S);
  S->GetRecord("elevation","overland_flow")->set_initialized();

  // declare secondary variables initialized, as they will be done by
  // the commit_state call
  S->GetRecord("overland_conductivity","overland_flow")->set_initialized();
  S->GetRecord("overland_flux",           "overland_flow")->set_initialized();
  S->GetRecord("overland_velocity",       "overland_flow")->set_initialized();
  S->GetRecord("upwind_overland_conductivity",  "overland_flow")->set_initialized();
  S->GetRecord("pres_elev","overland_flow")->set_initialized();

  // set the abs perm tensor to 1 -- it is not needed
  for (int c=0; c!=K_.size(); ++c) {
    K_[c](0,0) = 1.0;
  }

  // initialize operators
  matrix_->CreateMFDmassMatrices(K_);
  preconditioner_->CreateMFDmassMatrices(K_);

  // initialize the timestepper
  state_to_solution(S, solution_);
  atol_ = flow_plist_.get<double>("Absolute error tolerance",1e0);
  rtol_ = flow_plist_.get<double>("Relative error tolerance",1e0);

  if (!flow_plist_.get<bool>("Strongly Coupled PK", false)) {
    // -- instantiate time stepper
    Teuchos::RCP<Teuchos::ParameterList> bdf1_plist_p =
      Teuchos::rcp(new Teuchos::ParameterList(flow_plist_.sublist("Time integrator")));
    time_stepper_ = Teuchos::rcp(new BDF1TimeIntegrator(this, bdf1_plist_p, solution_));
    time_step_reduction_factor_ = bdf1_plist_p->get<double>("time step reduction factor");

    // -- initialize time derivative
    Teuchos::RCP<TreeVector> solution_dot = Teuchos::rcp( new TreeVector(*solution_));
    solution_dot->PutScalar(0.0);

    // -- set initial state
    time_stepper_->set_initial_state(S->time(), solution_, solution_dot);
  }
};

// Pointer copy of state to solution
void OverlandFlow::state_to_solution(const Teuchos::RCP<State>& S,
                                     const Teuchos::RCP<TreeVector>& solution) {
  solution->set_data(S->GetFieldData("overland_pressure", "overland_flow"));
};

// Pointer copy concentration fields from solution vector into state.  Used within
// compute_f() of strong couplers to set the current iterate in the state for
// use with other PKs.
void OverlandFlow::solution_to_state(const Teuchos::RCP<TreeVector>& solution,
                                     const Teuchos::RCP<State>& S) {
  S->SetData("overland_pressure", "overland_flow", solution->data());
};

// -- Advance from state S0 to state S1 at time S0.time + dt.
bool OverlandFlow::advance(double dt) {
  {
    //string str_label = string("advance, before time_step----") ;
    //print_pressure(S_next_,str_label);
  }
    
  state_to_solution(S_next_, solution_);

  // take a bdf timestep
  double h = dt;

  try {
    //    dt_ = time_stepper_->time_step(h, solution_);
    time_stepper_->time_step(h, solution_);
  } catch (Exceptions::Amanzi_exception &error) {
    std::cout << "Timestepper called error: " << error.what() << std::endl;
    if (error.what() == std::string("BDF time step failed")) {
      // try cutting the timestep
      dt_ = dt_*time_step_reduction_factor_;
      return true;
    } else {
      throw error;
    }
  }

  {
    string str_label = string("advance, after time_step----") ;
    print_pressure(S_next_,str_label);
    //exit(0) ;
  }

  time_stepper_->commit_solution(h, solution_);
  solution_to_state(solution_, S_next_);
  commit_state(h,S_next_);
  return false;
};

// -- commit the state... 
void OverlandFlow::commit_state(double dt, const Teuchos::RCP<State>& S) {
  // update secondary variables for IC
  UpdateSecondaryVariables_(S);

  // update the flux
  UpdatePermeabilityData_(S);
  Teuchos::RCP<const CompositeVector> upwind_conductivity =
    S->GetFieldData("upwind_overland_conductivity");

  Teuchos::RCP<const CompositeVector> pres_elev = S->GetFieldData("pres_elev");
  Teuchos::RCP<CompositeVector> darcy_flux =
    S->GetFieldData("overland_flux", "overland_flow");

  matrix_->CreateMFDstiffnessMatrices(*upwind_conductivity);
  matrix_->DeriveFlux(*pres_elev, darcy_flux);

  // cruft to dump overland flow solution
  std::stringstream filename;
  filename << "overland_pressure_" << S->time();
  std::string fname = filename.str();

  Teuchos::RCP<const CompositeVector> pres = S->GetFieldData("overland_pressure");
  Teuchos::RCP<const Epetra_MultiVector> pres_mv = pres->ViewComponent("cell", false);
  EpetraExt::MultiVectorToMatrixMarketFile(fname.c_str(), *pres_mv,"abc","abc",false);
};

// -- update diagnostics -- used prior to vis
void OverlandFlow::calculate_diagnostics(const Teuchos::RCP<State>& S) {
  // update the cell velocities
  Teuchos::RCP<CompositeVector> velocity = S->GetFieldData("overland_velocity", "overland_flow");
  Teuchos::RCP<const CompositeVector> flux = S->GetFieldData("overland_flux");
  matrix_->DeriveCellVelocity(*flux, velocity);
};

// relative permeability methods
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

/* ******************************************************************
* Defines upwinded relative permeabilities for faces using a given flux.
****************************************************************** */
void OverlandFlow::CalculateRelativePermeabilityUpwindFlux_(
        const Teuchos::RCP<State>& S,
        const CompositeVector& flux,
        const CompositeVector& conductivity,
        const Teuchos::RCP<CompositeVector>& upwind_conductivity ) {
  AmanziMesh::Entity_ID_List faces;
  std::vector<int> dirs;

  upwind_conductivity->ViewComponent("cell",true)->PutScalar(1.0);

  int c_owned = S->mesh()->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  for (int c=0; c!=c_owned; ++c) {
    S->mesh()->cell_get_faces_and_dirs(c, &faces, &dirs);

    for (int n=0; n!=faces.size(); ++n) {
      int f = faces[n];
      if ((flux("face",0,f) * dirs[n] >= 0.0) ||
          (bc_markers_[f] != Operators::MFD_BC_NULL)) {
        (*upwind_conductivity)("face",0,f) = conductivity("cell",0,c);
      }
    }
  }

  //print_vector(S, upwind_conductivity, "face conductivity");
}

void OverlandFlow::DeriveFaceValuesFromCellValues_(const Teuchos::RCP<State>& S,
                                                   const Teuchos::RCP<CompositeVector> & pres) {
  AmanziMesh::Entity_ID_List cells;

  int f_owned = S->mesh()->num_entities(AmanziMesh::FACE, AmanziMesh::OWNED);
  for (int f=0; f!=f_owned; ++f) {
    cells.clear();
    S->mesh()->face_get_cells(f, AmanziMesh::USED, &cells);
    int ncells = cells.size();
    // ---
    double face_value = 0.0;
    for (int n=0; n!=ncells; ++n) {
      face_value += (*pres)("cell",0,cells[n]);
    }
    (*pres)("face",0,f) = face_value / ncells;
  }
};

/* ******************************************************************
 * Add a boundary marker to used faces.
 ****************************************************************** */
void OverlandFlow::UpdateBoundaryConditions_( const Teuchos::RCP<State>& S,
                                              const CompositeVector & pres ) {

  const CompositeVector & elevation = *(S->GetFieldData("elevation"));

  for (int n=0; n!=bc_markers_.size(); ++n) {
    bc_markers_[n] = Operators::MFD_BC_NULL;
    bc_values_[n] = 0.0;
  }

  BoundaryFunction::Iterator bc;
  for (bc=bc_pressure_->begin(); bc!=bc_pressure_->end(); ++bc) {
    int f = bc->first;
    bc_markers_[f] = Operators::MFD_BC_DIRICHLET;
    bc_values_[f] = bc->second + elevation("face",0,f);
  }

  AmanziMesh::Entity_ID_List cells;
  for (bc=bc_zero_gradient_->begin(); bc!=bc_zero_gradient_->end(); ++bc) {
    int f = bc->first;

    cells.clear();
    S->mesh()->face_get_cells(f, AmanziMesh::USED, &cells);
    int ncells = cells.size();
    ASSERT( ncells==1 ) ;

    bc_markers_[f] = Operators::MFD_BC_DIRICHLET;
    bc_values_[f] = pres("cell",0,cells[0]) + elevation("face",0,f);
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
                                            CompositeVector & pres) {
  int nfaces = S->mesh()->num_entities(AmanziMesh::FACE, AmanziMesh::OWNED);
  for (int f=0; f!=nfaces; ++f) {
    if (bc_markers_[f] == Operators::MFD_BC_DIRICHLET) {
      pres("face",0,f) = bc_values_[f];
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
  manning[0] = 0.25 ;
}
double OverlandFlow::TestOneElevation_( double x, double y ) {
  double Lx = 600./3.28084 ;
  return slope_x[0]*( Lx - x ) ;
}

void OverlandFlow::TestOneSetElevation_( const Teuchos::RCP<State>& S ) {
  CompositeVector& elev = *(S ->GetFieldData("elevation","overland_flow"));

  // set the elevation vector's cells
  int c_owned = S->mesh()->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  for (int c=0; c!=c_owned; ++c) {
    Amanzi::AmanziGeometry::Point pc = S->mesh()->cell_centroid(c);
    elev("cell",0,c) = TestOneElevation_(pc[0], pc[1]);
  }

  // set the elevation vector's faces
  int f_owned = S->mesh()->num_entities(AmanziMesh::FACE, AmanziMesh::OWNED);
  for (int f=0; f!=f_owned; ++f) {
    Amanzi::AmanziGeometry::Point pf = S->mesh()->face_centroid(f);
    elev("face",0,f) = TestOneElevation_(pf[0], pf[1]);
  }
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
  slope_y[2] = 0.02  ;
  // ---
  manning.resize(3) ;
  manning[0] = 0.015 ;
  manning[1] = 0.15  ;
  manning[2] = 0.015 ; 
}

double OverlandFlow::TestTwoElevation_( double x, double y ) { 
  int izn = TestTwoZoneFlag_( x, y ) ;
  double retval = 0. ;
  switch ( izn ) {
  case 0: retval = -slope_x[0]*(x-zone_x[1]) + slope_y[0]*y ; break ;
  case 1: retval =                           + slope_y[1]*y ; break ;
  case 2: retval = +slope_x[2]*(x-zone_x[2]) + slope_y[2]*y ; break ;
  default: assert(false) ;
  }
  return retval ;
}

void OverlandFlow::TestTwoSetElevation_( const Teuchos::RCP<State>& S ) {

  CompositeVector & elev = *(S ->GetFieldData("elevation","overland_flow"));
  
  // add elevation to pres_elev, cell dofs
  int c_owned = S->mesh()->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  for (int c=0; c!=c_owned; ++c) {
    Amanzi::AmanziGeometry::Point pc = S->mesh()->cell_centroid( c ) ;
    elev("cell",0,c) = TestTwoElevation_( pc[0], pc[1] ) ;
  }  

  // add elevation to pres_elev, face dofs
  int f_owned = S->mesh()->num_entities(AmanziMesh::FACE, AmanziMesh::OWNED);
  for (int f=0; f!=f_owned; ++f) {
    Amanzi::AmanziGeometry::Point pf = S->mesh()->face_centroid( f ) ;    
    elev("face",0,f) = TestTwoElevation_( pf[0], pf[1] ) ;
  }  
}

void OverlandFlow::SetUpElevation_( const Teuchos::RCP<State>& S ) {

#if 1

  // TEST 1
  load_value = 2. * 0.0254 / 3600. ;
  //load_value = 1.;
  t_rain     = 1800.;
  TestOneSetUpElevationPars_();
  TestOneSetElevation_( S );

#else
  
  // TEST 2
  load_value = 3.e-6; // [m/s] == 10.8e-3 / 3600.
  t_rain     = 5400.; // [s], duration of raining event

  TestTwoSetUpElevationPars_();
  TestTwoSetElevation_( S );

#endif

}  

// loading term for the raining events
double OverlandFlow::rhs_load_value( double t ) {
  return t<t_rain ? load_value : 0. ; 
}


// OTHER DEBUGGING STUFF

void OverlandFlow::print_pressure( const Teuchos::RCP<State>& S, string prt_str ) const {

  MSGF("begin  ---> "<<prt_str) ;

  int c_owned = S->mesh()->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  int f_owned = S->mesh()->num_entities(AmanziMesh::FACE, AmanziMesh::OWNED);

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


void OverlandFlow::print_vector( const Teuchos::RCP<State>& S, const Teuchos::RCP<const CompositeVector>& pressure, string prt_str ) const {

  MSGF("begin  ---> "<<prt_str) ;

  int c_owned = S->mesh()->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  int f_owned = S->mesh()->num_entities(AmanziMesh::FACE, AmanziMesh::OWNED);

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






} // namespace
} // namespace
