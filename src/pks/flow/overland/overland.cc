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
  subfield_names[0].resize(1);
  subfield_names[0][0] = "cell";
  subfield_names[1].resize(1);
  subfield_names[1][0] = "face";
  S->RequireField("overland_pressure", "overland_flow", names2, locations2, 1, true);
  Teuchos::RCP<CompositeVector> pressure = S->GetFieldData("overland_pressure", "overland_flow");
  pressure->set_subfield_names(subfield_names);

  // update the tree-vector solution with flow's primary variable
  solution->set_data(pressure);
  solution_ = solution;
  
  // --------------------------------------------------------------
  // SET THE TWO MAIN SECONDARY VARS
  // --------------------------------------------------------------
  // -- secondary variable: elevation on both cells and faces, ghosted, with 1 dof
  //std::vector< std::vector<std::string> > subfield_names(2);
  subfield_names[0].resize(1);
  subfield_names[0][0] = "elevation_cell";
  subfield_names[1].resize(1);
  subfield_names[1][0] = "elevation_face";
  S->RequireField("elevation", "overland_flow", names2, locations2, 1, true);
  Teuchos::RCP<CompositeVector> elevation = S->GetFieldData("elevation", "overland_flow");
  pressure->set_subfield_names(subfield_names);
  
  // -- secondary variable: pres_elev on both cells and faces, ghosted, with 1 dof
  //std::vector< std::vector<std::string> > subfield_names(2);
  subfield_names[0].resize(1);
  subfield_names[0][0] = "pres_elev_cell";
  subfield_names[1].resize(1);
  subfield_names[1][0] = "pres_elev_face";
  S->RequireField("pres_elev", "overland_flow", names2, locations2, 1, true);
  Teuchos::RCP<CompositeVector> pres_elev = S->GetFieldData("pres_elev", "overland_flow");
  pressure->set_subfield_names(subfield_names);
  // --------------------------------------------------------------
  
  // -- other secondary variables
  S->RequireField("darcy_flux",            "overland_flow", AmanziMesh::FACE, 1, true);
  S->RequireField("darcy_velocity",        "overland_flow", AmanziMesh::CELL, 3, false);
  S->RequireField("relative_permeability", "overland_flow", AmanziMesh::CELL, 1, true);
  
  // scalar factors
  //S->RequireScalar("atmospheric_pressure", "overland_flow");
  //S->RequireScalar("porosity",             "overland_flow");
  //S->RequireScalar("molar_density_liquid", "overland_flow");
  //S->RequireScalar("density_liquid",       "overland_flow");
  //S->RequireScalar("viscosity_liquid",     "overland_flow");

  // -- independent variables not owned by this PK
  S->RequireField("cell_volume",          AmanziMesh::CELL, 1, true);
  //S->RequireGravity();
  
  // -- work vectors
  S->RequireField("rel_perm_faces", "overland_flow", AmanziMesh::FACE, 1, true);
  S->GetRecord("rel_perm_faces","overland_flow")->set_io_vis(false);

  // abs perm tensor
  variable_abs_perm_ = false;
  int c_owned = S->mesh()->count_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  K_.resize(c_owned);
  for (int c=0; c!=c_owned; ++c) { 
    K_[c].init( S->mesh()->space_dimension(), 1 );
  }
  
  // // -- water retention models
  // Teuchos::ParameterList wrm_plist = flow_plist_.sublist("Water Retention Models");
  // // count the number of region-model pairs
  // int wrm_count = 0;
  // for (Teuchos::ParameterList::ConstIterator i=wrm_plist.begin(); 
  //      i!=wrm_plist.end(); ++i) {
  //   if ( wrm_plist.isSublist(wrm_plist.name(i)) ) {
  //     wrm_count++;
  //   } else {
  //     std::string message("OverlandFlow: water retention model contains an entry that is not a sublist.");
  //     Exceptions::amanzi_throw(message);
  //   }
  // }
  // wrm_.resize(wrm_count);

  // // instantiate the region-model pairs
  // FlowRelations::WRMFactory wrm_factory;
  // int iblock = 0;
  // for (Teuchos::ParameterList::ConstIterator i=wrm_plist.begin(); i!=wrm_plist.end(); ++i) {
  //   Teuchos::ParameterList wrm_region_list = wrm_plist.sublist(wrm_plist.name(i));
  //   std::string region = wrm_region_list.get<std::string>("Region");
  //   wrm_[iblock] = Teuchos::rcp(new WRMRegionPair(region, wrm_factory.createWRM(wrm_region_list)));
  //   iblock++;
  // }

  // boundary conditions
  Teuchos::ParameterList bc_plist = flow_plist_.sublist("boundary conditions", true);
  FlowBCFactory bc_factory(S->mesh(), bc_plist);
  bc_pressure_ = bc_factory.CreatePressure();
  bc_zero_gradient_ = bc_factory.CreateZeroGradient();
  bc_flux_ = bc_factory.CreateMassFlux();
  // head boundary conditions not yet implemented, as they depend on constant density
  // bc_head_ = bc_factory.CreateStaticHead(0.0, 0.0, g);

  // // Relative permeability method
  // string method_name = flow_plist_.get<string>("Relative permeability method", "Upwind with gravity");
  // bool symmetric = false;
  // if (method_name == "Upwind with gravity") {
  //   Krel_method_ = FLOW_RELATIVE_PERM_UPWIND_GRAVITY;
  // } else if (method_name == "Cell centered") {
  //   Krel_method_ = FLOW_RELATIVE_PERM_CENTERED;
  //   symmetric = true;
  // } else if (method_name == "Upwind with Darcy flux") {
  //   Krel_method_ = FLOW_RELATIVE_PERM_UPWIND_DARCY_FLUX;
  // } else if (method_name == "Arithmetic mean") {
  //   Krel_method_ = FLOW_RELATIVE_PERM_ARITHMETIC_MEAN;
  // }

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
};

// -- Initialize owned (dependent) variables.
void OverlandFlow::initialize(const Teuchos::RCP<State>& S) {
  // initial timestep size
  dt_ = flow_plist_.get<double>("Initial time step", 1.);

  // initialize boundary conditions
  int nfaces = S->mesh()->num_entities(AmanziMesh::FACE, AmanziMesh::USED);
  bc_markers_.resize(nfaces, Operators::MFD_BC_NULL);
  bc_values_.resize(nfaces, 0.0);

  // update face pressures as a hint?
  Teuchos::RCP<CompositeVector> pres = S->GetFieldData("overland_pressure", "overland_flow") ;
  DeriveFaceValuesFromCellValues_(S, pres ) ;

  // declare secondary variables initialized, as they will be done by
  // the commit_state call
  //S->GetRecord("saturation_liquid",    "overland_flow")->set_initialized();
  S->GetRecord("relative_permeability","overland_flow")->set_initialized();
  S->GetRecord("darcy_flux",           "overland_flow")->set_initialized();
  S->GetRecord("darcy_velocity",       "overland_flow")->set_initialized();

  // rel perm is special -- if the mode is symmetric, it needs to be
  // initialized to 1
  S->GetFieldData("rel_perm_faces","overland_flow")->PutScalar(1.0);
  S->GetRecord   ("rel_perm_faces","overland_flow")->set_initialized();

  // initialize the timestepper
  state_to_solution(S, solution_);
  atol_ = flow_plist_.get<double>("Absolute error tolerance",1e-5);
  rtol_ = flow_plist_.get<double>("Relative error tolerance",1e-5);

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

  // -- set the elevation field
  S->GetRecord("pres_elev","overland_flow")->set_initialized();
  S->GetRecord("elevation","overland_flow")->set_initialized();
  SetUpElevation_( S ) ;
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
  state_to_solution(S_next_, solution_);

  // take a bdf timestep
  double h = dt;
  
  try {
    dt_ = time_stepper_->time_step(h, solution_);
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

  time_stepper_->commit_solution(h, solution_);
  solution_to_state(solution_, S_next_);
  commit_state(h,S_next_);
  return false;
};

// -- commit the state... not sure how/when this gets used
void OverlandFlow::commit_state(double dt, const Teuchos::RCP<State>& S) {
  // update secondary variables for IC
  UpdateSecondaryVariables_(S);

  // update the flux
  UpdatePermeabilityData_(S);
  Teuchos::RCP<const CompositeVector> rel_perm_faces = S->GetFieldData("rel_perm_faces");

  Teuchos::RCP<const CompositeVector> pres_elev      = S->GetFieldData("pres_elev");
  Teuchos::RCP<      CompositeVector> darcy_flux     = S->GetFieldData("darcy_flux", "overland_flow");
  
  matrix_->CreateMFDstiffnessMatrices(K_, rel_perm_faces);
  matrix_->DeriveFlux(*pres_elev, darcy_flux);
  //matrix_->DeriveFlux(*pres, darcy_flux);
  //AddGravityFluxesToVector_(S, darcy_flux);


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
  Teuchos::RCP<CompositeVector> velocity = S->GetFieldData("darcy_velocity", "overland_flow");
  Teuchos::RCP<const CompositeVector> flux = S->GetFieldData("darcy_flux");
  matrix_->DeriveCellVelocity(*flux, velocity);
};

// relative permeability methods
void OverlandFlow::UpdatePermeabilityData_(const Teuchos::RCP<State>& S) {
  
  // set absolute permeability
  SetAbsolutePermeabilityTensor_(S);

  // get reference to relative permeability
  Teuchos::RCP<const CompositeVector> rcp_rel_perm       = S->GetFieldData("relative_permeability");
  Teuchos::RCP<      CompositeVector> rcp_rel_perm_faces = S->GetFieldData("rel_perm_faces", "overland_flow");
  
  const CompositeVector & rel_perm       = *rcp_rel_perm ;
  CompositeVector       & rel_perm_faces = *rcp_rel_perm_faces ;
  
  rel_perm . ScatterMasterToGhosted();
  
  Teuchos::RCP<const CompositeVector> rcp_flux = S->GetFieldData("darcy_flux");
  const CompositeVector & flux = *rcp_flux ;
  Calculate_Relative_Permeability_Upwind_Flux_(S, flux, rel_perm, rel_perm_faces);
};

// // relative permeability methods
// void OverlandFlow::UpdatePermeabilityData_(const Teuchos::RCP<State>& S) {
  
//   // set absolute permeability
//   SetAbsolutePermeabilityTensor_(S);

//   // get reference to relative permeability
//   Teuchos::RCP<const CompositeVector> rcp_rel_perm       = S->GetFieldData("relative_permeability");
//   Teuchos::RCP<      CompositeVector> rcp_rel_perm_faces = S->GetFieldData("rel_perm_faces", "overland_flow");
  
//   const CompositeVector & rel_perm       = *rcp_rel_perm ;
//   CompositeVector       & rel_perm_faces = *rcp_rel_perm_faces ;
  
//   rel_perm . ScatterMasterToGhosted();
  
//   double rho  = *(S->GetScalarData("molar_density_liquid"));
//   double visc = *(S->GetScalarData("viscosity_liquid"));
  
//   if (Krel_method_ == FLOW_RELATIVE_PERM_CENTERED) {
    
//     // symmetric method, no faces needed
//     double scaling = rho / visc;
//     for (int c=0; c!=K_.size(); ++c) {
//       K_[c] *= rel_perm(c) * scaling ;
//     }
    
//   } else {
    
//     // faces needed
//     if (Krel_method_ == FLOW_RELATIVE_PERM_UPWIND_GRAVITY) {

//       Calculate_Relative_Permeability_Upwind_Gravity_(S, rel_perm, rel_perm_faces);
      
//     } else if (Krel_method_ == FLOW_RELATIVE_PERM_UPWIND_DARCY_FLUX) {
      
//       Teuchos::RCP<const CompositeVector> rcp_flux = S_->GetFieldData("darcy_flux");
//       const CompositeVector & flux = *rcp_flux ;
//       Calculate_Relative_Permeability_Upwind_Flux_(S, flux, rel_perm, rel_perm_faces);
      
//     } else if (Krel_method_ == FLOW_RELATIVE_PERM_ARITHMETIC_MEAN) {
      
//       Calculate_Relative_Permeability_Arithmetic_Mean_(S, rel_perm, rel_perm_faces);

//     }
    
//     // update K with just rho/mu
//     for (int c=0; c!=K_.size(); ++c) {
//       K_[c] *= rho / visc;
//     }
    
//   }
// };


// /* ******************************************************************
// * Defines upwinded relative permeabilities for faces using gravity.
// ****************************************************************** */
// void OverlandFlow::Calculate_Relative_Permeability_Upwind_Gravity_( const Teuchos::RCP<State> & S,
//                                                                     const CompositeVector & rel_perm_cells,
//                                                                     CompositeVector & rel_perm_faces ) {
//   AmanziMesh::Entity_ID_List faces;
//   std::vector<int> dirs;
  
//   Teuchos::RCP<const Epetra_Vector> g_vec = S->GetConstantVectorData("gravity");
//   AmanziGeometry::Point gravity(g_vec->MyLength());
//   for (int i=0; i!=g_vec->MyLength(); ++i) { 
//     gravity[i] = (*g_vec)[i];
//   }
  
//   rel_perm_faces.PutScalar(0.0);
//   int c_owned = S->mesh()->count_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
//   for (int c=0; c!=c_owned; ++c) {
//     S->mesh()->cell_get_faces_and_dirs(c, &faces, &dirs);

//     AmanziGeometry::Point Kgravity = K_[c] * gravity;

//     for (int n=0; n!=faces.size(); ++n) {
//       int f = faces[n];
//       const AmanziGeometry::Point& normal = S->mesh()->face_normal(f);
//       if ((normal * Kgravity) * dirs[n] >= 0.0) {
//         rel_perm_faces(f) = rel_perm_cells(c);
//       } else if (bc_markers_[f] != Operators::MFD_BC_NULL) {
//         rel_perm_faces(f) = rel_perm_cells(c);
//       }
//     }
//   }
// }

/* ******************************************************************
* Defines upwinded relative permeabilities for faces using a given flux.
****************************************************************** */
void OverlandFlow::Calculate_Relative_Permeability_Upwind_Flux_( const Teuchos::RCP<State>& S,
                                                                 const CompositeVector & flux, 
                                                                 const CompositeVector & rel_perm_cells,
                                                                 CompositeVector & rel_perm_faces ) {
  AmanziMesh::Entity_ID_List faces;
  std::vector<int> dirs;

  rel_perm_faces.PutScalar(0.0);
  int c_owned = S->mesh()->count_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  for (int c=0; c!=c_owned; ++c) {
    S->mesh()->cell_get_faces_and_dirs(c, &faces, &dirs);

    for (int n=0; n!=faces.size(); ++n) {
      int f = faces[n];
      if (flux(n) * dirs[n] >= 0.0) {
        rel_perm_faces(f) = rel_perm_cells(c);
      } else if (bc_markers_[f] != Operators::MFD_BC_NULL) {
        rel_perm_faces(f) = rel_perm_cells(c);
      }
    }
  }
}


// /* ******************************************************************
// * Defines relative permeabilities for faces via arithmetic averaging.
// ****************************************************************** */
// void OverlandFlow::Calculate_Relative_Permeability_Arithmetic_Mean_( const Teuchos::RCP<State>& S,
//                                                                      const CompositeVector & rel_perm_cells,
//                                                                      CompositeVector & rel_perm_faces ) {
//   AmanziMesh::Entity_ID_List cells;
  
//   rel_perm_faces.PutScalar(0.0);
//   int f_owned = S->mesh()->count_entities(AmanziMesh::FACE, AmanziMesh::OWNED);
//   for (int f=0; f!=f_owned; ++f) {
//     S->mesh()->face_get_cells(f, AmanziMesh::USED, &cells);
//     for (int n=0; n!=cells.size(); ++n) { 
//       rel_perm_faces(f) += rel_perm_cells(cells[n]);
//     }
//     rel_perm_faces(f) /= cells.size();
//   }
// }

void OverlandFlow::DeriveFaceValuesFromCellValues_(const Teuchos::RCP<State>& S,
                                                   const Teuchos::RCP<CompositeVector> & pres) {
  AmanziMesh::Entity_ID_List cells;
  
  int f_owned = S->mesh()->count_entities(AmanziMesh::FACE, AmanziMesh::OWNED);
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
    bc_values_[f] = pres("cell",0,cells[0]) ;
  }

  // for (bc=bc_head_->begin(); bc!=bc_head_->end(); ++bc) {
  //   int f = bc->first;
  //   bc_markers_[f] = Operators::MFD_BC_DIRICHLET;
  //   bc_values_[f] = bc->second;
  // }

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
  return 0.; //slope_x[0]*( Lx - x ) ;
}

void OverlandFlow::TestOneSetElevation_( const Teuchos::RCP<State>& S ) {

  CompositeVector & elev = *(S ->GetFieldData("elevation","overland_flow"));
  
  // add elevation to pres_elev, cell dofs
  int c_owned = S->mesh()->count_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  for (int c=0; c!=c_owned; ++c) {
    Amanzi::AmanziGeometry::Point pc = S->mesh()->cell_centroid( c ) ;
    elev("cell",0,c) = TestOneElevation_( pc[0], pc[1] ) ;
  }  

  // add elevation to pres_elev, cell dofs
  int f_owned = S->mesh()->count_entities(AmanziMesh::FACE, AmanziMesh::OWNED);
  for (int f=0; f!=f_owned; ++f) {
    Amanzi::AmanziGeometry::Point pf = S->mesh()->face_centroid( f ) ;    
    elev("face",0,f) = TestOneElevation_( pf[0], pf[1] ) ;
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
  int c_owned = S->mesh()->count_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  for (int c=0; c!=c_owned; ++c) {
    Amanzi::AmanziGeometry::Point pc = S->mesh()->cell_centroid( c ) ;
    elev("cell",0,c) = TestTwoElevation_( pc[0], pc[1] ) ;
  }  

  // add elevation to pres_elev, cell dofs
  int f_owned = S->mesh()->count_entities(AmanziMesh::FACE, AmanziMesh::OWNED);
  for (int f=0; f!=f_owned; ++f) {
    Amanzi::AmanziGeometry::Point pf = S->mesh()->face_centroid( f ) ;    
    elev("face",0,f) = TestTwoElevation_( pf[0], pf[1] ) ;
  }  
}

void OverlandFlow::SetUpElevation_( const Teuchos::RCP<State>& S ) {

#if 1

  // TEST 1
  load_value = 2. * 0.0254 / 3600. ;
  t_rain     = 1800.  ;
  TestOneSetUpElevationPars_() ;
  TestOneSetElevation_( S ) ;

#else
  
  // TEST 2
  load_value = 3.e-6 ; // [m/s] == 10.8e-3 / 3600.
  t_rain     = 5400.  ;     // [s], duration of raining event

  TestTwoSetUpElevationPars_() ;
  TestTwoSetElevation_( S ) ;

#endif

}  


// loading term for the raining events
double OverlandFlow::rhs_load_value( double t ) {
  return t<t_rain ? load_value : 0. ; 
}

} // namespace
} // namespace
