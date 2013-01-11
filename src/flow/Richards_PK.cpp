/*
This is the flow component of the Amanzi code. 

Copyright 2010-2012 held jointly by LANS/LANL, LBNL, and PNNL. 
Amanzi is released under the three-clause BSD License. 
The terms of use and "as is" disclaimer for this license are 
provided Reconstruction.cppin the top-level COPYRIGHT file.

Authors: Neil Carlson (nnc@lanl.gov), 
         Konstantin Lipnikov (lipnikov@lanl.gov)

Usage: Richards_PK FPK(parameter_list, flow_state);
       FPK->InitPK();
       FPK->Initialize();  // optional
       FPK->InitSteadyState(T, dT);
       FPK->InitTransient(T, dT);
*/

#include <vector>

#include "Epetra_IntVector.h"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

#include "dbc.hh"
#include "exceptions.hh"

#include "Mesh.hh"
#include "Point.hh"
// #include "Matrix_Audit.hpp"
#include "Matrix_MFD_PLambda.hpp"
#include "Matrix_MFD_TPFA.hpp"

#include "Flow_BC_Factory.hpp"
#include "Flow_State.hpp"
#include "boundary_function.hh"
#include "Richards_PK.hpp"


namespace Amanzi {
namespace AmanziFlow {

/* ******************************************************************
* We set up only default values and call Init() routine to complete
* each variable initialization
****************************************************************** */
Richards_PK::Richards_PK(Teuchos::ParameterList& global_list, Teuchos::RCP<Flow_State> FS_MPC)
{
  // initialize pointers (Do we need smart pointers here? lipnikov@lanl.gov)
  bc_pressure = NULL;
  bc_head = NULL;
  bc_flux = NULL;
  bc_seepage = NULL;

  super_map_ = NULL;
  solver = NULL;
  matrix_ = NULL;
  preconditioner_ = NULL;

  bdf2_dae = NULL;
  bdf1_dae = NULL;

  Flow_PK::Init(FS_MPC);
  FS = FS_MPC;

  // extract two critical sublists
  Teuchos::ParameterList flow_list;
  if (global_list.isSublist("Flow")) {
    flow_list = global_list.sublist("Flow");
  } else {
    Errors::Message msg("Flow PK: input parameter list does not have <Flow> sublist.");
    Exceptions::amanzi_throw(msg);
  }

  if (flow_list.isSublist("Richards Problem")) {
    rp_list_ = flow_list.sublist("Richards Problem");
  } else {
    Errors::Message msg("Flow PK: input parameter list does not have <Richards Problem> sublist.");
    Exceptions::amanzi_throw(msg);
  }

  if (global_list.isSublist("Preconditioners")) {
    preconditioner_list_ = global_list.sublist("Preconditioners");
  } else {
    Errors::Message msg("Flow PK: input parameter list does not have <Preconditioners> sublist.");
    Exceptions::amanzi_throw(msg);
  }

  if (global_list.isSublist("Solvers")) {
    solver_list_ = global_list.sublist("Solvers");
  } else {
    Errors::Message msg("Flow PK: input parameter list does not have <Solvers> sublist.");
    Exceptions::amanzi_throw(msg);
  }

  mesh_ = FS->mesh();
  dim = mesh_->space_dimension();

  // Create the combined cell/face DoF map.
  super_map_ = CreateSuperMap();

  // Other fundamental physical quantaties
  rho = *(FS->fluid_density());
  mu = *(FS->fluid_viscosity());
  gravity_.init(dim);
  for (int k = 0; k < dim; k++) gravity_[k] = (*(FS->gravity()))[k];

#ifdef HAVE_MPI
  const Epetra_Comm& comm = mesh_->cell_map(false).Comm();
  MyPID = comm.MyPID();

  const Epetra_Map& source_cmap = mesh_->cell_map(false);
  const Epetra_Map& target_cmap = mesh_->cell_map(true);

  cell_importer_ = Teuchos::rcp(new Epetra_Import(target_cmap, source_cmap));

  const Epetra_Map& source_fmap = mesh_->face_map(false);
  const Epetra_Map& target_fmap = mesh_->face_map(true);

  face_importer_ = Teuchos::rcp(new Epetra_Import(target_fmap, source_fmap));
#endif

  // miscalleneous default parameters
  ti_specs = NULL;
  block_picard = 1;
  error_control_ = FLOW_TI_ERROR_CONTROL_PRESSURE;

  mfd3d_method_ = FLOW_MFD3D_OPTIMIZED;
  mfd3d_method_preconditioner_ = FLOW_MFD3D_OPTIMIZED;

  Krel_method = FLOW_RELATIVE_PERM_UPWIND_GRAVITY;

  verbosity = FLOW_VERBOSITY_HIGH;
  src_sink_distribution = FLOW_SOURCE_DISTRIBUTION_NONE;
  internal_tests = 0;
  experimental_solver_ = FLOW_SOLVER_NKA;
  IsPureNewton = false;
}


/* ******************************************************************
* Clean memory.
****************************************************************** */
Richards_PK::~Richards_PK()
{
  delete super_map_;
  if (solver) delete solver;
  delete matrix_;
  delete preconditioner_;

  delete bdf2_dae;
  delete bc_pressure;
  delete bc_flux;
  delete bc_head;
  delete bc_seepage;

  if (ti_specs != NULL) OutputTimeHistory(ti_specs->dT_history);
}


/* ******************************************************************
* Extract information from Richards Problem parameter list.
****************************************************************** */
void Richards_PK::InitPK()
{
  // Allocate memory for boundary data. It must go first.
  bc_tuple zero = {0.0, 0.0};
  bc_values.resize(nfaces_wghost, zero);
  bc_model.resize(nfaces_wghost, 0);
  bc_submodel.resize(nfaces_wghost, 0);

  rainfall_factor.resize(nfaces_owned, 1.0);

  // Read flow list and populate various structures. 
  ProcessParameterList();

  // Select a proper matrix class
  if (experimental_solver_ == FLOW_SOLVER_PICARD_NEWTON) {
    matrix_ = new Matrix_MFD_PLambda(FS, *super_map_);
    preconditioner_ = new Matrix_MFD_PLambda(FS, *super_map_);
  } else if (experimental_solver_ == FLOW_SOLVER_NEWTON) {
    matrix_ = new Matrix_MFD_TPFA(FS, *super_map_);
    preconditioner_ = new Matrix_MFD_TPFA(FS, *super_map_);
    IsPureNewton = true;
  } else {
    matrix_ = new Matrix_MFD(FS, *super_map_);
    preconditioner_ = new Matrix_MFD(FS, *super_map_);
  }

  // Create the solution (pressure) vector.
  solution = Teuchos::rcp(new Epetra_Vector(*super_map_));
  solution_cells = Teuchos::rcp(FS->CreateCellView(*solution));
  solution_faces = Teuchos::rcp(FS->CreateFaceView(*solution));
  rhs = Teuchos::rcp(new Epetra_Vector(*super_map_));
  rhs = matrix_->rhs();  // import rhs from the matrix

  const Epetra_BlockMap& cmap = mesh_->cell_map(false);
  pdot_cells_prev = Teuchos::rcp(new Epetra_Vector(cmap));
  pdot_cells = Teuchos::rcp(new Epetra_Vector(cmap));

  // Initialize times.
  double time = FS->get_time();
  if (time >= 0.0) T_physics = time;

  // Process other fundamental structures.
  K.resize(ncells_wghost);
  is_matrix_symmetric = SetSymmetryProperty();
  matrix_->SetSymmetryProperty(is_matrix_symmetric);
  matrix_->SymbolicAssembleGlobalMatrices(*super_map_);

  // Allocate data for relative permeability
  const Epetra_Map& cmap_wghost = mesh_->cell_map(true);
  const Epetra_Map& fmap_wghost = mesh_->face_map(true);

  Krel_cells = Teuchos::rcp(new Epetra_Vector(cmap_wghost));
  Krel_faces = Teuchos::rcp(new Epetra_Vector(fmap_wghost));

  dKdP_cells = Teuchos::rcp(new Epetra_Vector(cmap_wghost));  // for P-N and Newton
  dKdP_faces = Teuchos::rcp(new Epetra_Vector(fmap_wghost));

  Krel_cells->PutScalar(1.0);  // we start with fully saturated media
  Krel_faces->PutScalar(1.0);

  if (Krel_method == FLOW_RELATIVE_PERM_UPWIND_GRAVITY || 
      Krel_method == FLOW_RELATIVE_PERM_EXPERIMENTAL) {
    SetAbsolutePermeabilityTensor(K);
    CalculateKVectorUnit(gravity_, Kgravity_unit);
  }

  // Allocate memory for wells
  if (src_sink_distribution == FLOW_SOURCE_DISTRIBUTION_PERMEABILITY) {
    Kxy = Teuchos::rcp(new Epetra_Vector(mesh_->cell_map(false)));
  }

  // injected water mass
  mass_bc = 0.0;

  // miscalleneous maps 
  map_c2mb = Teuchos::rcp(new Epetra_Vector(cmap_wghost));
  PopulateMapC2MB();

  flow_status_ = FLOW_STATUS_INIT;
}


/* ******************************************************************
* Initialization of auxiliary variables (lambda and two saturations).
* WARNING: Flow_PK may use complex initialization of the remaining 
* state variables.
****************************************************************** */
void Richards_PK::InitializeAuxiliaryData()
{
  // pressures
  Epetra_Vector& pressure = FS->ref_pressure();
  Epetra_Vector& lambda = FS->ref_lambda();
  DeriveFaceValuesFromCellValues(pressure, lambda);

  double time = T_physics;
  UpdateSourceBoundaryData(time, lambda);

  // saturations
  Epetra_Vector& ws = FS->ref_water_saturation();
  Epetra_Vector& ws_prev = FS->ref_prev_water_saturation();

  DeriveSaturationFromPressure(pressure, ws);
  ws_prev = ws;

}


/* ******************************************************************
* Initial pressure is set to the pressure for fully saturated rock.
****************************************************************** */
void Richards_PK::InitializeSteadySaturated()
{ 
  double T = FS->get_time();
  SolveFullySaturatedProblem(T, *solution);
}


/* ******************************************************************
* Specific initialization of the initial pressure guess calculation.
****************************************************************** */
void Richards_PK::InitPicard(double T0)
{
  ti_specs = &ti_specs_igs_;
  error_control_ = FLOW_TI_ERROR_CONTROL_PRESSURE;
  error_control_ |= ti_specs->error_control_options;

  InitNextTI(T0, 0.0, ti_specs_igs_);

  // calculate initial guess: cleaning is required (lipnikov@lanl.gov)
  T_physics = ti_specs_igs_.T0;
  dT = ti_specs_igs_.dT0;
  if (experimental_solver_ == FLOW_SOLVER_PICARD_NEWTON) 
    AdvanceToSteadyState_PicardNewton(ti_specs_igs_);
  else
    AdvanceToSteadyState_Picard(ti_specs_igs_);

  Epetra_Vector& ws = FS->ref_water_saturation();
  Epetra_Vector& ws_prev = FS->ref_prev_water_saturation();
  DeriveSaturationFromPressure(*solution_cells, ws);
  ws_prev = ws;

  flow_status_ = FLOW_STATUS_INITIAL_GUESS;
}


/* ******************************************************************
* Specific initialization of a steady state time integration phase.
****************************************************************** */
void Richards_PK::InitSteadyState(double T0, double dT0)
{
  if (ti_specs != NULL) OutputTimeHistory(ti_specs->dT_history);
  ti_specs = &ti_specs_sss_;

  error_control_ = FLOW_TI_ERROR_CONTROL_PRESSURE +  // usually 1 [Pa]
                   FLOW_TI_ERROR_CONTROL_SATURATION;  // usually 1e-4;
  error_control_ |= ti_specs->error_control_options;

  InitNextTI(T0, dT0, ti_specs_sss_);

  flow_status_ = FLOW_STATUS_STEADY_STATE;

  // DEBUG
  // AdvanceToSteadyState();
  // CommitState(FS); WriteGMVfile(FS); exit(0);
}


/* ******************************************************************
* Specific initialization of a transient time integration phase.  
* WARNING: Initialization of lambda is done in MPC and may be 
* erroneous in pure transient mode.
****************************************************************** */
void Richards_PK::InitTransient(double T0, double dT0)
{
  if (ti_specs != NULL) OutputTimeHistory(ti_specs->dT_history);
  ti_specs = &ti_specs_trs_;

  error_control_ = FLOW_TI_ERROR_CONTROL_PRESSURE +  // usually 1 [Pa]
                   FLOW_TI_ERROR_CONTROL_SATURATION;  // usually 1e-4
  error_control_ |= ti_specs->error_control_options;

  InitNextTI(T0, dT0, ti_specs_trs_);

  flow_status_ = FLOW_STATUS_TRANSIENT_STATE;
}


/* ******************************************************************
* Generic initialization of a next time integration phase.
****************************************************************** */
void Richards_PK::InitNextTI(double T0, double dT0, TI_Specs ti_specs)
{
  set_time(T0, dT0);
  bool ini_with_darcy = ti_specs.initialize_with_darcy;

  if (MyPID == 0 && verbosity >= FLOW_VERBOSITY_MEDIUM) {
    std::printf("***********************************************************\n");
    std::printf("Flow PK: TI phase: \"%s\"\n", ti_specs.ti_method_name.c_str());
    std::printf("%5s starts at T=%9.4e [y] with dT=%9.4e [sec]\n", "", T0 / FLOW_YEAR, dT0);
    std::printf("%5s time stepping strategy id %2d\n", "", ti_specs.dT_method);
    std::printf("%5s source/sink distribution method id %2d\n", "", src_sink_distribution);
    std::printf("%5s error control options: %X\n", "", error_control_);
    std::printf("%5s linear solver criteria: ||r||< %9.3e  #itr < %d\n", "", 
        ti_specs.ls_specs.convergence_tol, ti_specs.ls_specs.max_itrs);
    std::printf("%7s preconditioner: \"%s\"\n", " ", ti_specs.preconditioner_name.c_str());

    if (ini_with_darcy) {
      std::printf("%5s initial pressure guess: \"saturated solution\"\n", "");
      if (ti_specs.clip_saturation > 0.0) {
        std::printf("%7s clipping saturation value: %9.4g [-]\n", "", ti_specs.clip_saturation);
      } else if (ti_specs.clip_pressure > -5 * atm_pressure) {
        std::printf("%7s clipping pressure value: %9.4g [Pa]\n", "", ti_specs.clip_pressure);
      }
    }
  }

  // set up new preconditioner
  int method = ti_specs.preconditioner_method;
  Teuchos::ParameterList& tmp_list = preconditioner_list_.sublist(ti_specs.preconditioner_name);
  Teuchos::ParameterList ML_list;
  if (method == FLOW_PRECONDITIONER_TRILINOS_ML) {
    ML_list = tmp_list.sublist("ML Parameters"); 
  } else if (method == FLOW_PRECONDITIONER_HYPRE_AMG) {
    ML_list = tmp_list.sublist("BoomerAMG Parameters"); 
  } else if (method == FLOW_PRECONDITIONER_TRILINOS_BLOCK_ILU) {
    ML_list = tmp_list.sublist("Block ILU Parameters");
  }

  string mfd3d_method_name = tmp_list.get<string>("discretization method", "optimized mfd");
  ProcessStringMFD3D(mfd3d_method_name, &mfd3d_method_preconditioner_); 

  preconditioner_->SetSymmetryProperty(is_matrix_symmetric);
  preconditioner_->SymbolicAssembleGlobalMatrices(*super_map_);
  preconditioner_->InitPreconditioner(method, ML_list);

  // set up new time integration or solver
  std::string ti_method_name(ti_specs.ti_method_name);
  if (ti_specs.ti_method == FLOW_TIME_INTEGRATION_BDF2) {
    Teuchos::ParameterList tmp_list = rp_list_.sublist(ti_method_name).sublist("BDF2").sublist("BDF2 parameters");
    if (!tmp_list.isSublist("VerboseObject"))
        tmp_list.sublist("VerboseObject") = rp_list_.sublist("VerboseObject");

    Teuchos::RCP<Teuchos::ParameterList> bdf2_list(new Teuchos::ParameterList(tmp_list));
    if (bdf2_dae == NULL) bdf2_dae = new BDF2::Dae(*this, *super_map_);
    bdf2_dae->setParameterList(bdf2_list);

  } else if (ti_specs.ti_method == FLOW_TIME_INTEGRATION_BDF1) {
    Teuchos::ParameterList tmp_list = rp_list_.sublist(ti_method_name).sublist("BDF1").sublist("BDF1 parameters");
    if (! tmp_list.isSublist("VerboseObject"))
        tmp_list.sublist("VerboseObject") = rp_list_.sublist("VerboseObject");

    Teuchos::RCP<Teuchos::ParameterList> bdf1_list(new Teuchos::ParameterList(tmp_list));
    if (bdf1_dae == NULL) bdf1_dae = new BDF1Dae(*this, *super_map_);
    bdf1_dae->setParameterList(bdf1_list);

  } else if (solver == NULL) {
    solver = new AztecOO;
    solver->SetUserOperator(matrix_);
    solver->SetPrecOperator(preconditioner_);
    solver->SetAztecOption(AZ_solver, AZ_gmres);
  }

  // initialize mass matrices
  SetAbsolutePermeabilityTensor(K);
  for (int c = 0; c < ncells_wghost; c++) K[c] *= rho / mu;
  matrix_->CreateMFDmassMatrices(mfd3d_method_, K);
  preconditioner_->CreateMFDmassMatrices(mfd3d_method_preconditioner_, K);

  if (MyPID == 0 && verbosity >= FLOW_VERBOSITY_MEDIUM) {
    int nokay = matrix_->nokay();
    int npassed = matrix_->npassed();
    std::printf("%5s discretization method: \"%s\"\n", "", mfd3d_method_name.c_str());
    std::printf("%7s successful and passed elemental matrices: %d %d\n", "", nokay, npassed);   
    std::printf("***********************************************************\n");
  }

  // Well modeling
  if (src_sink_distribution == FLOW_SOURCE_DISTRIBUTION_PERMEABILITY) {
    CalculatePermeabilityFactorInWell(K, *Kxy);
  }

  // linear solver control options
  max_itrs_linear = ti_specs.ls_specs.max_itrs;
  convergence_tol_linear = ti_specs.ls_specs.convergence_tol;

  // initialize pressure and lambda
  Epetra_Vector& pressure = FS->ref_pressure();
  Epetra_Vector& lambda = FS->ref_lambda();

  *solution_cells = pressure;
  *solution_faces = lambda;

  if (ini_with_darcy) {
    SolveFullySaturatedProblem(T0, *solution);  // It gives consistent hydrostatic solution.

    if (ti_specs.clip_saturation > 0.0) {
      double pmin = atm_pressure;
      ClipHydrostaticPressure(pmin, ti_specs.clip_saturation, *solution_cells);
      DeriveFaceValuesFromCellValues(*solution_cells, *solution_faces);
    } else if (ti_specs.clip_pressure > -5 * atm_pressure) {
      ClipHydrostaticPressure(ti_specs.clip_pressure, *solution_cells);
      DeriveFaceValuesFromCellValues(*solution_cells, *solution_faces);
    }
    pressure = *solution_cells;
  }

  // initialize saturation
  Epetra_Vector& ws = FS->ref_water_saturation();
  DeriveSaturationFromPressure(pressure, ws);
 
  // re-initialize lambda (experimental)
  if (ti_specs.pressure_lambda_constraints) {
    double Tp = T0 + dT0;
    EnforceConstraints_MFD(Tp, *solution);
  }

  // nonlinear solver control options
  ti_specs.num_itrs = 0;
  block_picard = 0;
}


/* ******************************************************************
* this routine avoid limitations of MPC by bumping the time step.                                          
****************************************************************** */
double Richards_PK::CalculateFlowDt()
{
  if (ti_specs->ti_method == FLOW_TIME_INTEGRATION_PICARD && block_picard == 1) dT *= 1e+4;
  return dT;
}


/* ******************************************************************* 
* Performs one time step of size dT_MPC either for steady-state or 
* transient calculations.
* Warning: BDF2 and BDF1 will merge eventually.
******************************************************************* */
int Richards_PK::Advance(double dT_MPC)
{
  dT = dT_MPC;
  double time = FS->get_time();
  if (time >= 0.0) T_physics = time;

  // predict water mass change during time step
  time = T_physics;
  if (ti_specs->num_itrs == 0) {  // initialization
    Epetra_Vector udot(*super_map_);
    ComputeUDot(time, *solution, udot);

    if (ti_specs->ti_method == FLOW_TIME_INTEGRATION_BDF2) {
      bdf2_dae->set_initial_state(time, *solution, udot);

    } else if (ti_specs->ti_method == FLOW_TIME_INTEGRATION_BDF1) {
      bdf1_dae->set_initial_state(time, *solution, udot);

    } else if (ti_specs->ti_method == FLOW_TIME_INTEGRATION_PICARD) {
      if (flow_status_ == FLOW_STATUS_STEADY_STATE) {
        AdvanceToSteadyState();
        block_picard = 1;  // We will wait for transient initialization.
      }
    }

    int ierr;
    update_precon(time, *solution, dT, ierr);
    ti_specs->num_itrs++;
  }

  if (ti_specs->ti_method == FLOW_TIME_INTEGRATION_BDF2) {
    bdf2_dae->bdf2_step(dT, 0.0, *solution, dTnext);
    bdf2_dae->commit_solution(dT, *solution);
    bdf2_dae->write_bdf2_stepping_statistics();

    T_physics = bdf2_dae->most_recent_time();

  } else if (ti_specs->ti_method == FLOW_TIME_INTEGRATION_BDF1) {
    bdf1_dae->bdf1_step(dT, *solution, dTnext);
    bdf1_dae->commit_solution(dT, *solution);
    bdf1_dae->write_bdf1_stepping_statistics();

    T_physics = bdf1_dae->most_recent_time();

  } else if (ti_specs->ti_method == FLOW_TIME_INTEGRATION_PICARD) {
    if (block_picard == 0) {
      PicardTimeStep(time, dT, dTnext);  // Updates solution vector.
      //AndersonMixingTimeStep(time, dT, dTnext);
    } else {
      dTnext = dT;
    }
  }

  // Calculate time derivative and 2nd-order solution approximation.
  // Estimate of a time step multiplier overrides the above estimate.
  if (ti_specs->dT_method == FLOW_DT_ADAPTIVE) {
    Epetra_Vector& pressure = FS->ref_pressure();  // pressure at t^n

    for (int c = 0; c < ncells_owned; c++) {
      (*pdot_cells)[c] = ((*solution)[c] - pressure[c]) / dT; 
      (*solution)[c] = pressure[c] + ((*pdot_cells_prev)[c] + (*pdot_cells)[c]) * dT / 2;
    }

    double err, dTfactor;
    err = AdaptiveTimeStepEstimate(&dTfactor);
    if (err > 0.0) throw 1000;  // fix (lipnikov@lan.gov)
    dTnext = std::min<double>(dT_MPC * dTfactor, ti_specs_trs_.dTmax);
  }

  dt_tuple times(time, dT_MPC);
  ti_specs->dT_history.push_back(times);

  ti_specs->num_itrs++;
  return 0;
}


/* ******************************************************************
* Transfer part of the internal data needed by transport to the 
* flow state FS_MPC. MPC may request to populate the original FS.
* The consistency condition is improved by adjusting saturation while
* preserving its LED property.
****************************************************************** */
void Richards_PK::CommitState(Teuchos::RCP<Flow_State> FS_MPC)
{
  // save cell-based and face-based pressures 
  Epetra_Vector& pressure = FS_MPC->ref_pressure();
  pressure = *solution_cells;
  Epetra_Vector& lambda = FS_MPC->ref_lambda();
  lambda = *solution_faces;

  // update saturations
  Epetra_Vector& ws = FS_MPC->ref_water_saturation();
  Epetra_Vector& ws_prev = FS_MPC->ref_prev_water_saturation();

  ws_prev = ws;
  DeriveSaturationFromPressure(pressure, ws);

  // calculate Darcy flux as diffusive part + advective part.
  Epetra_Vector& flux = FS_MPC->ref_darcy_flux();
  matrix_->CreateMFDstiffnessMatrices(*Krel_cells, *Krel_faces, Krel_method);  // We remove dT from mass matrices.
  matrix_->DeriveDarcyMassFlux(*solution, *face_importer_, flux);
  AddGravityFluxes_DarcyFlux(K, *Krel_cells, *Krel_faces, Krel_method, flux);
  for (int c = 0; c < nfaces_owned; c++) flux[c] /= rho;

  // update time derivative
  *pdot_cells_prev = *pdot_cells;

  // update mass balance
  // ImproveAlgebraicConsistency(flux, ws_prev, ws);

  if (verbosity >= FLOW_VERBOSITY_HIGH && flow_status_ >= FLOW_STATUS_TRANSIENT_STATE) {
    Epetra_Vector& phi = FS_MPC->ref_porosity();
    double mass_bc_dT = WaterVolumeChangePerSecond(bc_model, flux) * rho * dT;

    mass_amanzi = 0.0;
    for (int c = 0; c < ncells_owned; c++) {
      mass_amanzi += ws[c] * rho * phi[c] * mesh_->cell_volume(c);
    }

    double mass_amanzi_tmp = mass_amanzi, mass_bc_tmp = mass_bc_dT;
    mesh_->get_comm()->SumAll(&mass_amanzi_tmp, &mass_amanzi, 1);
    mesh_->get_comm()->SumAll(&mass_bc_tmp, &mass_bc_dT, 1);

    mass_bc += mass_bc_dT;
    if (MyPID == 0)
        std::printf("Flow PK: water mass =%10.5e [kg], total boundary flux = %10.5e [kg]\n", mass_amanzi, mass_bc);
  }

  dT = dTnext;
}


/* ******************************************************************
* Estimate du/dt from the pressure equations, a du/dt = g - A*u.
****************************************************************** */
double Richards_PK::ComputeUDot(double T, const Epetra_Vector& u, Epetra_Vector& udot)
{
  AssembleMatrixMFD(u, T);
  double norm_udot = matrix_->ComputeNegativeResidual(u, udot);

  Epetra_Vector* udot_cells = FS->CreateCellView(udot);
  const Epetra_Vector& phi = FS->ref_porosity();
  Epetra_Vector dSdP(mesh_->cell_map(false));
  DerivedSdP(u, dSdP);
 
  for (int c = 0; c < ncells_owned; c++) {
    double volume = mesh_->cell_volume(c);
    if (dSdP[c] > 0.0) (*udot_cells)[c] /= volume * dSdP[c] * phi[c] * rho;
  }

  Epetra_Vector* udot_faces = FS->CreateFaceView(udot);
  DeriveFaceValuesFromCellValues(*udot_cells, *udot_faces);
  // udot_faces->PutScalar(0.0);

  return norm_udot;
}


/* ******************************************************************
* Temporary convertion from double to tensor.                                               
****************************************************************** */
void Richards_PK::SetAbsolutePermeabilityTensor(std::vector<WhetStone::Tensor>& K)
{
  const Epetra_Vector& vertical_permeability = FS->ref_vertical_permeability();
  const Epetra_Vector& horizontal_permeability = FS->ref_horizontal_permeability();

  const Epetra_BlockMap& cmap_wghost = mesh_->cell_map(true);
  Epetra_Vector perm_vert_gh(cmap_wghost);
  Epetra_Vector perm_horz_gh(cmap_wghost);

  FS->CopyMasterCell2GhostCell(vertical_permeability, perm_vert_gh);
  FS->CopyMasterCell2GhostCell(horizontal_permeability, perm_horz_gh);

  for (int c = 0; c < ncells_wghost; c++) {
    if (perm_vert_gh[c] == perm_horz_gh[c]) {
      K[c].init(dim, 1);
      K[c](0, 0) = perm_vert_gh[c];
    } else {
      K[c].init(dim, 2);
      for (int i = 0; i < dim-1; i++) K[c](i, i) = perm_horz_gh[c];
      K[c](dim-1, dim-1) = perm_vert_gh[c];
    }
  }
}


/* ******************************************************************
* Adds time derivative to the cell-based part of MFD algebraic system.                                               
****************************************************************** */
void Richards_PK::AddTimeDerivative_MFD(
    Epetra_Vector& pressure_cells, double dT_prec, Matrix_MFD* matrix_operator)
{
  Epetra_Vector dSdP(mesh_->cell_map(false));
  DerivedSdP(pressure_cells, dSdP);

  const Epetra_Vector& phi = FS->ref_porosity();
  std::vector<double>& Acc_cells = matrix_operator->Acc_cells();
  std::vector<double>& Fc_cells = matrix_operator->Fc_cells();

  int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);

  for (int c = 0; c < ncells; c++) {
    double volume = mesh_->cell_volume(c);
    double factor = rho * phi[c] * dSdP[c] * volume / dT_prec;
    Acc_cells[c] += factor;
    Fc_cells[c] += factor * pressure_cells[c];
  }
}


/* ******************************************************************
* Estimate dT increase factor by comparing the 1st and 2nd order
* time approximations. 
* WARNING: it is implmented for transient phase only.
****************************************************************** */
double Richards_PK::AdaptiveTimeStepEstimate(double* dTfactor)
{
  double tol, error, error_max = 0.0;
  double dTfactor_cell;

  *dTfactor = 100.0;
  for (int c = 0; c < ncells_owned; c++) {
    error = fabs((*pdot_cells)[c] - (*pdot_cells_prev)[c]) * dT / 2;
    tol = ti_specs_trs_.rtol * fabs((*solution)[c]) + ti_specs_trs_.atol;

    dTfactor_cell = sqrt(tol / std::max<double>(error, FLOW_DT_ADAPTIVE_ERROR_TOLERANCE));
    *dTfactor = std::min<double>(*dTfactor, dTfactor_cell);

    error_max = std::max<double>(error_max, error - tol);
  }

  *dTfactor *= FLOW_DT_ADAPTIVE_SAFETY_FACTOR;
  *dTfactor = std::min<double>(*dTfactor, FLOW_DT_ADAPTIVE_INCREASE);
  *dTfactor = std::max<double>(*dTfactor, FLOW_DT_ADAPTIVE_REDUCTION);

#ifdef HAVE_MPI
    double dT_tmp = *dTfactor;
    solution->Comm().MinAll(&dT_tmp, dTfactor, 1);  // find the global minimum
 
    double error_tmp = error_max;
    solution->Comm().MaxAll(&error_tmp, &error_max, 1);  // find the global maximum
#endif
  return error_max;
}


/* ******************************************************************
* Identify type of the matrix.
****************************************************************** */
bool Richards_PK::SetSymmetryProperty()
{
  bool sym = false;
  if ((Krel_method == FLOW_RELATIVE_PERM_CENTERED) ||
      (Krel_method == FLOW_RELATIVE_PERM_EXPERIMENTAL)) sym = true;
  if (experimental_solver_ == FLOW_SOLVER_NEWTON || 
      experimental_solver_ == FLOW_SOLVER_PICARD_NEWTON) sym = false;

  return sym;
}


/* ******************************************************************
* Saturation should be in exact balance with Darcy fluxes in order to
* have local extrema dimishing (LED) property for concentrations. 
* WARNING: we can enforce it strictly only in some cells.
****************************************************************** */
void Richards_PK::ImproveAlgebraicConsistency(const Epetra_Vector& flux, 
                                              const Epetra_Vector& ws_prev, Epetra_Vector& ws)
{
  // create a disctributed flux vector
  Epetra_Vector flux_d(mesh_->face_map(true));
  for (int f = 0; f < nfaces_owned; f++) flux_d[f] = flux[f];
  FS->CopyMasterFace2GhostFace(flux_d);

  // create a distributed saturation vector
  Epetra_Vector ws_d(mesh_->cell_map(true));
  for (int c = 0; c < ncells_owned; c++) ws_d[c] = ws[c];
  FS->CopyMasterCell2GhostCell(ws_d);

  WhetStone::MFD3D mfd(mesh_);
  const Epetra_Vector& phi = FS->ref_porosity();
  AmanziMesh::Entity_ID_List faces;
  std::vector<int> dirs;

  for (int c = 0; c < ncells_owned; c++) {
    mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
    int nfaces = faces.size();

    // calculate min/max values
    double wsmin, wsmax;
    wsmin = wsmax = ws[c];
    for (int n = 0; n < nfaces; n++) {
      int f = faces[n];
      int c2 = mfd.cell_get_face_adj_cell(c, f);
      if (c2 >= 0) {
        wsmin = std::min<double>(wsmin, ws[c2]);
        wsmax = std::max<double>(wsmax, ws[c2]);
      }
    }

    // predict new saturation
    ws[c] = ws_prev[c];
    double factor = dT / (phi[c] * mesh_->cell_volume(c));
    for (int n = 0; n < nfaces; n++) {
      int f = faces[n];
      ws[c] -= factor * flux_d[f] * dirs[n]; 
    }

    // limit new saturation
    ws[c] = std::max<double>(ws[c], wsmin);
    ws[c] = std::min<double>(ws[c], wsmax);
  }
}


}  // namespace AmanziFlow
}  // namespace Amanzi

