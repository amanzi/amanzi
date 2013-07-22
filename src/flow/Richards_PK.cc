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
#include "BoundaryFunction.hh"

#include "Mesh.hh"
#include "Point.hh"
#include "mfd3d_diffusion.hh"

#include "Matrix_MFD_PLambda.hh"
#include "Matrix_MFD_TPFA.hh"

#include "Flow_BC_Factory.hh"
#include "Flow_State.hh"
#include "Richards_PK.hh"


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
  matrix_ = NULL;
  preconditioner_ = NULL;

  bdf2_dae = NULL;
  bdf1_dae = NULL;

  Flow_PK::Init(global_list, FS_MPC);
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

  mesh_ = FS->mesh();
  dim = mesh_->space_dimension();

  // Create the combined cell/face DoF map.
  super_map_ = CreateSuperMap();

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

  verbosity = FLOW_VERBOSITY_NONE;
  src_sink = NULL;
  src_sink_distribution = 0;
  experimental_solver_ = FLOW_SOLVER_NKA;
}


/* ******************************************************************
* Clean memory.
****************************************************************** */
Richards_PK::~Richards_PK()
{
  delete super_map_;
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

  // Create water retention models
  Teuchos::ParameterList& wrm_list = rp_list_.sublist("Water retention models");
  rel_perm = Teuchos::rcp(new RelativePermeability(mesh_, wrm_list));
  rel_perm->Init(atm_pressure, FS_aux);
  rel_perm->ProcessParameterList();
  rel_perm->PopulateMapC2MB();
  rel_perm->set_experimental_solver(experimental_solver_);

  std::string krel_method_name = rp_list_.get<string>("relative permeability");
  rel_perm->ProcessStringRelativePermeability(krel_method_name);

  // Select a proper matrix class
  if (experimental_solver_ == FLOW_SOLVER_PICARD_NEWTON) {
    matrix_ = new Matrix_MFD_PLambda(FS, *super_map_);
    preconditioner_ = new Matrix_MFD_PLambda(FS, *super_map_);
  } else if (experimental_solver_ == FLOW_SOLVER_NEWTON) {
    matrix_ = new Matrix_MFD_TPFA(FS, *super_map_);
    preconditioner_ = new Matrix_MFD_TPFA(FS, *super_map_);
    matrix_->AddActionProperty(AmanziFlow::FLOW_MATRIX_ACTION_PRECONDITIONER);
    preconditioner_->AddActionProperty(AmanziFlow::FLOW_MATRIX_ACTION_MATRIX);
  } else {
    matrix_ = new Matrix_MFD(FS, *super_map_);
    preconditioner_ = new Matrix_MFD(FS, *super_map_);
  }
  matrix_->AddActionProperty(AmanziFlow::FLOW_MATRIX_ACTION_MATRIX);
  preconditioner_->AddActionProperty(AmanziFlow::FLOW_MATRIX_ACTION_PRECONDITIONER);

  // Create the solution (pressure) vector.
  solution = Teuchos::rcp(new Epetra_Vector(*super_map_));
  solution_cells = Teuchos::rcp(FS->CreateCellView(*solution));
  solution_faces = Teuchos::rcp(FS->CreateFaceView(*solution));
  rhs = Teuchos::rcp(new Epetra_Vector(*super_map_));
  rhs = matrix_->rhs();  // import rhs from the matrix

  const Epetra_BlockMap& cmap = mesh_->cell_map(false);
  const Epetra_BlockMap& fmap_ghost = mesh_->face_map(true);

  pdot_cells_prev = Teuchos::rcp(new Epetra_Vector(cmap));
  pdot_cells = Teuchos::rcp(new Epetra_Vector(cmap));

  // Initialize times.
  double time = FS->get_time();
  if (time >= 0.0) T_physics = time;

  // Initialize actions on boundary condtions. 
  ProcessShiftWaterTableList(rp_list_, bc_head, shift_water_table_);

  // Process other fundamental structures.
  K.resize(ncells_wghost);
  SetAbsolutePermeabilityTensor(K);
  rel_perm->CalculateKVectorUnit(K, gravity_);

  is_matrix_symmetric = SetSymmetryProperty();
  matrix_->SetSymmetryProperty(is_matrix_symmetric);
  matrix_->SymbolicAssembleGlobalMatrices(*super_map_);

  /// Initialize Transmisibillities and Gravity terms contribution
  if (experimental_solver_ == FLOW_SOLVER_NEWTON) {
    Transmis_faces = Teuchos::rcp(new Epetra_Vector(fmap_ghost));
    Grav_term_faces = Teuchos::rcp(new Epetra_Vector(fmap_ghost));
    ComputeTransmissibilities( *Transmis_faces, *Grav_term_faces);
  }

  // Allocate memory for wells
  if (src_sink_distribution & Amanzi::DOMAIN_FUNCTION_ACTION_DISTRIBUTE_PERMEABILITY) {
    const Epetra_Map& cmap_owned = mesh_->cell_map(false);
    Kxy = Teuchos::rcp(new Epetra_Vector(cmap_owned));
  }

  // initialize boundary and source data 
  Epetra_Vector& pressure = FS->ref_pressure();
  Epetra_Vector& lambda = FS->ref_lambda();
  UpdateSourceBoundaryData(time, pressure, lambda);

  // injected water mass
  mass_bc = 0.0;

  // CPU statistcs
  if (verbosity >= FLOW_VERBOSITY_HIGH) {
    timer.add("Mass matrix generation", Amanzi::Timer::ACCUMULATE);
    timer.add("AztecOO solver", Amanzi::Timer::ACCUMULATE);
  }

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
  UpdateSourceBoundaryData(time, pressure, lambda);

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
  if (MyPID == 0 && verbosity >= FLOW_VERBOSITY_MEDIUM) {
    std::printf("Flow PK: initializing with a saturated steady state...\n");
  }
  double T = FS->get_time();
  SolveFullySaturatedProblem(T, *solution, ti_specs->ls_specs_ini);

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

  PrintStatisticsCPU();
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
void Richards_PK::InitNextTI(double T0, double dT0, TI_Specs& ti_specs)
{
  ResetPKtimes(T0, dT0);

  if (MyPID == 0 && verbosity >= FLOW_VERBOSITY_MEDIUM) {
    LinearSolver_Specs& ls = ti_specs.ls_specs;
    std::printf("*************************************************************\n");
    std::printf("Flow PK: TI phase: \"%s\"\n", ti_specs.ti_method_name.c_str());
    std::printf("%5s starts at T=%9.4e [y] with dT=%9.4e [sec]\n", "", T0 / FLOW_YEAR, dT0);
    std::printf("%5s time stepping strategy id %2d\n", "", ti_specs.dT_method);
    std::printf("%5s source/sink distribution method id %2d\n", "", src_sink_distribution);
    std::printf("%5s error control options: %X\n", "", error_control_);
    std::printf("%5s linear solver criteria: ||r||< %9.3e  #itr < %d\n",
        "", ls.convergence_tol, ls.max_itrs);
    std::printf("%7s iterative method AztecOO id %d (gmres=1)\n", "", ls.method);
    std::printf("%7s preconditioner: \"%s\"\n", " ", ti_specs.preconditioner_name.c_str());

    if (ti_specs.initialize_with_darcy) {
      LinearSolver_Specs& ls_ini = ti_specs.ls_specs_ini;
      std::printf("%5s pressure re-initialization (saturated solution)\n", "");
      std::printf("%7s linear solver criteria: ||r||< %9.3e  #itr < %d\n",
          "", ls_ini.convergence_tol, ls_ini.max_itrs);
      std::printf("%7s iterative method AztecOO id %d (cg=0)\n", "", ls_ini.method);
      std::printf("%7s preconditioner: \"%s\"\n", " ", ls_ini.preconditioner_name.c_str());
      if (ti_specs.clip_saturation > 0.0) {
        std::printf("%7s clipping saturation value: %9.4g [-]\n", "", ti_specs.clip_saturation);
      } else if (ti_specs.clip_pressure > -5 * FLOW_PRESSURE_ATMOSPHERIC) {
        std::printf("%7s clipping pressure value: %9.4g [Pa]\n", "", ti_specs.clip_pressure);
      }
    }
  }

  // set up new preconditioner
  int method = ti_specs.preconditioner_method;
  Teuchos::ParameterList& tmp_list = preconditioner_list_.sublist(ti_specs.preconditioner_name);
  Teuchos::ParameterList prec_list;
  if (method == FLOW_PRECONDITIONER_TRILINOS_ML) {
    prec_list = tmp_list.sublist("ML Parameters"); 
  } else if (method == FLOW_PRECONDITIONER_HYPRE_AMG) {
    prec_list = tmp_list.sublist("BoomerAMG Parameters"); 
  } else if (method == FLOW_PRECONDITIONER_TRILINOS_BLOCK_ILU) {
    prec_list = tmp_list.sublist("Block ILU Parameters");
  }

  string mfd3d_method_name = tmp_list.get<string>("discretization method", "monotone mfd");
  ProcessStringMFD3D(mfd3d_method_name, &mfd3d_method_preconditioner_); 

  preconditioner_->SetSymmetryProperty(is_matrix_symmetric);
  preconditioner_->SymbolicAssembleGlobalMatrices(*super_map_);
  preconditioner_->InitPreconditioner(method, prec_list);

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
  }

  // initialize mass matrices
  SetAbsolutePermeabilityTensor(K);
  for (int c = 0; c < ncells_wghost; c++) K[c] *= rho_ / mu_;

  if (experimental_solver_ != FLOW_SOLVER_NEWTON) {
    if (verbosity >= FLOW_VERBOSITY_HIGH) timer.start("Mass matrix generation");  
    matrix_->CreateMFDmassMatrices(mfd3d_method_, K);
    if (verbosity >= FLOW_VERBOSITY_HIGH) timer.stop("Mass matrix generation");  
    preconditioner_->CreateMFDmassMatrices(mfd3d_method_preconditioner_, K);
  }

  if (verbosity >= FLOW_VERBOSITY_MEDIUM) {
    int missed_tmp = missed_bc_faces_;
#ifdef HAVE_MPI
    mesh_->get_comm()->SumAll(&missed_tmp, &missed_bc_faces_, 1);
#endif

    if (MyPID == 0) {
      int nokay = matrix_->nokay();
      int npassed = matrix_->npassed();

      std::printf("%5s discretization method (prec): \"%s\"\n", "", mfd3d_method_name.c_str());
      std::printf("%7s assign default zero flux BC to %d faces\n", "", missed_bc_faces_);
      std::printf("%7s successful and passed elemental matrices: %d %d\n", "", nokay, npassed);   
      std::printf("*************************************************************\n");
    }
  }

  // Well modeling
  if (src_sink_distribution & Amanzi::DOMAIN_FUNCTION_ACTION_DISTRIBUTE_PERMEABILITY) {
    CalculatePermeabilityFactorInWell(K, *Kxy);
  }

  // initialize pressure and lambda
  Epetra_Vector& pressure = FS->ref_pressure();
  Epetra_Vector& lambda = FS->ref_lambda();

  *solution_cells = pressure;
  *solution_faces = lambda;

  if (ti_specs.initialize_with_darcy) {
    // Get a hydrostatic solution consistent with b.c.
    SolveFullySaturatedProblem(T0, *solution, ti_specs.ls_specs_ini);

    if (ti_specs.clip_saturation > 0.0) {
      double pmin = FLOW_PRESSURE_ATMOSPHERIC;
      ClipHydrostaticPressure(pmin, ti_specs.clip_saturation, *solution_cells);
      DeriveFaceValuesFromCellValues(*solution_cells, *solution_faces);
    } else if (ti_specs.clip_pressure > -5 * FLOW_PRESSURE_ATMOSPHERIC) {
      ClipHydrostaticPressure(ti_specs.clip_pressure, *solution_cells);
      DeriveFaceValuesFromCellValues(*solution_cells, *solution_faces);
    }
    pressure = *solution_cells;
  }

  // initialize saturation
  Epetra_Vector& ws = FS->ref_water_saturation();
  DeriveSaturationFromPressure(pressure, ws);
 
  // re-initialize lambda (experimental)
  if (ti_specs.pressure_lambda_constraints && experimental_solver_ == FLOW_SOLVER_NKA) {
    double Tp = T0 + dT0;
    EnforceConstraints_MFD(Tp, *solution);
  }

  // nonlinear solver control options
  ti_specs.num_itrs = 0;
  block_picard = 0;

  Epetra_Vector& Krel_faces = rel_perm->Krel_faces();

  Epetra_Vector& flux = FS->ref_darcy_flux();
  if (experimental_solver_ != FLOW_SOLVER_NEWTON) {
    matrix_->CreateMFDstiffnessMatrices(*rel_perm);  // We remove dT from mass matrices.
    matrix_->DeriveDarcyMassFlux(*solution, *face_importer_, flux);

    AddGravityFluxes_DarcyFlux(K, flux, *rel_perm);
  } else {
    Matrix_MFD_TPFA* matrix_tpfa = dynamic_cast<Matrix_MFD_TPFA*>(matrix_);
    if (matrix_tpfa == 0) {
      Errors::Message msg;
      msg << "Flow PK: cannot cast pointer to class Matrix_MFD_TPFA\n";
      Exceptions::amanzi_throw(msg);
    }
    matrix_tpfa->DeriveDarcyMassFlux(
        *solution, Krel_faces, *Transmis_faces, *Grav_term_faces, bc_model, bc_values, flux);
  }

  for (int f = 0; f < nfaces_owned; f++) flux[f] /= rho_;
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
  // update the axiliary flow state
  FS->CopyMasterFace2GhostFace(FS->ref_darcy_flux(), FS_aux->ref_darcy_flux());

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
        AdvanceToSteadyState(time, dT_MPC);
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
    dTnext = std::min(dT_MPC * dTfactor, ti_specs->dTmax);
  }

  dt_tuple times(time, dT_MPC);
  ti_specs->dT_history.push_back(times);

  ti_specs->num_itrs++;
  if (ti_specs->num_itrs % 100 == 0) PrintStatisticsCPU();

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

  Epetra_Vector& ws = FS_MPC->ref_water_saturation();
  Epetra_Vector& ws_prev = FS_MPC->ref_prev_water_saturation();
  Epetra_Vector& Krel_faces = rel_perm->Krel_faces();

  ws_prev = ws;
  DeriveSaturationFromPressure(pressure, ws);

  // calculate Darcy flux as diffusive part + advective part.
  Epetra_Vector& flux = FS_MPC->ref_darcy_flux();
  if (experimental_solver_ != FLOW_SOLVER_NEWTON) {
    matrix_->CreateMFDstiffnessMatrices(*rel_perm);  // We remove dT from mass matrices.
    matrix_->DeriveDarcyMassFlux(*solution, *face_importer_, flux);

    AddGravityFluxes_DarcyFlux(K, flux, *rel_perm);
  } else {
    Matrix_MFD_TPFA* matrix_tpfa = dynamic_cast<Matrix_MFD_TPFA*>(matrix_);
    if (matrix_tpfa == 0) {
      Errors::Message msg;
      msg << "Flow PK: cannot cast pointer to class Matrix_MFD_TPFA\n";
      Exceptions::amanzi_throw(msg);
    }
    matrix_tpfa->DeriveDarcyMassFlux(
        *solution, Krel_faces, *Transmis_faces, *Grav_term_faces, bc_model, bc_values, flux);
  }

  for (int f = 0; f < nfaces_owned; f++) flux[f] /= rho_;

  // update time derivative
  *pdot_cells_prev = *pdot_cells;

  // update mass balance
  // ImproveAlgebraicConsistency(flux, ws_prev, ws);
  
  if (verbosity >= FLOW_VERBOSITY_HIGH && flow_status_ >= FLOW_STATUS_TRANSIENT_STATE) {
    Epetra_Vector& phi = FS_MPC->ref_porosity();
    double mass_bc_dT = WaterVolumeChangePerSecond(bc_model, flux) * rho_ * dT;

    mass_amanzi = 0.0;
    for (int c = 0; c < ncells_owned; c++) {
      mass_amanzi += ws[c] * rho_ * phi[c] * mesh_->cell_volume(c);
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
  rel_perm->DerivedSdP(u, dSdP);
 
  for (int c = 0; c < ncells_owned; c++) {
    double volume = mesh_->cell_volume(c);
    if (dSdP[c] > 0.0) (*udot_cells)[c] /= volume * dSdP[c] * phi[c] * rho_;
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
  rel_perm->DerivedSdP(pressure_cells, dSdP);

  const Epetra_Vector& phi = FS->ref_porosity();
  std::vector<double>& Acc_cells = matrix_operator->Acc_cells();
  std::vector<double>& Fc_cells = matrix_operator->Fc_cells();

  int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);

  for (int c = 0; c < ncells; c++) {
    double volume = mesh_->cell_volume(c);
    double factor = rho_ * phi[c] * dSdP[c] * volume / dT_prec;
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

    dTfactor_cell = sqrt(tol / std::max(error, FLOW_DT_ADAPTIVE_ERROR_TOLERANCE));
    *dTfactor = std::min(*dTfactor, dTfactor_cell);

    error_max = std::max(error_max, error - tol);
  }

  *dTfactor *= FLOW_DT_ADAPTIVE_SAFETY_FACTOR;
  *dTfactor = std::min(*dTfactor, FLOW_DT_ADAPTIVE_INCREASE);
  *dTfactor = std::max(*dTfactor, FLOW_DT_ADAPTIVE_REDUCTION);

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
  int method = rel_perm->method();
  bool sym = false;
  if ((method == FLOW_RELATIVE_PERM_CENTERED) ||
      (method == FLOW_RELATIVE_PERM_AMANZI)) sym = true;
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

  WhetStone::MFD3D_Diffusion mfd(mesh_);
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
        wsmin = std::min(wsmin, ws[c2]);
        wsmax = std::max(wsmax, ws[c2]);
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
    ws[c] = std::max(ws[c], wsmin);
    ws[c] = std::min(ws[c], wsmax);
  }
}


}  // namespace AmanziFlow
}  // namespace Amanzi

