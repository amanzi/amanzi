/*
This is the flow component of the Amanzi code. 

Copyright 2010-2012 held jointly by LANS/LANL, LBNL, and PNNL. 
Amanzi is released under the three-clause BSD License. 
The terms of use and "as is" disclaimer for this license are 
provided Reconstruction.cppin the top-level COPYRIGHT file.

Authors: Neil Carlson (nnc@lanl.gov), 
         Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <vector>

#include "Epetra_IntVector.h"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

#include "dbc.hh"
#include "errors.hh"
#include "exceptions.hh"

#include "Mesh.hh"
#include "Point.hh"
#include "Matrix_Audit.hpp"

#include "Flow_BC_Factory.hpp"
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
  Flow_PK::Init(FS_MPC);

  FS = FS_MPC;

  // extract two critical sublists 
  Teuchos::ParameterList flow_list;
  if (global_list.isSublist("Flow")) {
    flow_list = global_list.sublist("Flow");
  } else {
    Errors::Message msg("Richards PK: input parameter list does not have <Flow> sublist.");
    Exceptions::amanzi_throw(msg);
  }

  if (flow_list.isSublist("Richards Problem")) {
    rp_list_ = flow_list.sublist("Richards Problem");
  } else {
    Errors::Message msg("Richards PK: input parameter list does not have <Richards Problem> sublist.");
    Exceptions::amanzi_throw(msg);
  }

  if (global_list.isSublist("Preconditioners")) {
    preconditioner_list_ = global_list.sublist("Preconditioners");
  } else {
    Errors::Message msg("Richards PK: input parameter list does not have <Preconditioners> sublist.");
    Exceptions::amanzi_throw(msg);
  }

  if (global_list.isSublist("Solvers")) {
    solver_list_ = global_list.sublist("Solvers");
  } else {
    Errors::Message msg("Richards PK: input parameter list does not have <Solvers> sublist.");
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

  // miscalleneous
  solver = NULL;
  bdf2_dae = NULL;
  bdf1_dae = NULL;
  block_picard = 1;
  error_control_ = FLOW_TI_ERROR_CONTROL_PRESSURE;

  ti_method_sss = FLOW_TIME_INTEGRATION_BDF1;  // time integration (TI) parameters
  ti_method_trs = FLOW_TIME_INTEGRATION_BDF2;
  num_itrs_trs = 0;

  absolute_tol_sss = absolute_tol_trs = 1.0;
  relative_tol_sss = relative_tol_trs = 1e-5;
  initialize_with_darcy = 0;

  mfd3d_method_ = FLOW_MFD3D_HEXAHEDRA_MONOTONE;  // will be changed (lipnikov@lanl.gov)
  mfd3d_method_preconditioner_ = FLOW_MFD3D_TWO_POINT_FLUX;

  Krel_method = FLOW_RELATIVE_PERM_UPWIND_GRAVITY;

  verbosity = FLOW_VERBOSITY_HIGH;
  internal_tests = 0;

  num_nonlinear_steps = 0;
}


/* ******************************************************************
* Clean memory.
****************************************************************** */
Richards_PK::~Richards_PK()
{
  delete super_map_;
  if (solver) delete solver;
  if (matrix == preconditioner) {
    delete matrix;
  } else {
    delete matrix;
    delete preconditioner;
  }
  delete bdf2_dae;
  delete bc_pressure;
  delete bc_flux;
  delete bc_head;
  delete bc_seepage;
}


/* ******************************************************************
* Extract information from Richards Problem parameter list.
****************************************************************** */
void Richards_PK::InitPK(Matrix_MFD* matrix_, Matrix_MFD* preconditioner_)
{
  if (matrix_ == NULL)
      matrix = new Matrix_MFD(FS, *super_map_);
  else
      matrix = matrix_;

  if (preconditioner_ == NULL)
      preconditioner = new Matrix_MFD(FS, *super_map_);
  else
      preconditioner = preconditioner_;

  // Create the solution (pressure) vector.
  solution = Teuchos::rcp(new Epetra_Vector(*super_map_));
  solution_cells = Teuchos::rcp(FS->CreateCellView(*solution));
  solution_faces = Teuchos::rcp(FS->CreateFaceView(*solution));
  rhs = Teuchos::rcp(new Epetra_Vector(*super_map_));
  rhs = matrix->rhs();  // import rhs from the matrix

  // Get solver parameters from the flow parameter list.
  ProcessParameterList();

  // Process boundary data (state may be incomplete at this moment)
  int nfaces = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::USED);
  bc_markers.resize(nfaces, FLOW_BC_FACE_NULL);
  bc_values.resize(nfaces, 0.0);

  double T_physical = FS->get_time();
  T_internal = (standalone_mode) ? T_internal : T_physical;

  // Process other fundamental structures
  K.resize(ncells_owned);
  is_matrix_symmetric = (Krel_method == FLOW_RELATIVE_PERM_CENTERED);
  matrix->SetSymmetryProperty(is_matrix_symmetric);
  matrix->SymbolicAssembleGlobalMatrices(*super_map_);

  // Allocate data for relative permeability
  const Epetra_Map& cmap = mesh_->cell_map(true);
  const Epetra_Map& fmap = mesh_->face_map(true);

  Krel_cells = Teuchos::rcp(new Epetra_Vector(cmap));
  Krel_faces = Teuchos::rcp(new Epetra_Vector(fmap));

  Krel_cells->PutScalar(1.0);  // we start with fully saturated media
  Krel_faces->PutScalar(1.0);

  if (Krel_method == FLOW_RELATIVE_PERM_UPWIND_GRAVITY) {
    // Kgravity_unit.resize(ncells_wghost);  Resize does not work properly.
    SetAbsolutePermeabilityTensor(K);
    CalculateKVectorUnit(gravity_, Kgravity_unit);
  }

  flow_status_ = FLOW_STATUS_INIT;
}


/* ******************************************************************
* Separate initialization of solver may be required for steady state
* and transient runs. BDF2 and BDF1 will eventually merge but are 
* separated strictly (no code optimization) for the moment.
****************************************************************** */
void Richards_PK::InitSteadyState(double T0, double dT0)
{
  if (MyPID == 0 && verbosity >= FLOW_VERBOSITY_MEDIUM) {
    std::printf("Richards PK: initializing steady-state at T(sec)=%9.4e dT(sec)=%9.4e \n", T0, dT0);
    if (initialize_with_darcy) {
       std::printf("Richards PK: initializing with a clipped Darcy pressure\n");
       std::printf("Richards PK: clipping saturation value =%5.2g\n", clip_saturation);
     }
  }

  // set up new preconditioner
  Teuchos::ParameterList& tmp_list = preconditioner_list_.sublist(preconditioner_name_sss_);
  Teuchos::ParameterList ML_list = tmp_list.sublist("ML Parameters");

  string mfd3d_method_name = tmp_list.get<string>("discretization method", "two point flux approximation");
  ProcessStringMFD3D(mfd3d_method_name, &mfd3d_method_preconditioner_); 

  preconditioner->SetSymmetryProperty(is_matrix_symmetric);
  preconditioner->SymbolicAssembleGlobalMatrices(*super_map_);
  preconditioner->InitML_Preconditioner(ML_list);

  // set up new time integration or solver
  if (ti_method_sss == FLOW_TIME_INTEGRATION_BDF2) {
    Teuchos::ParameterList solver_list = rp_list_.sublist("steady state time integrator").sublist("nonlinear solver BDF2");
    if (solver_list.isSublist("VerboseObject"))
       solver_list.sublist("VerboseObject") = rp_list_.sublist("VerboseObject");

    Teuchos::RCP<Teuchos::ParameterList> bdf2_list(new Teuchos::ParameterList(solver_list));
    if (bdf2_dae == NULL) bdf2_dae = new BDF2::Dae(*this, *super_map_);
    bdf2_dae->setParameterList(bdf2_list);

  } else if (ti_method_sss == FLOW_TIME_INTEGRATION_BDF1) {
    Teuchos::ParameterList solver_list = rp_list_.sublist("steady state time integrator").sublist("nonlinear solver BDF1");
    if (solver_list.isSublist("VerboseObject"))
        solver_list.sublist("VerboseObject") = rp_list_.sublist("VerboseObject");

    Teuchos::RCP<Teuchos::ParameterList> bdf1_list(new Teuchos::ParameterList(solver_list));
    if (bdf1_dae == NULL) bdf1_dae = new BDF1Dae(*this, *super_map_);
    bdf1_dae->setParameterList(bdf1_list);

  } else {
    solver = new AztecOO;
    solver->SetUserOperator(matrix);
    solver->SetPrecOperator(preconditioner);
    solver->SetAztecOption(AZ_solver, AZ_cg);  // symmetry is required
  }

  // initialize mass matrices
  SetAbsolutePermeabilityTensor(K);
  for (int c = 0; c < K.size(); c++) K[c] *= rho / mu;
  matrix->CreateMFDmassMatrices(mfd3d_method_, K);
  preconditioner->CreateMFDmassMatrices(mfd3d_method_preconditioner_, K);

  if (MyPID == 0 && verbosity >= FLOW_VERBOSITY_HIGH) {
    int nokay = matrix->nokay();
    int npassed = matrix->npassed();
    std::printf("Richards PK: successful and passed matrices: %8d %8d\n", nokay, npassed);   
  }

  // (re)initialize pressure and saturation
  Epetra_Vector& pressure = FS->ref_pressure();
  Epetra_Vector& lambda = FS->ref_lambda();
  Epetra_Vector& water_saturation = FS->ref_water_saturation();

  *solution_cells = pressure;
  *solution_faces = lambda;  // useless due to logic below (lipnikov@lanl.gov)

  if (initialize_with_darcy) {
    double pmin = atm_pressure;
    InitializePressureHydrostatic(T0);
    ClipHydrostaticPressure(pmin, clip_saturation, *solution);
    pressure = *solution_cells;
  }

  DeriveFaceValuesFromCellValues(*solution_cells, *solution_faces);
  DeriveSaturationFromPressure(pressure, water_saturation);

  // control options
  absolute_tol = absolute_tol_sss;
  relative_tol = relative_tol_sss;

  set_time(T0, dT0);  // overrides data provided in the input file
  ti_method = ti_method_sss;
  num_itrs = 0;
  block_picard = 0;
  error_control_ = FLOW_TI_ERROR_CONTROL_PRESSURE +  // usually 1 [Pa]
                   FLOW_TI_ERROR_CONTROL_SATURATION;  // usually 1e-4

  flow_status_ = FLOW_STATUS_STEADY_STATE_INIT;

  // DEBUG
  // AdvanceToSteadyState();
  // CommitStateForTransport(FS); CommitState(FS); WriteGMVfile(FS); exit(0);
}


/* ******************************************************************
* Initialization analyzes status of matrix/preconditioner pair. 
* BDF2 and BDF1 will eventually merge but are separated strictly 
* (no code optimization) for the moment.  
* WARNING: Initialization of lambda is done in MPC and may be 
* erroneous in pure transient mode.
****************************************************************** */
void Richards_PK::InitTransient(double T0, double dT0)
{
  if (MyPID == 0 && verbosity >= FLOW_VERBOSITY_MEDIUM) {
     std::printf("Richards PK: initializing transient flow: T(sec)=%9.4e dT(sec)=%9.4e\n", T0, dT0);
  }

  // set up new preconditioner
  Teuchos::ParameterList& tmp_list = preconditioner_list_.sublist(preconditioner_name_trs_);
  Teuchos::ParameterList ML_list = tmp_list.sublist("ML Parameters");

  string mfd3d_method_name = tmp_list.get<string>("discretization method", "two point flux approximation");
  ProcessStringMFD3D(mfd3d_method_name, &mfd3d_method_preconditioner_); 

  preconditioner->SetSymmetryProperty(is_matrix_symmetric);
  preconditioner->SymbolicAssembleGlobalMatrices(*super_map_);
  preconditioner->InitML_Preconditioner(ML_list);

  if (ti_method_trs == FLOW_TIME_INTEGRATION_BDF2) {
    if (bdf2_dae != NULL) delete bdf2_dae;  // The only way to reset BDF2 is to delete it.

    Teuchos::ParameterList solver_list = rp_list_.sublist("transient time integrator").sublist("nonlinear solver BDF2");
    if (solver_list.isSublist("VerboseObject"))
        solver_list.sublist("VerboseObject") = rp_list_.sublist("VerboseObject");

    Teuchos::RCP<Teuchos::ParameterList> bdf2_list(new Teuchos::ParameterList(solver_list));
    bdf2_dae = new BDF2::Dae(*this, *super_map_);
    bdf2_dae->setParameterList(bdf2_list);

  } else if (ti_method_trs == FLOW_TIME_INTEGRATION_BDF1) {
    if (bdf1_dae != NULL) delete bdf1_dae;  // the only way to reset BDF1 is to delete it

    Teuchos::ParameterList solver_list = rp_list_.sublist("transient time integrator").sublist("nonlinear solver BDF1");
    if (solver_list.isSublist("VerboseObject"))
        solver_list.sublist("VerboseObject") = rp_list_.sublist("VerboseObject");

    Teuchos::RCP<Teuchos::ParameterList> bdf1_list(new Teuchos::ParameterList(solver_list));
    bdf1_dae = new BDF1Dae(*this, *super_map_);
    bdf1_dae->setParameterList(bdf1_list);

  } else if (solver == NULL) {
    solver = new AztecOO;
    solver->SetUserOperator(matrix);
    solver->SetPrecOperator(preconditioner);
    solver->SetAztecOption(AZ_solver, AZ_cg);  // symmetry is required
  }

  // initialize mass matrices
  SetAbsolutePermeabilityTensor(K);
  for (int c = 0; c < K.size(); c++) K[c] *= rho / mu;
  matrix->CreateMFDmassMatrices(mfd3d_method_, K);
  preconditioner->CreateMFDmassMatrices(mfd3d_method_preconditioner_, K);

  // initialize pressure and saturation
  Epetra_Vector& pressure = FS->ref_pressure();
  Epetra_Vector& lambda = FS->ref_lambda();
  Epetra_Vector& water_saturation = FS->ref_water_saturation();
  *solution_cells = pressure;
  *solution_faces = lambda;

  //DeriveFaceValuesFromCellValues(*solution_cells, *solution_faces);
  DeriveSaturationFromPressure(pressure, water_saturation);

  // control options
  absolute_tol = absolute_tol_trs;
  relative_tol = relative_tol_trs;

  dT = dT0_trs;
  T_internal = T0_trs;

  set_time(T0, dT0);  // overrides data in input file
  ti_method = ti_method_trs;
  num_itrs = 0;
  block_picard = 0;
  error_control_ = FLOW_TI_ERROR_CONTROL_PRESSURE +  // usually 1 [Pa]
                   FLOW_TI_ERROR_CONTROL_SATURATION;  // usually 1e-4

  flow_status_ = FLOW_STATUS_TRANSIENT_STATE_INIT;
}


/* ******************************************************************
* this routine avoid limitations of MPC by bumping the time step.                                          
****************************************************************** */
double Richards_PK::CalculateFlowDt()
{
  if (ti_method == FLOW_TIME_INTEGRATION_PICARD && block_picard == 1) dT *= 1e+4;
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

  if (num_itrs == 0) {  // set-up internal clock
    double T_physical = FS->get_time();
    T_internal = (standalone_mode) ? T_internal : T_physical;
  }

  double time = T_internal;
  if (num_itrs == 0) {  // initialization
    Epetra_Vector udot(*super_map_);
    ComputeUDot(time, *solution, udot);
    if (ti_method == FLOW_TIME_INTEGRATION_BDF2) {
      bdf2_dae->set_initial_state(time, *solution, udot);

    } else if (ti_method == FLOW_TIME_INTEGRATION_BDF1) {
      bdf1_dae->set_initial_state(time, *solution, udot);

    } else if (ti_method == FLOW_TIME_INTEGRATION_PICARD) {
      if (flow_status_ == FLOW_STATUS_STEADY_STATE_INIT) {
        AdvanceToSteadyState();
        block_picard = 1;  // We will wait for transient initialization.
      }
    }

    int ierr;
    update_precon(time, *solution, dT, ierr);
  }

  if (ti_method == FLOW_TIME_INTEGRATION_BDF2) {
    bdf2_dae->bdf2_step(dT, 0.0, *solution, dTnext);
    bdf2_dae->commit_solution(dT, *solution);
    bdf2_dae->write_bdf2_stepping_statistics();

    T_internal = bdf2_dae->most_recent_time();

  } else if (ti_method == FLOW_TIME_INTEGRATION_BDF1) {
    bdf1_dae->bdf1_step(dT, *solution, dTnext);
    bdf1_dae->commit_solution(dT, *solution);
    bdf1_dae->write_bdf1_stepping_statistics();

    T_internal = bdf1_dae->most_recent_time();

  } else if (ti_method == FLOW_TIME_INTEGRATION_PICARD) {
    if (block_picard == 0) {
      PicardStep(time, dT, dTnext);  // Updates solution vector.
    } else {
      dTnext = dT;
    }
  }
  num_itrs++;

  flow_status_ = FLOW_STATUS_TRANSIENT_STATE_COMPLETE;
  return 0;
}


/* ******************************************************************
* Transfer part of the internal data needed by transport to the 
* flow state FS_MPC. MPC may request to populate the original FS. 
****************************************************************** */
void Richards_PK::CommitStateForTransport(Teuchos::RCP<Flow_State> FS_MPC)
{
  Epetra_Vector& flux = FS_MPC->ref_darcy_flux();
  matrix->CreateMFDstiffnessMatrices(*Krel_cells, *Krel_faces);
  matrix->DeriveDarcyMassFlux(*solution, *face_importer_, flux);
  AddGravityFluxes_DarcyFlux(K, *Krel_cells, *Krel_faces, flux);
  for (int c = 0; c < nfaces_owned; c++) flux[c] /= rho;

  Epetra_Vector& pressure = FS_MPC->ref_pressure();  // save pressure
  pressure = *solution_cells;
  Epetra_Vector& lambda = FS_MPC->ref_lambda(); // save lambda
  lambda = *solution_faces;

  // update saturations
  Epetra_Vector& ws = FS_MPC->ref_water_saturation();
  Epetra_Vector& ws_prev = FS_MPC->ref_prev_water_saturation();
  ws_prev = ws;

  // CalculateConsistentSaturation(flux, ws_prev, ws);
  DeriveSaturationFromPressure(pressure, ws);

  // DEBUG
  // writeGMVfile(FS_MPC);
}


/* ******************************************************************
* Transfer internal data to flow state FS_MPC. MPC may request
* to populate the original state FS. 
****************************************************************** */
void Richards_PK::CommitState(Teuchos::RCP<Flow_State> FS_MPC)
{
  dT = dTnext;

  Epetra_Vector& flux = FS_MPC->ref_darcy_flux();
  Epetra_MultiVector& velocity = FS_MPC->ref_darcy_velocity();
  DeriveDarcyVelocity(flux, velocity);
}


/* ******************************************************************
* Saturation should be in exact balance with Darcy fluxes in order to
* have extrema dimishing property for concentrations. 
****************************************************************** */
void Richards_PK::CalculateConsistentSaturation(const Epetra_Vector& flux, 
                                                const Epetra_Vector& ws_prev, Epetra_Vector& ws)
{
  // create a disctributed flux vector
  Epetra_Vector flux_d(mesh_->face_map(true));
  for (int f = 0; f < nfaces_owned; f++) flux_d[f] = flux[f];
  FS->CombineGhostFace2MasterFace(flux_d);

  const Epetra_Vector& phi = FS->ref_porosity();
  AmanziMesh::Entity_ID_List faces;
  std::vector<int> dirs;

  for (int c = 0; c < ncells_owned; c++) {
    mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
    int nfaces = faces.size();

    ws[c] = ws_prev[c];
    double factor = dT / (phi[c] * mesh_->cell_volume(c));
    for (int n = 0; n < nfaces; n++) {
      int f = faces[n];
      ws[c] -= factor * flux_d[f] * dirs[n]; 
    }
  }
}


/* ******************************************************************
* Estimate du/dt from the pressure equations, du/dt = g - A*u.
****************************************************************** */
double Richards_PK::ComputeUDot(double T, const Epetra_Vector& u, Epetra_Vector& udot)
{
  ComputePreconditionerMFD(u, matrix, T, 0.0, false);  // Calculate only stiffness matrix.
  double norm_udot = matrix->ComputeNegativeResidual(u, udot);

  Epetra_Vector* udot_faces = FS->CreateFaceView(udot);
  udot_faces->PutScalar(0.0);

  return norm_udot;
}


/* ******************************************************************
* Gathers together routines to compute MFD matrices Axx(u) and 
* preconditioner Sff(u) using internal time step dT.                             
****************************************************************** */
void Richards_PK::ComputePreconditionerMFD(
    const Epetra_Vector& u, Matrix_MFD* matrix,
    double Tp, double dTp, bool flag_update_ML)
{
  // setup absolute and compute relative permeabilities
  Epetra_Vector* u_cells = FS->CreateCellView(u);
  Epetra_Vector* u_faces = FS->CreateFaceView(u);

  if (!is_matrix_symmetric) {
    CalculateRelativePermeabilityFace(*u_cells);
    Krel_cells->PutScalar(1.0);
  } else {
    CalculateRelativePermeabilityCell(*u_cells);
    Krel_faces->PutScalar(1.0);
  }

  // update boundary conditions
  bc_pressure->Compute(Tp);
  bc_flux->Compute(Tp);
  bc_head->Compute(Tp);
  bc_seepage->Compute(Tp);
  UpdateBoundaryConditions(
      bc_pressure, bc_head, bc_flux, bc_seepage,
      *u_faces, atm_pressure,
      bc_markers, bc_values);

  // setup a new algebraic problem
  matrix->CreateMFDstiffnessMatrices(*Krel_cells, *Krel_faces);
  matrix->CreateMFDrhsVectors();
  AddGravityFluxes_MFD(K, *Krel_cells, *Krel_faces, matrix);
  if (flag_update_ML) AddTimeDerivative_MFD(*u_cells, dTp, matrix);
  matrix->ApplyBoundaryConditions(bc_markers, bc_values);
  matrix->AssembleGlobalMatrices();
  if (flag_update_ML) {
    matrix->ComputeSchurComplement(bc_markers, bc_values);
    matrix->UpdateML_Preconditioner();
  }

  // DEBUG
  // Matrix_Audit audit(mesh_, matrix);
  // audit.InitAudit();
  // audit.CheckSpectralBounds();
}


/* ******************************************************************
* BDF methods need a good initial guess.
* This method gives a less smoother solution than in Flow 1.0.
* WARNING: Each owned face must have at least one owned cell. 
* Probability that this assumption is violated is close to zero. 
* Even when it happens, the code will not crash.
****************************************************************** */
void Richards_PK::DeriveFaceValuesFromCellValues(const Epetra_Vector& ucells, Epetra_Vector& ufaces)
{
  AmanziMesh::Entity_ID_List cells;

  for (int f = 0; f < nfaces_owned; f++) {
    cells.clear();
    mesh_->face_get_cells(f, AmanziMesh::OWNED, &cells);
    int ncells = cells.size();

    if (ncells > 0) {
      double face_value = 0.0;
      for (int n = 0; n < ncells; n++) face_value += ucells[cells[n]];
      ufaces[f] = face_value / ncells;
    } else {
      ufaces[f] = atm_pressure;
    }
  }
}


/* ******************************************************************
* Temporary convertion from double to tensor.                                               
****************************************************************** */
void Richards_PK::SetAbsolutePermeabilityTensor(std::vector<WhetStone::Tensor>& K)
{
  const Epetra_Vector& vertical_permeability = FS->ref_vertical_permeability();
  const Epetra_Vector& horizontal_permeability = FS->ref_horizontal_permeability();

  for (int c = 0; c < K.size(); c++) {
    if (vertical_permeability[c] == horizontal_permeability[c]) {
      K[c].init(dim, 1);
      K[c](0, 0) = vertical_permeability[c];
    } else {
      K[c].init(dim, 2);
      for (int i = 0; i < dim-1; i++) K[c](i, i) = horizontal_permeability[c];
      K[c](dim-1, dim-1) = vertical_permeability[c];
    }
  }
}


/* ******************************************************************
* Adds time derivative to the cell-based part of MFD algebraic system.                                               
****************************************************************** */
void Richards_PK::AddTimeDerivative_MFD(
    Epetra_Vector& pressure_cells, double dT_prec, Matrix_MFD* matrix)
{
  Epetra_Vector dSdP(mesh_->cell_map(false));
  DerivedSdP(pressure_cells, dSdP);

  const Epetra_Vector& phi = FS->ref_porosity();
  std::vector<double>& Acc_cells = matrix->Acc_cells();
  std::vector<double>& Fc_cells = matrix->Fc_cells();

  int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);

  for (int c = 0; c < ncells; c++) {
    double volume = mesh_->cell_volume(c);
    double factor = rho * phi[c] * dSdP[c] * volume / dT_prec;
    Acc_cells[c] += factor;
    Fc_cells[c] += factor * pressure_cells[c];
  }
}


/* ******************************************************************
* Initialize saturated pressure solutions using boundary conditions 
* at time Tp.
* WARNING: data in vectors Krel and rhs are destroyed.
****************************************************************** */
void Richards_PK::InitializePressureHydrostatic(const double Tp)
{
  AztecOO* solver_tmp = new AztecOO;

  // update boundary conditions
  bc_pressure->Compute(Tp);
  bc_head->Compute(Tp);
  bc_flux->Compute(Tp);
  bc_seepage->Compute(Tp);
  UpdateBoundaryConditions(
      bc_pressure, bc_head, bc_flux, bc_seepage,
      *solution_faces, atm_pressure,
      bc_markers, bc_values);

  // set fully saturated media
  Krel_cells->PutScalar(1.0);
  Krel_faces->PutScalar(1.0);

  // calculate and assemble elemental stifness matrices
  matrix->CreateMFDstiffnessMatrices(*Krel_cells, *Krel_faces);
  matrix->CreateMFDrhsVectors();
  AddGravityFluxes_MFD(K, *Krel_cells, *Krel_faces, matrix);
  matrix->ApplyBoundaryConditions(bc_markers, bc_values);
  matrix->AssembleGlobalMatrices();

  preconditioner->CreateMFDstiffnessMatrices(*Krel_cells, *Krel_faces);
  preconditioner->CreateMFDrhsVectors();
  AddGravityFluxes_MFD(K, *Krel_cells, *Krel_faces, preconditioner);
  preconditioner->ApplyBoundaryConditions(bc_markers, bc_values);
  preconditioner->AssembleGlobalMatrices();
  preconditioner->ComputeSchurComplement(bc_markers, bc_values);
  preconditioner->UpdateML_Preconditioner();

  // solve symmetric problem
  solver_tmp->SetUserOperator(matrix);
  solver_tmp->SetPrecOperator(preconditioner);
  solver_tmp->SetAztecOption(AZ_solver, AZ_cg);
  solver_tmp->SetAztecOption(AZ_output, AZ_none);
  solver_tmp->SetAztecOption(AZ_conv, AZ_rhs);

  Epetra_Vector b(*(matrix->rhs()));
  solver_tmp->SetRHS(&b);

  solver_tmp->SetLHS(&*solution);
  solver_tmp->Iterate(max_itrs, convergence_tol);

  if (MyPID == 0 && verbosity >= FLOW_VERBOSITY_HIGH) {
    int num_itrs = solver_tmp->NumIters();
    double linear_residual = solver_tmp->TrueResidual();
    std::printf("Richards PK: initial pressure solver(%8.3e, %4d)\n", linear_residual, num_itrs);
  }

  delete solver_tmp;
}


/* ******************************************************************
* A wrapper for a similar matrix call.
****************************************************************** */
void Richards_PK::DeriveDarcyVelocity(const Epetra_Vector& flux, Epetra_MultiVector& velocity)
{
  matrix->DeriveDarcyVelocity(flux, *face_importer_, velocity);
}


}  // namespace AmanziFlow
}  // namespace Amanzi

