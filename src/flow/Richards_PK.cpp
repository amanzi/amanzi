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
Richards_PK::Richards_PK(Teuchos::ParameterList& flow_list, Teuchos::RCP<Flow_State> FS_MPC)
{
  Flow_PK::Init(FS_MPC);

  FS = FS_MPC;
  if (flow_list.isSublist("Richards Problem")) {
    rp_list = flow_list.sublist("Richards Problem");
  } else {
    Errors::Message msg("Flow does not have <Richards Problem> sublist.");
    Exceptions::amanzi_throw(msg);
  }

  mesh_ = FS->mesh();
  dim = mesh_->space_dimension();

  // Create the combined cell/face DoF map.
  super_map_ = createSuperMap();

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

  ti_method_sss = FLOW_TIME_INTEGRATION_BDF1;  // time integration (TI) parameters
  ti_method_trs = FLOW_TIME_INTEGRATION_BDF2;
  num_itrs_trs = 0;

  absolute_tol_sss = absolute_tol_trs = 1.0;
  relative_tol_sss = relative_tol_trs = 1e-5;
  initialize_with_darcy = 0;

  mfd3d_method = FLOW_MFD3D_HEXAHEDRA_MONOTONE;  // will be changed (lipnikov@lanl.gov)
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
  solution_cells = Teuchos::rcp(FS->createCellView(*solution));
  solution_faces = Teuchos::rcp(FS->createFaceView(*solution));
  rhs = Teuchos::rcp(new Epetra_Vector(*super_map_));
  rhs = matrix->get_rhs();  // import rhs from the matrix

  // Get solver parameters from the flow parameter list.
  ProcessParameterList();

  // Process boundary data (state may be incomplete at this moment)
  int nfaces = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::USED);
  bc_markers.resize(nfaces, FLOW_BC_FACE_NULL);
  bc_values.resize(nfaces, 0.0);

  double T_physical = FS->get_time();
  T_internal = (standalone_mode) ? T_internal : T_physical;

  double time = T_internal;
  bc_pressure->Compute(time);
  bc_flux->Compute(time);
  bc_head->Compute(time);
  bc_seepage->Compute(time);
  UpdateBoundaryConditions(
      bc_pressure, bc_head, bc_flux, bc_seepage,
      *solution_cells, atm_pressure,
      bc_markers, bc_values);

  // Process other fundamental structures
  K.resize(ncells_owned);
  is_matrix_symmetric = (Krel_method == FLOW_RELATIVE_PERM_CENTERED);
  matrix->setSymmetryProperty(is_matrix_symmetric);
  matrix->symbolicAssembleGlobalMatrices(*super_map_);

  // Allocate data for relative permeability
  const Epetra_Map& cmap = mesh_->cell_map(true);
  const Epetra_Map& fmap = mesh_->face_map(true);

  Krel_cells = Teuchos::rcp(new Epetra_Vector(cmap));
  Krel_faces = Teuchos::rcp(new Epetra_Vector(fmap));

  Krel_cells->PutScalar(1.0);  // we start with fully saturated media
  Krel_faces->PutScalar(1.0);

  if (Krel_method == FLOW_RELATIVE_PERM_UPWIND_GRAVITY) {
    // Kgravity_unit.resize(ncells_wghost);  Resize does not work properly.
    CalculateKVectorUnit(gravity_, Kgravity_unit);
  }

  if (MyPID == 0 && verbosity >= FLOW_VERBOSITY_MEDIUM) {
    if (mfd3d_method == FLOW_MFD3D_HEXAHEDRA_MONOTONE) {
      std::printf("Richards Flow: discretization method is for orthogonal hexes.\n");
    } else if (mfd3d_method == FLOW_MFD3D_POLYHEDRA) {
      std::printf("Richards Flow: discretization method is for generic polyhedra.\n");
    }
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
     std::printf("Richards Flow: initializing steady-state at T(sec)=%9.4e dT(sec)=%9.4e \n", T0, dT0);
     if (initialize_with_darcy) {
       std::printf("Richards Flow: initializing with a clipped Darcy pressure\n");
       std::printf("Richards Flow: clipping saturation value =%5.2g\n", clip_saturation);
     }
  }

  Teuchos::ParameterList ML_list = rp_list.sublist("Diffusion Preconditioner").sublist("ML Parameters");

  preconditioner->setSymmetryProperty(is_matrix_symmetric);
  preconditioner->symbolicAssembleGlobalMatrices(*super_map_);
  preconditioner->init_ML_preconditioner(ML_list);

  if (ti_method_sss == FLOW_TIME_INTEGRATION_BDF2) {
    Teuchos::ParameterList solver_list = rp_list.sublist("steady state time integrator").sublist("nonlinear solver BDF2");
    if (solver_list.isSublist("VerboseObject"))
       solver_list.sublist("VerboseObject") = rp_list.sublist("VerboseObject");

    Teuchos::RCP<Teuchos::ParameterList> bdf2_list(new Teuchos::ParameterList(solver_list));
    if (bdf2_dae == NULL) bdf2_dae = new BDF2::Dae(*this, *super_map_);
    bdf2_dae->setParameterList(bdf2_list);

  } else if (ti_method_sss == FLOW_TIME_INTEGRATION_BDF1) {
    Teuchos::ParameterList solver_list = rp_list.sublist("steady state time integrator").sublist("nonlinear solver BDF1");
    if (solver_list.isSublist("VerboseObject"))
        solver_list.sublist("VerboseObject") = rp_list.sublist("VerboseObject");

    Teuchos::RCP<Teuchos::ParameterList> bdf1_list(new Teuchos::ParameterList(solver_list));
    if (bdf1_dae == NULL) bdf1_dae = new BDF1Dae(*this, *super_map_);
    bdf1_dae->setParameterList(bdf1_list);

  } else {
    solver = new AztecOO;
    solver->SetUserOperator(matrix);
    solver->SetPrecOperator(preconditioner);
    solver->SetAztecOption(AZ_solver, AZ_cg);  // symmetry is required
  }

  // (re)initialize pressure and saturation
  Epetra_Vector& pressure = FS->ref_pressure();
  Epetra_Vector& water_saturation = FS->ref_water_saturation();

  *solution_cells = pressure;

  if (initialize_with_darcy) {
    double pmin = atm_pressure;
    InitializePressureHydrostatic(T0, *solution);
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

  flow_status_ = FLOW_STATUS_STEADY_STATE_INIT;

  // DEBUG
  // AdvanceToSteadyState();
  // CommitStateForTransport(FS); CommitState(FS); writeGMVfile(FS); exit(0);
}


/* ******************************************************************
* Initialization analyzes status of matrix/preconditioner pair. 
* BDF2 and BDF1 will eventually merge but are separated strictly 
* (no code optimization) for the moment.     
****************************************************************** */
void Richards_PK::InitTransient(double T0, double dT0)
{
  if (MyPID == 0 && verbosity >= FLOW_VERBOSITY_MEDIUM) {
     std::printf("Initializing Transient Flow: T(sec)=%9.4e dT(sec)=%9.4e\n", T0, dT0);
  }
  Teuchos::ParameterList ML_list = rp_list.sublist("Diffusion Preconditioner").sublist("ML Parameters");

  preconditioner->setSymmetryProperty(is_matrix_symmetric);
  preconditioner->symbolicAssembleGlobalMatrices(*super_map_);
  preconditioner->init_ML_preconditioner(ML_list);

  if (ti_method_trs == FLOW_TIME_INTEGRATION_BDF2) {
    if (bdf2_dae != NULL) delete bdf2_dae;  // The only way to reset BDF2 is to delete it.

    Teuchos::ParameterList solver_list = rp_list.sublist("transient time integrator").sublist("nonlinear solver BDF2");
    if (solver_list.isSublist("VerboseObject"))
        solver_list.sublist("VerboseObject") = rp_list.sublist("VerboseObject");

    Teuchos::RCP<Teuchos::ParameterList> bdf2_list(new Teuchos::ParameterList(solver_list));
    bdf2_dae = new BDF2::Dae(*this, *super_map_);
    bdf2_dae->setParameterList(bdf2_list);

  } else if (ti_method_trs == FLOW_TIME_INTEGRATION_BDF1) {
    if (bdf1_dae != NULL) delete bdf1_dae;  // the only way to reset BDF1 is to delete it

    Teuchos::ParameterList solver_list = rp_list.sublist("transient time integrator").sublist("nonlinear solver BDF1");
    if (solver_list.isSublist("VerboseObject"))
        solver_list.sublist("VerboseObject") = rp_list.sublist("VerboseObject");

    Teuchos::RCP<Teuchos::ParameterList> bdf1_list(new Teuchos::ParameterList(solver_list));
    bdf1_dae = new BDF1Dae(*this, *super_map_);
    bdf1_dae->setParameterList(bdf1_list);

  } else if (solver == NULL) {
    solver = new AztecOO;
    solver->SetUserOperator(matrix);
    solver->SetPrecOperator(preconditioner);
    solver->SetAztecOption(AZ_solver, AZ_cg);  // symmetry is required
  }

  // initialize pressure and saturation
  Epetra_Vector& pressure = FS->ref_pressure();
  Epetra_Vector& water_saturation = FS->ref_water_saturation();
  *solution_cells = pressure;

  DeriveFaceValuesFromCellValues(*solution_cells, *solution_faces);
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

  double dTnext;
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
  dT = dTnext;
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
  Epetra_Vector& pressure = FS_MPC->ref_pressure();
  pressure = *solution_cells;

  Epetra_Vector& water_saturation = FS_MPC->ref_water_saturation();
  FS_MPC->ref_prev_water_saturation() = water_saturation;
  DeriveSaturationFromPressure(pressure, water_saturation);

  Epetra_Vector& flux = FS_MPC->ref_darcy_flux();
  matrix->createMFDstiffnessMatrices(mfd3d_method, K, *Krel_faces);  // Should be improved. (lipnikov@lanl.gov)
  matrix->deriveDarcyMassFlux(*solution, *face_importer_, flux);
  addGravityFluxes_DarcyFlux(K, *Krel_faces, flux);
  for (int c = 0; c < nfaces_owned; c++) flux[c] /= rho;

  // DEBUG
  // if (MyPID == 0 && verbosity >= FLOW_VERBOSITY_HIGH) {
  //   std::printf("Commited flow state for transport: see flow.gmv\n");
  // }
  // writeGMVfile(FS_MPC);
}


/* ******************************************************************
* Transfer internal data to flow state FS_MPC. MPC may request
* to populate the original state FS. 
****************************************************************** */
void Richards_PK::CommitState(Teuchos::RCP<Flow_State> FS_MPC)
{
  // Epetra_Vector& pressure = FS_MPC->ref_pressure();
  // pressure = *solution_cells;

  Epetra_Vector& flux = FS_MPC->ref_darcy_flux();
  Epetra_MultiVector& velocity = FS_MPC->ref_darcy_velocity();
  DeriveDarcyVelocity(flux, velocity);
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
* Estimate du/dt from the pressure equations, -du/dt = g - A*u.
****************************************************************** */
double Richards_PK::ComputeUDot(double T, const Epetra_Vector& u, Epetra_Vector& udot)
{
  double norm_udot;
  ComputePreconditionerMFD(u, matrix, mfd3d_method, T, 0.0, false);  // Calculate only stiffness matrix.
  norm_udot = matrix->computeNegativeResidual(u, udot);

  Epetra_Vector* udot_faces = FS->createFaceView(udot);
  udot_faces->PutScalar(0.0);

  return norm_udot;
}


/* ******************************************************************
* Gathers together routines to compute MFD matrices Axx(u) and 
* preconditioner Sff(u) using internal time step dT.                             
****************************************************************** */
void Richards_PK::ComputePreconditionerMFD(
    const Epetra_Vector& u, Matrix_MFD* matrix, int disc_method, 
    double Tp, double dTp, bool flag_update_ML)
{
  // setup absolute and compute relative permeabilities
  Epetra_Vector* u_cells = FS->createCellView(u);

  SetAbsolutePermeabilityTensor(K);
  if (!is_matrix_symmetric) {  // Define K and Krel_faces
    CalculateRelativePermeabilityFace(*u_cells);
    for (int c = 0; c < K.size(); c++) K[c] *= rho / mu;
  } else {  // Define K and Krel_cells, Krel_faces is always one
    CalculateRelativePermeabilityCell(*u_cells);
    for (int c = 0; c < K.size(); c++) K[c] *= (*Krel_cells)[c] * rho / mu;
  }

  // update boundary conditions
  bc_pressure->Compute(Tp);
  bc_flux->Compute(Tp);
  bc_head->Compute(Tp);
  bc_seepage->Compute(Tp);
  UpdateBoundaryConditions(
      bc_pressure, bc_head, bc_flux, bc_seepage,
      *u_cells, atm_pressure,
      bc_markers, bc_values);

  // setup a new algebraic problem
  matrix->createMFDstiffnessMatrices(disc_method, K, *Krel_faces);
  matrix->createMFDrhsVectors();
  addGravityFluxes_MFD(K, *Krel_faces, matrix);
  if (flag_update_ML) AddTimeDerivative_MFD(*u_cells, dTp, matrix);
  matrix->applyBoundaryConditions(bc_markers, bc_values);
  matrix->assembleGlobalMatrices();
  if (flag_update_ML) {
    matrix->computeSchurComplement(bc_markers, bc_values);
    matrix->update_ML_preconditioner();
  }

  // DEBUG
  //Matrix_Audit audit(mesh_, matrix);
  //audit.InitAudit();
  //audit.CheckSpectralBounds();
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
  std::vector<double>& Acc_cells = matrix->get_Acc_cells();
  std::vector<double>& Fc_cells = matrix->get_Fc_cells();

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
void Richards_PK::InitializePressureHydrostatic(const double Tp, Epetra_Vector& u)
{
  AztecOO* solver_tmp = new AztecOO;

  // update boundary conditions
  bc_pressure->Compute(Tp);
  bc_head->Compute(Tp);
  bc_flux->Compute(Tp);
  bc_seepage->Compute(Tp);
  UpdateBoundaryConditions(
      bc_pressure, bc_head, bc_flux, bc_seepage,
      *solution_cells, atm_pressure,
      bc_markers, bc_values);

  // work-around limited support for tensors
  SetAbsolutePermeabilityTensor(K);
  for (int c = 0; c < K.size(); c++) K[c] *= rho / mu;
  Krel_faces->PutScalar(1.0);

  // calculate and assemble elemental stifness matrices
  matrix->createMFDstiffnessMatrices(mfd3d_method, K, *Krel_faces);
  matrix->createMFDrhsVectors();
  addGravityFluxes_MFD(K, *Krel_faces, matrix);
  matrix->applyBoundaryConditions(bc_markers, bc_values);
  matrix->assembleGlobalMatrices();

  int disc_method = AmanziFlow::FLOW_MFD3D_TWO_POINT_FLUX;
  preconditioner->createMFDstiffnessMatrices(disc_method, K, *Krel_faces);
  preconditioner->createMFDrhsVectors();
  addGravityFluxes_MFD(K, *Krel_faces, preconditioner);
  preconditioner->applyBoundaryConditions(bc_markers, bc_values);
  preconditioner->assembleGlobalMatrices();
  preconditioner->computeSchurComplement(bc_markers, bc_values);
  preconditioner->update_ML_preconditioner();

  // solve symmetric problem
  solver_tmp->SetUserOperator(preconditioner);
  solver_tmp->SetPrecOperator(preconditioner);
  solver_tmp->SetAztecOption(AZ_solver, AZ_cg);
  solver_tmp->SetAztecOption(AZ_output, AZ_none);
  solver_tmp->SetAztecOption(AZ_conv, AZ_rhs);

  rhs = matrix->get_rhs();
  solver_tmp->SetRHS(&*rhs);

  u.PutScalar(0.0);
  solver_tmp->SetLHS(&u);
  solver_tmp->Iterate(max_itrs, convergence_tol);

  if (MyPID == 0 && verbosity >= FLOW_VERBOSITY_HIGH) {
    int num_itrs = solver_tmp->NumIters();
    double linear_residual = solver_tmp->TrueResidual();
    std::printf("Initial pressure: linear solver(%8.3e, %4d)\n", linear_residual, num_itrs);
  }

  delete solver_tmp;
}


/* ******************************************************************
* A wrapper for a similar matrix call.
****************************************************************** */
void Richards_PK::DeriveDarcyVelocity(const Epetra_Vector& flux, Epetra_MultiVector& velocity)
{
  matrix->deriveDarcyVelocity(flux, *face_importer_, velocity);
}


}  // namespace AmanziFlow
}  // namespace Amanzi

