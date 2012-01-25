/*
This is the flow component of the Amanzi code. 
License: BSD
Authors: Neil Carlson (nnc@lanl.gov), 
         Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include "Epetra_IntVector.h"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

#include "dbc.hh"
#include "errors.hh"
#include "exceptions.hh"

#include "Mesh.hh"
#include "Point.hh"
#include "gmv_mesh.hh"

#include "Flow_BC_Factory.hpp"
#include "boundary-function.hh"
#include "Richards_PK.hpp"

namespace Amanzi {
namespace AmanziFlow {

/* ******************************************************************
* We set up only default values and call Init() routine to complete
* each variable initialization
****************************************************************** */
Richards_PK::Richards_PK(Teuchos::ParameterList& rp_list_, Teuchos::RCP<Flow_State> FS_MPC)
{
  Flow_PK::Init(FS_MPC);

  FS = FS_MPC;
  rp_list = rp_list_;

  mesh_ = FS->get_mesh();
  dim = mesh_->space_dimension();

  // Create the combined cell/face DoF map.
  super_map_ = createSuperMap();

  // Other fundamental physical quantaties
  rho = *(FS->get_fluid_density());
  mu = *(FS->get_fluid_viscosity()); 
  gravity.init(dim);
  for (int k=0; k<dim; k++) gravity[k] = (*(FS->get_gravity()))[k];

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

  method_sss = FLOW_STEADY_STATE_BACKWARD_EULER;
  method_trs = FLOW_STEADY_STATE_BDF2;
  num_itrs_trs = 0;

  absolute_tol_sss = absolute_tol_trs = 1.0; 
  relative_tol_sss = relative_tol_trs = 1e-5;

  mfd3d_method = FLOW_MFD3D_HEXAHEDRA_MONOTONE;  // will be changed (lipnikov@lanl.gov)
  flag_upwind = true;

  verbosity = FLOW_VERBOSITY_HIGH;
  internal_tests = 0;
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
}


/* ******************************************************************
* Extract information from Richards Problem parameter list.
****************************************************************** */
void Richards_PK::Init(Matrix_MFD* matrix_, Matrix_MFD* preconditioner_)
{
  if (matrix_ == NULL) matrix = new Matrix_MFD(FS, *super_map_);
  else matrix = matrix_;

  if (preconditioner_ == NULL) preconditioner = matrix;
  else preconditioner = preconditioner_;

  // Create the solution (pressure) vector.
  solution = Teuchos::rcp(new Epetra_Vector(*super_map_));
  solution_cells = Teuchos::rcp(FS->createCellView(*solution));
  solution_faces = Teuchos::rcp(FS->createFaceView(*solution));
  rhs = Teuchos::rcp(new Epetra_Vector(*super_map_));
  rhs = matrix->get_rhs();  // import rhs from the matrix 

  // Get some solver parameters from the flow parameter list.
  processParameterList();

  // Process boundary data
  int nfaces = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::USED);
  bc_markers.resize(nfaces, FLOW_BC_FACE_NULL);
  bc_values.resize(nfaces, 0.0);

  T_physical = FS->get_time();
  double time = (standalone_mode) ? T_internal : T_physical;

  bc_pressure->Compute(time);
  bc_flux->Compute(time);
  bc_head->Compute(time);
  updateBoundaryConditions(bc_pressure, bc_head, bc_flux, bc_markers, bc_values);

  // Process other fundamental structures
  K.resize(ncells_owned);
  matrix->setSymmetryProperty(!flag_upwind);
  matrix->symbolicAssembleGlobalMatrices(*super_map_);

  // Create the BDF2 time integrator
  Teuchos::ParameterList solver_list = rp_list.sublist("Steady state solution").sublist("Nonlinear solvers");
  Teuchos::RCP<Teuchos::ParameterList> bdf2_list(new Teuchos::ParameterList(solver_list));
  bdf2_dae = new BDF2::Dae(*this, *super_map_);
  bdf2_dae->setParameterList(bdf2_list);

  // Preconditioner
  Teuchos::ParameterList ML_list = rp_list.sublist("Diffusion Preconditioner").sublist("ML Parameters");

  if (method_sss == FLOW_STEADY_STATE_BDF2) {
    preconditioner = new Matrix_MFD(FS, *super_map_);
    preconditioner->setSymmetryProperty(!flag_upwind);
    preconditioner->symbolicAssembleGlobalMatrices(*super_map_);
    preconditioner->init_ML_preconditioner(ML_list); 
  } else {
    preconditioner->init_ML_preconditioner(ML_list);
    solver = new AztecOO;
    solver->SetUserOperator(matrix);
    solver->SetPrecOperator(preconditioner);
    solver->SetAztecOption(AZ_solver, AZ_cg);  // symmetry is required
  }

  // Allocate data for relative permeability
  const Epetra_Map& cmap = mesh_->cell_map(true);
  const Epetra_Map& fmap = mesh_->face_map(true);
 
  Krel_cells = Teuchos::rcp(new Epetra_Vector(cmap));
  Krel_faces = Teuchos::rcp(new Epetra_Vector(fmap));

  Krel_cells->PutScalar(1.0);  // we start with fully saturated media
  Krel_faces->PutScalar(1.0);
}


/* ******************************************************************
*  Wrapper for advance to steady-state routines.                                                    
****************************************************************** */
int Richards_PK::advance_to_steady_state()
{
  for (int c=0; c<ncells_owned; c++) (*solution_cells)[c] = FS->ref_pressure()[c];
  deriveFaceValuesFromCellValues(*solution_cells, *solution_faces);
  applyBoundaryConditions(bc_markers, bc_values, *solution_faces);

  int ierr = 0;
  if (method_sss == FLOW_STEADY_STATE_PICARD) {
    ierr = advanceSteadyState_Picard();
  } else if (method_sss == FLOW_STEADY_STATE_BACKWARD_EULER) {
    ierr = advanceSteadyState_BackwardEuler();
  } else if (method_sss == FLOW_STEADY_STATE_BDF2) {
    ierr = advanceSteadyState_BDF2();
  }

  //create flow state to be returned
  FS_next->ref_pressure() = *solution_cells;

  // calculate darcy mass flux
  Epetra_Vector& darcy_mass_flux = FS_next->ref_darcy_mass_flux();
  matrix->createMFDstiffnessMatrices(mfd3d_method, K, *Krel_faces);  // Should be improved. (lipnikov@lanl.gov)
  matrix->deriveDarcyFlux(*solution, *face_importer_, darcy_mass_flux);
  addGravityFluxes_DarcyFlux(K, *Krel_faces, darcy_mass_flux);

  status = FLOW_NEXT_STATE_COMPLETE;
  return ierr;
}


/* ******************************************************************* 
* Performs one time step of size dT (under *development*). 
******************************************************************* */
int Richards_PK::advance(double dT) 
{
  T_physical = FS->get_time();
  double time = (standalone_mode) ? T_internal : T_physical;
  double dTnext;

  if (num_itrs_trs == 0) {  // initialization
    Epetra_Vector udot(*super_map_);
    computeUDot(time, *solution, udot);
    bdf2_dae->set_initial_state(time, *solution, udot);

    int ierr;
    update_precon(time, *solution, dT, ierr);
  }

  bdf2_dae->bdf2_step(dT, 0.0, *solution, dTnext);
  bdf2_dae->commit_solution(dT, *solution);
  bdf2_dae->write_bdf2_stepping_statistics();

  // calculate darcy mass flux
  Epetra_Vector& darcy_mass_flux = FS_next->ref_darcy_mass_flux();
  matrix->createMFDstiffnessMatrices(mfd3d_method, K, *Krel_faces);  // Should be improved. (lipnikov@lanl.gov)
  matrix->deriveDarcyFlux(*solution, *face_importer_, darcy_mass_flux);
  addGravityFluxes_DarcyFlux(K, *Krel_faces, darcy_mass_flux);

  num_itrs_trs++;
  return 0;
}


/* ******************************************************************* 
* Performs one time step of size dT. 
******************************************************************* */
int Richards_PK::advanceSteadyState_BDF2() 
{
  T_internal = T0_sss;
  dT = dT0_sss;

  int itrs = 0;
  while (itrs < max_itrs_sss && T_internal < T1_sss) {
    if (itrs == 0) {  // initialization of BDF2
      Epetra_Vector udot(*super_map_);
      computeUDot(T0_sss, *solution, udot);
      bdf2_dae->set_initial_state(T0_sss, *solution, udot);

      int ierr;
      update_precon(T0_sss, *solution, dT0_sss, ierr);
    }

    double dTnext;
    bdf2_dae->bdf2_step(dT, 0.0, *solution, dTnext);
    bdf2_dae->commit_solution(dT, *solution);
    bdf2_dae->write_bdf2_stepping_statistics();

    T_internal = bdf2_dae->most_recent_time();
    dT = dTnext;
    itrs++;
  }
  return 0;
}


/* ******************************************************************
* Calculates steady-state solution assuming that abosolute and
* relative permeabilities do not depend explicitly on time.
* This is the experimental method.                                                 
****************************************************************** */
int Richards_PK::advanceSteadyState_Picard()
{
  Epetra_Vector  solution_old(*solution);
  Epetra_Vector& solution_new = *solution;
  Epetra_Vector  residual(*solution);

  if (flag_upwind) solver->SetAztecOption(AZ_solver, AZ_bicgstab);  // symmetry is NOT required
  solver->SetAztecOption(AZ_output, AZ_none);

  int itrs = 0;
  double L2norm, L2error = 1.0;
  double relaxation = 1e-1;

  while (L2error > convergence_tol_sss && itrs < max_itrs_sss) {
    calculateRelativePermeability(*solution_cells);
    setAbsolutePermeabilityTensor(K);
    if (flag_upwind) {  // Define K and Krel_faces
      calculateRelativePermeabilityUpwindGravity(*solution_cells);
      for (int c=0; c<K.size(); c++) K[c] *= rho / mu;
    } else {  // Define K and Krel_cells, Krel_faces is always one
      for (int c=0; c<K.size(); c++) K[c] *= (*Krel_cells)[c] * rho / mu;  
    }

    // create algebraic problem (matrix = preconditioner)
    matrix->createMFDstiffnessMatrices(mfd3d_method, K, *Krel_faces);
    matrix->createMFDrhsVectors();
    addGravityFluxes_MFD(K, *Krel_faces, matrix);
    matrix->applyBoundaryConditions(bc_markers, bc_values);
    matrix->assembleGlobalMatrices();
    matrix->computeSchurComplement(bc_markers, bc_values);
    matrix->update_ML_preconditioner();
    rhs = matrix->get_rhs();  // export rhs from matrix class

    // check convergence of non-linear residual
    L2error = matrix->computeResidual(solution_new, residual);
    residual.Norm2(&L2error);
    rhs->Norm2(&L2norm);
    L2error /= L2norm;

    // call AztecOO solver
    Epetra_Vector b(*rhs);
    solver->SetRHS(&b);  // AztecOO modifies the right-hand-side.
    solver->SetLHS(&*solution);  // initial solution guess 

    solver->Iterate(max_itrs, convergence_tol);
    int num_itrs = solver->NumIters();
    double linear_residual = solver->TrueResidual();

    if (MyPID == 0 && verbosity >= FLOW_VERBOSITY_HIGH) {
      std::printf("Picard:%4d   Pressure(res=%9.4e, rhs=%9.4e, relax=%8.3e)  solver(%8.3e, %4d)\n", 
          itrs, L2error, L2norm, relaxation, linear_residual, num_itrs);
    }

    for (int c=0; c<ncells_owned; c++) {
      solution_new[c] = (1.0 - relaxation) * solution_old[c] + relaxation * solution_new[c];
      solution_old[c] = solution_new[c];
    }

    T_internal += dT;
    itrs++;
  }

  return 0;
}


/* ******************************************************************
* BDF methods need a good initial guess.
****************************************************************** */
void Richards_PK::deriveFaceValuesFromCellValues(const Epetra_Vector& ucells, Epetra_Vector& ufaces)
{
  AmanziMesh::Entity_ID_List cells; 

  for (int f=0; f<nfaces_owned; f++) {
    cells.clear();
    mesh_->face_get_cells(f, AmanziMesh::USED, &cells);
    int ncells = cells.size();
         
    double face_value = 0.0;
    for (int n=0; n<ncells; n++) face_value += ucells[cells[n]];
	
    ufaces[f] = face_value / ncells;
  }
}


/* ******************************************************************
* Estimate du/dt from the pressure equations, du/dt = g - A*u.
****************************************************************** */
double Richards_PK::computeUDot(double T, const Epetra_Vector& u, Epetra_Vector& udot)
{
  double norm_udot;
  computePreconditionerMFD(u, matrix, T, 0.0, false);  // Calculate only stiffness matrix.
  norm_udot = matrix->computeResidual(u, udot);
  udot.Update(-1.0, udot, 0.0);

  Epetra_Vector* udot_faces = FS->createFaceView(udot);
  udot_faces->PutScalar(0.0);

  return norm_udot;
}


/* ******************************************************************
* Gathers together routines to compute MFD matrices Axx(u) and 
* preconditioner Sff(u) using internal time step dT.                             
****************************************************************** */
void Richards_PK::computePreconditionerMFD(
    const Epetra_Vector& u, Matrix_MFD* matrix, double Tp, double dTp, bool flag_update_ML)
{
  // setup absolute and compute relative permeabilities
  Epetra_Vector* u_cells = FS->createCellView(u);
  calculateRelativePermeability(*u_cells);
  setAbsolutePermeabilityTensor(K);

  if (flag_upwind) {  // Define K and Krel_faces
    calculateRelativePermeabilityUpwindGravity(*u_cells);
    for (int c=0; c<K.size(); c++) K[c] *= rho / mu;
  } else {  // Define K and Krel_cells, Krel_faces is always one
    for (int c=0; c<K.size(); c++) K[c] *= (*Krel_cells)[c] * rho / mu;  
  }

  // update boundary conditions
  bc_pressure->Compute(Tp);
  bc_flux->Compute(Tp);
  bc_head->Compute(Tp);
  updateBoundaryConditions(bc_pressure, bc_head, bc_flux, bc_markers, bc_values);

  // setup new algebraic problems
  matrix->createMFDstiffnessMatrices(mfd3d_method, K, *Krel_faces);
  matrix->createMFDrhsVectors();
  addGravityFluxes_MFD(K, *Krel_faces, matrix);
  if (flag_update_ML) addTimeDerivative_MFD(*u_cells, dTp, matrix);
  matrix->applyBoundaryConditions(bc_markers, bc_values);
  matrix->assembleGlobalMatrices();
  if (flag_update_ML) {
    matrix->computeSchurComplement(bc_markers, bc_values);
    matrix->update_ML_preconditioner();  
  }
}


/* ******************************************************************
* Temporary convertion from double to tensor.                                               
****************************************************************** */
void Richards_PK::setAbsolutePermeabilityTensor(std::vector<WhetStone::Tensor>& K)
{
  const Epetra_Vector& vertical_permeability = FS->ref_vertical_permeability();
  const Epetra_Vector& horizontal_permeability = FS->ref_horizontal_permeability();

  for (int c=0; c<K.size(); c++) {
    if (vertical_permeability[c] == horizontal_permeability[c]) {
      K[c].init(dim, 1);
      K[c](0, 0) = vertical_permeability[c];
    } else {
      K[c].init(dim, 2);
      for (int i=0; i<dim-1; i++) K[c](i, i) = horizontal_permeability[c];
      K[c](dim-1, dim-1) = vertical_permeability[c];
    }
  }
}


/* ******************************************************************
* Adds time derivative to cell-based part of MFD algebraic system.                                               
****************************************************************** */
void Richards_PK::addTimeDerivative_MFD(
    Epetra_Vector& pressure_cells, double dT_prec, Matrix_MFD* matrix)
{
  Epetra_Vector dSdP(mesh_->cell_map(false));
  derivedSdP(pressure_cells, dSdP);
 
  const Epetra_Vector& phi = FS->ref_porosity();
  std::vector<double>& Acc_cells = matrix->get_Acc_cells();
  std::vector<double>& Fc_cells = matrix->get_Fc_cells();

  int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);

  for (int c=0; c<ncells; c++) {
    double volume = mesh_->cell_volume(c);
    double factor = rho * phi[c] * dSdP[c] * volume / dT_prec;
    Acc_cells[c] += factor;
    Fc_cells[c] += factor * pressure_cells[c];
  }
}

}  // namespace AmanziFlow
}  // namespace Amanzi

