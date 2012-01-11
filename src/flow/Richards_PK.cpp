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
  super_map_ = create_super_map();

  // Other fundamental physical quantaties
  rho = *(FS->get_fluid_density());
  mu = *(FS->get_fluid_viscosity()); 
  gravity.init(dim);
  for (int k=0; k<dim; k++) gravity[k] = (*(FS->get_gravity()))[k];

#ifdef HAVE_MPI
  const  Epetra_Comm & comm = mesh_->cell_map(false).Comm(); 
  MyPID = comm.MyPID();

  const Epetra_Map& source_cmap = mesh_->cell_map(false);
  const Epetra_Map& target_cmap = mesh_->cell_map(true);

  cell_importer_ = Teuchos::rcp(new Epetra_Import(target_cmap, source_cmap));

  const Epetra_Map& source_fmap = mesh_->face_map(false);
  const Epetra_Map& target_fmap = mesh_->face_map(true);

  face_importer_ = Teuchos::rcp(new Epetra_Import(target_fmap, source_fmap));
#endif

  // miscalleneous
  flag_upwind = true;
  verbosity = FLOW_VERBOSITY_HIGH;
}


/* ******************************************************************
* Clean memory.
****************************************************************** */
Richards_PK::~Richards_PK() 
{ 
  delete super_map_; 
  delete solver; 
  if (matrix == preconditioner) {
    delete matrix; 
  } else {
    delete matrix;
    delete preconditioner;
  }
  delete bdf2_dae;
  delete bc_pressure;
  delete bc_head;
  delete bc_flux;
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

  // Create the solution vectors.
  solution = Teuchos::rcp(new Epetra_Vector(*super_map_));
  solution_cells = Teuchos::rcp(FS->createCellView(*solution));
  solution_faces = Teuchos::rcp(FS->createFaceView(*solution));
  rhs = Teuchos::rcp(new Epetra_Vector(*super_map_));

  solver = new AztecOO;
  solver->SetUserOperator(matrix);
  solver->SetPrecOperator(preconditioner);
  solver->SetAztecOption(AZ_solver, AZ_cg);  // symmetry is required

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
  K.resize(number_owned_cells);
  matrix->setSymmetryProperty(!flag_upwind);
  matrix->symbolicAssembleGlobalMatrices(*super_map_);

  // Preconditioner
  Teuchos::ParameterList ML_list = rp_list.sublist("ML Parameters");
  preconditioner->init_ML_preconditioner(ML_list); 

  // Create the Richards model evaluator and time integrator
  Teuchos::ParameterList& rme_list = rp_list.sublist("Richards model evaluator");
  absolute_tol_bdf = rme_list.get<double>("Absolute error tolerance", 1.0);
  relative_tol_bdf = rme_list.get<double>("Relative error tolerance", 1e-5); 

  Teuchos::RCP<Teuchos::ParameterList> bdf2_list(new Teuchos::ParameterList(rp_list.sublist("Time integrator")));
  bdf2_dae = new BDF2::Dae(*this, *super_map_);
  bdf2_dae->setParameterList(bdf2_list);

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
  int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  for (int c=0; c<ncells; c++) (*solution_cells)[c] = FS->ref_pressure()[c];

  if (method_sss == FLOW_STEADY_STATE_PICARD) {
    return advanceSteadyState_Picard();
  } else if (method_sss == FLOW_STEADY_STATE_BACKWARD_EULER) {
    //return advanceSteadyState_ForwardEuler();
    return advanceSteadyState_BackwardEuler();
  } else if (method_sss == FLOW_STEADY_STATE_BDF2) {
    return advanceSteadyState_BDF2();
  }
}


/* ******************************************************************* 
* Performs one time step of size dT. 
******************************************************************* */
int Richards_PK::advanceSteadyState_BDF2() 
{
  if (flag_upwind) solver->SetAztecOption(AZ_solver, AZ_bicgstab);  // symmetry is NOT required
  solver->SetAztecOption(AZ_output, AZ_none);

  double& time = (standalone_mode) ? T_internal : T_physical;
  int itrs = 0;

  while (itrs < max_itrs_sss) {
    calculateRelativePermeability(*solution_cells);
    setAbsolutePermeabilityTensor(K);
    if (flag_upwind) {  // Define K and Krel_faces
      calculateRelativePermeabilityUpwindGravity(*solution_cells);
      for (int c=0; c<K.size(); c++) K[c] *= rho / mu;
    } else {  // Define K and Krel_cells, Krel_faces is always one
      for (int c=0; c<K.size(); c++) K[c] *= (*Krel_cells)[c] * rho / mu;  
    }

    // update boundary conditions
    bc_pressure->Compute(time);
    bc_flux->Compute(time);
    bc_head->Compute(time);
    updateBoundaryConditions(bc_pressure, bc_head, bc_flux, bc_markers, bc_values);

    // create algebraic problem
    matrix->createMFDstiffnessMatrices(K, *Krel_faces);
    matrix->createMFDrhsVectors();
    addGravityFluxes_MFD(K, *Krel_faces, matrix);
    addTimeDerivative_MFD(*solution_cells, matrix);
    matrix->applyBoundaryConditions(bc_markers, bc_values);
    matrix->assembleGlobalMatrices();
    matrix->computeSchurComplement(bc_markers, bc_values);
    matrix->update_ML_preconditioner();

    T_physical = FS->get_time();
    double dTnext;

    if (itrs == 0) {  // initialization of BDF2
      Epetra_Vector udot(*super_map_);
      computeUDot(time, *solution, udot);
      bdf2_dae->set_initial_state(time, *solution, udot);

      int ierr;
      update_precon(time, *solution, dT, ierr);
    }

    bdf2_dae->bdf2_step(dT, 0.0, *solution, dTnext);
    bdf2_dae->commit_solution(dT, *solution);
    bdf2_dae->write_bdf2_stepping_statistics();

    T_physical = T_internal += dT;
    dT = dTnext;
    itrs++;

    double& time = (standalone_mode) ? T_internal : T_physical;
    time = bdf2_dae->most_recent_time();

    // DEBUG
    GMV::open_data_file(*mesh_, (std::string)"flow.gmv");
    GMV::start_data();
    GMV::write_cell_data(*solution_cells, 0, "pressure");
    GMV::write_cell_data(*Krel_cells, 0, "rel_permeability");
    GMV::close_data_file();
  }
  return 0;
}


/* ******************************************************************
* Calculates steady-state solution assuming that abosolute and relative
* permeabilities do not depend explicitly on time.                                                    
****************************************************************** */
int Richards_PK::advanceSteadyState_Picard()
{
  Epetra_Vector  solution_old(*solution);
  Epetra_Vector& solution_new = *solution;
  Epetra_Vector  residual(*solution);

  if (flag_upwind) solver->SetAztecOption(AZ_solver, AZ_bicgstab);  // symmetry is NOT required
  solver->SetAztecOption(AZ_output, AZ_none);

  int itrs = 0;
  double L2norm, L2error = 1.0, L2error_prev = 1.0;
  double relaxation = 0.5;

  while (L2error > err_tol_sss && itrs < max_itrs_sss) {
    calculateRelativePermeability(*solution_cells);
    setAbsolutePermeabilityTensor(K);
    if (flag_upwind) {  // Define K and Krel_faces
      calculateRelativePermeabilityUpwindGravity(*solution_cells);
      for (int c=0; c<K.size(); c++) K[c] *= rho / mu;
    } else {  // Define K and Krel_cells, Krel_faces is always one
      for (int c=0; c<K.size(); c++) K[c] *= (*Krel_cells)[c] * rho / mu;  
    }

    // create algebraic problem
    matrix->createMFDstiffnessMatrices(K, *Krel_faces);
    matrix->createMFDrhsVectors();
    addGravityFluxes_MFD(K, *Krel_faces, matrix);
    matrix->applyBoundaryConditions(bc_markers, bc_values);
    matrix->assembleGlobalMatrices();
    matrix->computeSchurComplement(bc_markers, bc_values);
    matrix->update_ML_preconditioner();

    // check convergence of non-linear residual
    rhs = matrix->get_rhs();
    L2error = matrix->computeResidual(solution_new, residual);

    // call AztecOO solver
    Epetra_Vector b(*rhs);
    solver->SetRHS(&b);  // AztecOO modifies the right-hand-side.
    solver->SetLHS(&*solution);  // initial solution guess 

    solver->Iterate(max_itrs, err_tol);
    num_itrs = solver->NumIters();
    double res_norm = solver->TrueResidual();

    // update relaxation parameter
    solution_new.Norm2(&L2norm);
    L2error /= L2norm;

    if (L2error > L2error_prev) relaxation = std::min(0.95, relaxation * 1.5); 
    else relaxation = std::max(0.05, relaxation * 0.9);  

    // information output
    if (verbosity >= FLOW_VERBOSITY_HIGH) {
      std::printf("Picard step:%4d   Pressure(res=%9.4e, sol=%9.4e, relax=%5.3f)  CG info(%8.3e, %4d)\n", 
          itrs, L2error, L2norm, relaxation, res_norm, num_itrs);
    }

    for (int c=0; c<number_owned_cells; c++) {
      solution_new[c] = relaxation * solution_old[c] + (1.0 - relaxation) * solution_new[c];
      solution_old[c] = solution_new[c];
    }

    L2error_prev = L2error;
    T_internal += dT;
    itrs++;

    // DEBUG
    GMV::open_data_file(*mesh_, (std::string)"flow.gmv");
    GMV::start_data();
    GMV::write_cell_data(*solution_cells, 0, "pressure");
    GMV::write_cell_data(*Krel_cells, 0, "rel_permeability");
    GMV::close_data_file();
  }
  
  Epetra_Vector& darcy_flux = FS->ref_darcy_flux();
  matrix->createMFDstiffnessMatrices(K, *Krel_faces);  // Should be improved. (lipnikov@lanl.gov)
  matrix->deriveDarcyFlux(*solution, *face_importer_, darcy_flux);
  addGravityFluxes_DarcyFlux(darcy_flux);

  return 0;
}


/* ******************************************************************
* Calculates steady-state solution assuming that abosolute and relative
* permeabilities do not depend explicitly on time.                                                    
****************************************************************** */
int Richards_PK::advanceSteadyState_BackwardEuler()
{
  Epetra_Vector  solution_old(*solution);
  Epetra_Vector& solution_new = *solution;

  if (flag_upwind) solver->SetAztecOption(AZ_solver, AZ_cgs);  // symmetry is NOT required
  solver->SetAztecOption(AZ_output, AZ_none);

  T_internal = 0.0;
  dT = dT0;

  int itrs = 0;
  double L2error = 1.0;
  while (L2error > err_tol_sss && itrs < max_itrs_sss) {
    calculateRelativePermeability(*solution_cells);
    setAbsolutePermeabilityTensor(K);
    if (flag_upwind) {  // Define K and Krel_faces
      calculateRelativePermeabilityUpwindGravity(*solution_cells);
      for (int c=0; c<K.size(); c++) K[c] *= rho / mu;
    } else {  // Define K and Krel_cells, Krel_faces is always one
      for (int c=0; c<K.size(); c++) K[c] *= (*Krel_cells)[c] * rho / mu;  
    }

    // update boundary conditions
    double time = T_internal;
    bc_pressure->Compute(time);
    bc_flux->Compute(time);
    bc_head->Compute(time);
    updateBoundaryConditions(bc_pressure, bc_head, bc_flux, bc_markers, bc_values);

    // create algebraic problem
    matrix->createMFDstiffnessMatrices(K, *Krel_faces);
    matrix->createMFDrhsVectors();
    addGravityFluxes_MFD(K, *Krel_faces, matrix);
    addTimeDerivative_MFD(*solution_cells, matrix);
    matrix->applyBoundaryConditions(bc_markers, bc_values);
    matrix->assembleGlobalMatrices();
    matrix->computeSchurComplement(bc_markers, bc_values);
    matrix->update_ML_preconditioner();

    // call AztecOO solver
    rhs = matrix->get_rhs();
    Epetra_Vector b(*rhs);
    solver->SetRHS(&b);  // AztecOO modifies the right-hand-side.
    solver->SetLHS(&*solution);  // initial solution guess 

    solver->Iterate(max_itrs, err_tol);
    num_itrs = solver->NumIters();
    double residual = solver->TrueResidual();

    // error estimates
    double sol_error = FS->norm_cell(solution_old, solution_new);
    double sol_norm = FS->norm_cell(solution_new);
    L2error = sol_error / sol_norm;

    if (residual > sol_norm * 1e-3) {
      dT /= 4;
      solution_new = solution_old;
      if (verbosity >= FLOW_VERBOSITY_HIGH) {
        std::printf("Declined: %4d  Pressure(diff=%9.4e, sol=%9.4e)  CG info(%8.3e, %4d),  dT=%8.3e\n", 
            itrs, L2error, sol_norm, residual, num_itrs, dT);
      }
    } else { 
      dT = std::min(dT0, dT * 2);
      solution_old = solution_new;

      if (verbosity >= FLOW_VERBOSITY_HIGH) {
        std::printf("Time step:%4d  Pressure(diff=%9.4e, sol=%9.4e)  CG info(%8.3e, %4d),  dT=%8.3e\n", 
            itrs, L2error, sol_norm, residual, num_itrs, dT);
      }

      T_internal += dT;
      itrs++;
    }

    // DEBUG
    GMV::open_data_file(*mesh_, (std::string)"flow.gmv");
    GMV::start_data();
    GMV::write_cell_data(*solution_cells, 0, "pressure");
    GMV::write_cell_data(*(FS->get_absolute_permeability()), 0, "abs_permeability");
    GMV::write_cell_data(*Krel_cells, 0, "rel_permeability");
    GMV::close_data_file();
  }
  
  Epetra_Vector& darcy_flux = FS->ref_darcy_flux();
  matrix->createMFDstiffnessMatrices(K, *Krel_faces);  // Should be improved. (lipnikov@lanl.gov)
  matrix->deriveDarcyFlux(*solution, *face_importer_, darcy_flux);
  addGravityFluxes_DarcyFlux(darcy_flux);

  return 0;
}


/* ******************************************************************
* Calculates steady-state solution assuming that abosolute and relative
* permeabilities do not depend explicitly on time.                                                    
****************************************************************** */
int Richards_PK::advanceSteadyState_ForwardEuler()
{
  Epetra_Vector solution_new(*solution);

  T_internal = 0.0;
  dT = 1.0;

  int itrs = 0;
  double L2error = 1.0;
  while (L2error > err_tol_sss && itrs < 20) {
    calculateRelativePermeability(*solution_cells);
    setAbsolutePermeabilityTensor(K);
    if (flag_upwind) {  // Define K and Krel_faces
      calculateRelativePermeabilityUpwindGravity(*solution_cells);
      for (int c=0; c<K.size(); c++) K[c] *= rho / mu;
    } else {  // Define K and Krel_cells, Krel_faces is always one
      for (int c=0; c<K.size(); c++) K[c] *= (*Krel_cells)[c] * rho / mu;  
    }

    // create algebraic problem
    matrix->createMFDstiffnessMatrices(K, *Krel_faces);
    matrix->createMFDrhsVectors();
    addGravityFluxes_MFD(K, *Krel_faces, matrix);
    matrix->applyBoundaryConditions(bc_markers, bc_values);
    matrix->assembleGlobalMatrices();

    // calculate solution at ne time step
    rhs = matrix->get_rhs();
    L2error = matrix->computeResidual(*solution, solution_new);
    solution_new.Update(1.0, *solution, dT);

    // error estimates
    double sol_error = 0.0, sol_norm = 0.0;
    for (int n=0; n<solution->MyLength(); n++) {
      sol_error += std::pow(solution_new[n] - (*solution)[n], 2.0);
      sol_norm += std::pow(solution_new[n], 2.0);
    }
    sol_error = sqrt(sol_error / sol_norm);
    sol_norm = sqrt(sol_norm);
 
    *solution = solution_new;
    if (verbosity >= FLOW_VERBOSITY_HIGH) {
      std::printf("Time step:%6d   Pressure(diff=%9.4e, sol=%9.4e, res=%9.4e)  time=%8.3e\n", 
          itrs, sol_error, sol_norm, L2error, T_internal);
    }

    T_internal += dT;
    itrs++;

    // DEBUG
    GMV::open_data_file(*mesh_, (std::string)"flow.gmv");
    GMV::start_data();
    GMV::write_cell_data(*solution_cells, 0, "pressure");
    GMV::write_cell_data(*(FS->get_absolute_permeability()), 0, "abs_permeability");
    GMV::write_cell_data(*Krel_cells, 0, "rel_permeability");
    GMV::close_data_file();
  }
  return 0;
}


/* ******************************************************************* 
* Performs one time step of size dT. 
******************************************************************* */
int Richards_PK::advance(double dT) 
{
  solver->SetAztecOption(AZ_output, AZ_none);

  calculateRelativePermeability(*solution_cells);
  if (flag_upwind) calculateRelativePermeabilityUpwindGravity(*solution_cells);
 
  // work-around limited support for tensors
  setAbsolutePermeabilityTensor(K);  
  for (int c=0; c<K.size(); c++) K[c] *= (*Krel_cells)[c] * rho / mu;

  matrix->createMFDstiffnessMatrices(K, *Krel_faces);
  matrix->createMFDrhsVectors();
  addGravityFluxes_MFD(K, *Krel_faces, matrix);
  addTimeDerivative_MFD(*solution_cells, matrix);
  matrix->applyBoundaryConditions(bc_markers, bc_values);
  matrix->assembleGlobalMatrices();
  matrix->computeSchurComplement(bc_markers, bc_values);
  matrix->update_ML_preconditioner();

  T_physical = FS->get_time();
  double time = (standalone_mode) ? T_internal : T_physical;
  double dTnext;

  if (itrs_bdf2 == 0) {  // initialization
    Epetra_Vector udot(*super_map_);
    computeUDot(time, *solution, udot);
    bdf2_dae->set_initial_state(time, *solution, udot);

    int ierr;
    update_precon(time, *solution, dT, ierr);
  }

  bdf2_dae->bdf2_step(dT, 0.0, *solution, dTnext);
  bdf2_dae->commit_solution(dT, *solution);
  bdf2_dae->write_bdf2_stepping_statistics();

  itrs_bdf2++;
  return 0;
}


/* ******************************************************************
* Estimate du/dt from the pressure equations, u is the global
* solution vector.                                            
****************************************************************** */
double Richards_PK::computeUDot(const double T, const Epetra_Vector& u, Epetra_Vector &udot)
{
  double norm_udot;
  norm_udot = matrix->computeResidual(u, udot);

  Epetra_Vector* udot_faces = FS->createFaceView(udot);
  udot_faces->PutScalar(0.0);

  return norm_udot;
}


/* ******************************************************************
* Gathers together routines to compute MFD matrices Axx(u) and 
* preconditioner Sff(u) using internal time step dT.                             
****************************************************************** */
void Richards_PK::computePreconditionerMFD(const Epetra_Vector &u, Matrix_MFD* matrix)
{
  Epetra_Vector* u_cells = FS->createCellView(u);

  // setup absolute and compute relative permeabilities
  calculateRelativePermeability(*u_cells);
  setAbsolutePermeabilityTensor(K);
  if (flag_upwind) {  // Define K and Krel_faces
    calculateRelativePermeabilityUpwindGravity(*solution_cells);
    for (int c=0; c<K.size(); c++) K[c] *= rho / mu;
  } else {  // Define K and Krel_cells, Krel_faces is always one
    for (int c=0; c<K.size(); c++) K[c] *= (*Krel_cells)[c] * rho / mu;  
  }

  // setup new algebraic problems
  matrix->createMFDstiffnessMatrices(K, *Krel_faces);
  matrix->createMFDrhsVectors();
  addGravityFluxes_MFD(K, *Krel_faces, matrix);
  addTimeDerivative_MFD(*solution_cells, matrix);
  matrix->applyBoundaryConditions(bc_markers, bc_values);
  matrix->assembleGlobalMatrices();
  matrix->computeSchurComplement(bc_markers, bc_values);
  matrix->update_ML_preconditioner();  
}


/* ******************************************************************
* Temporary convertion from double to tensor.                                               
****************************************************************** */
void Richards_PK::setAbsolutePermeabilityTensor(std::vector<WhetStone::Tensor>& K)
{
  const Epetra_Vector& permeability = FS->ref_absolute_permeability();

  for (int c=cmin; c<=cmax; c++) {
    K[c].init(dim, 1);
    K[c](0, 0) = permeability[c];
  }
}


/* ******************************************************************
* Adds time derivative to cell-based part of MFD algebraic system.                                               
****************************************************************** */
void Richards_PK::addTimeDerivative_MFD(Epetra_Vector& pressure_cells, Matrix_MFD* matrix)
{
  Epetra_Vector dSdP(mesh_->cell_map(false));
  derivedSdP(pressure_cells, dSdP);
 
  const Epetra_Vector& phi = FS->ref_porosity();
  std::vector<double>& Acc_cells = matrix->get_Acc_cells();
  std::vector<double>& Fc_cells = matrix->get_Fc_cells();

  int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);

  for (int c=0; c<ncells; c++) {
    double volume = mesh_->cell_volume(c);
    double factor = dSdP[c] * phi[c] * volume / dT;
    Acc_cells[c] += factor;
    Fc_cells[c] += factor * pressure_cells[c];
  }
}

}  // namespace AmanziFlow
}  // namespace Amanzi

