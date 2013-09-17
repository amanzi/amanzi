/*
This is the flow component of the Amanzi code. 

Copyright 2010-2012 held jointly by LANS/LANL, LBNL, and PNNL. 
Amanzi is released under the three-clause BSD License. 
The terms of use and "as is" disclaimer for this license are 
provided in the top-level COPYRIGHT file.

Author:  Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <algorithm>
#include <vector>

#include "LinearOperatorFactory.hh"
#include "Richards_PK.hh"

namespace Amanzi {
namespace AmanziFlow {

/* ******************************************************************
*  Wrapper for advance to steady-state routines.                                                    
****************************************************************** */
int Richards_PK::AdvanceToSteadyState(double T0, double dT0)
{
  // update the axiliary flow state
  FS->CopyMasterFace2GhostFace(FS->ref_darcy_flux(), FS_aux->ref_darcy_flux());

  // override internal parameters
  T_physics = ti_specs_sss_.T0 = T0;
  dT = ti_specs_sss_.dT0 = dT0;

  int ierr = 0;
  int ti_method = ti_specs_sss_.ti_method;
  if (ti_method == FLOW_TIME_INTEGRATION_PICARD) {
    ierr = AdvanceToSteadyState_Picard(ti_specs_sss_);
  } else if (ti_method == FLOW_TIME_INTEGRATION_BACKWARD_EULER) {
    ierr = AdvanceToSteadyState_BackwardEuler(ti_specs_sss_);
  } else if (ti_method == FLOW_TIME_INTEGRATION_BDF1) {
    ierr = AdvanceToSteadyState_BDF1(ti_specs_sss_);
  } else if (ti_method == FLOW_TIME_INTEGRATION_BDF2) {
    ierr = AdvanceToSteadyState_BDF2(ti_specs_sss_);
  }

  Epetra_Vector& ws = FS->ref_water_saturation();
  Epetra_Vector& ws_prev = FS->ref_prev_water_saturation();
  DeriveSaturationFromPressure(*solution_cells, ws);
  ws_prev = ws;

  return ierr;
}


/* ******************************************************************* 
* Performs one time step of size dT using first-order time integrator.
******************************************************************* */
int Richards_PK::AdvanceToSteadyState_BDF1(TI_Specs& ti_specs)
{
  bool last_step = false;

  int max_itrs = ti_specs.max_itrs;
  double T0 = ti_specs.T0;
  double T1 = ti_specs.T1;
  double dT0 = ti_specs.dT0;

  int itrs = 0;
  while (itrs < max_itrs && T_physics < T1) {
    if (itrs == 0) {  // initialization of BDF1
      Epetra_Vector udot(*super_map_);
      ComputeUDot(T0, *solution, udot);
      bdf1_dae->set_initial_state(T0, *solution, udot);

      int ierr;
      update_precon(T0, *solution, dT0, ierr);
    }

    double dTnext;
    try { 
      bdf1_dae->bdf1_step(dT, *solution, dTnext);
    } catch (int i) {
      dT /= 2;
      continue;
    }
    bdf1_dae->commit_solution(dT, *solution);
    bdf1_dae->write_bdf1_stepping_statistics();

    T_physics = bdf1_dae->most_recent_time();
    dT = dTnext;
    itrs++;

    double Tdiff = T1 - T_physics;
    if (dTnext > Tdiff) {
      dT = Tdiff * 0.99999991;  // To avoid hitting the wrong BC
      last_step = true;
    }
    if (last_step && dT < 1e-3) break;
  }

  ti_specs.num_itrs = itrs;
  return 0;
}


/* ******************************************************************* 
* Performs one time step of size dT using second-order time integrator.
******************************************************************* */
int Richards_PK::AdvanceToSteadyState_BDF2(TI_Specs& ti_specs)
{
  bool last_step = false;

  int max_itrs_nonlinear = ti_specs.max_itrs;
  double T0 = ti_specs.T0;
  double T1 = ti_specs.T1;
  double dT0 = ti_specs.dT0;

  int itrs = 0;
  while (itrs < max_itrs_nonlinear && T_physics < T1) {
    if (itrs == 0) {  // initialization of BDF2
      Epetra_Vector udot(*super_map_);
      ComputeUDot(T0, *solution, udot);
      bdf2_dae->set_initial_state(T0, *solution, udot);

      int ierr;
      update_precon(T0, *solution, dT0, ierr);
    }

    double dTnext;
    try {
      bdf2_dae->bdf2_step(dT, 0.0, *solution, dTnext);
    } catch (int i) {
      dT /= 2;
      continue;
    }
    bdf2_dae->commit_solution(dT, *solution);
    bdf2_dae->write_bdf2_stepping_statistics();

    T_physics = bdf2_dae->most_recent_time();
    dT = dTnext;
    itrs++;

    double Tdiff = T1 - T_physics;
    if (dTnext > Tdiff) {
      dT = Tdiff * 0.99999991;  // To avoid hitting the wrong BC
      last_step = true;
    }
    if (last_step && dT < 1e-3) break;
  }

  ti_specs.num_itrs = itrs;
  return 0;
}


/* ******************************************************************
* Calculates steady-state solution assuming that absolute and
* relative permeabilities do not depend explicitly on time.
* This is the experimental method.                                                 
****************************************************************** */
int Richards_PK::AdvanceToSteadyState_Picard(TI_Specs& ti_specs)
{
  Epetra_Vector  solution_old(*solution);
  Epetra_Vector& solution_new = *solution;
  Epetra_Vector  residual(*solution);

  Epetra_Vector* solution_old_cells = FS->CreateCellView(solution_old);
  Epetra_Vector* solution_new_cells = FS->CreateCellView(solution_new);

  // update steady state boundary conditions
  double time = T_physics;
  bc_pressure->Compute(time);
  bc_flux->Compute(time);
  if (shift_water_table_.getRawPtr() == NULL)
    bc_head->Compute(time);
  else
    bc_head->ComputeShift(time, shift_water_table_->Values());

  // update steady state source conditons
  if (src_sink != NULL) {
    if (src_sink_distribution & Amanzi::Functions::DOMAIN_FUNCTION_ACTION_DISTRIBUTE_PERMEABILITY) {
      src_sink->ComputeDistribute(time, Kxy->Values()); 
    } else {
      src_sink->ComputeDistribute(time, NULL);
    }
  }

  int max_itrs_nonlinear = ti_specs.max_itrs;
  double residual_tol_nonlinear = ti_specs.residual_tol;
  LinearSolver_Specs& ls_specs = ti_specs.ls_specs;

  int itrs = 0;
  double L2norm, L2error = 1.0;

  while (L2error > residual_tol_nonlinear && itrs < max_itrs_nonlinear) {
    // update dynamic boundary conditions
    bc_seepage->Compute(time);
    ProcessBoundaryConditions(
        bc_pressure, bc_head, bc_flux, bc_seepage,
        *solution_cells, *solution_faces, atm_pressure,
        rainfall_factor, bc_submodel, bc_model, bc_values);

    // update permeabilities
    rel_perm->Compute(*solution, bc_model, bc_values);

    // create algebraic problem
    matrix_->CreateMFDstiffnessMatrices(*rel_perm);
    matrix_->CreateMFDrhsVectors();
    AddGravityFluxes_MFD(K, &*matrix_, *rel_perm);
    matrix_->ApplyBoundaryConditions(bc_model, bc_values);
    matrix_->AssembleGlobalMatrices();
    rhs = matrix_->rhs();  // export RHS from the matrix class
    if (src_sink != NULL) AddSourceTerms(src_sink, *rhs);

    // create preconditioner
    preconditioner_->CreateMFDstiffnessMatrices(*rel_perm);
    preconditioner_->CreateMFDrhsVectors();
    preconditioner_->ApplyBoundaryConditions(bc_model, bc_values);
    preconditioner_->AssembleSchurComplement(bc_model, bc_values);
    preconditioner_->UpdatePreconditioner();

    // check convergence of non-linear residual
    L2error = matrix_->ComputeResidual(solution_new, residual);
    residual.Norm2(&L2error);
    rhs->Norm2(&L2norm);
    L2error /= L2norm;

    // solve linear problem
    AmanziSolvers::LinearOperatorFactory<Matrix_MFD, Epetra_Vector, Epetra_Map> factory;
    Teuchos::RCP<AmanziSolvers::LinearOperator<Matrix_MFD, Epetra_Vector, Epetra_Map> >
       solver = factory.Create(ls_specs.solver_name, solver_list_, matrix_, preconditioner_);

    solver->ApplyInverse(*rhs, *solution);

    int num_itrs_linear = solver->num_itrs();
    double linear_residual = solver->residual();

    // update relaxation
    double relaxation;
    relaxation = CalculateRelaxationFactor(*solution_old_cells, *solution_new_cells);

    if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH) {
      Teuchos::OSTab tab = vo_->getOSTab();
      *(vo_->os()) << "Picard:" << itrs << " ||r||=" << L2error << " relax=" << relaxation 
                   << " solver(" << linear_residual << ", " << num_itrs_linear << ")" << endl;
    }

// Epetra_Vector& pressure = FS->ref_pressure();
// for (int c = 0; c < ncells_owned; c++) pressure[c] = solution_old[c] - solution_new[c];
// WriteGMVfile(FS); 
    int ndof = ncells_owned + nfaces_owned;
    for (int c = 0; c < ndof; c++) {
      solution_new[c] = (1.0 - relaxation) * solution_old[c] + relaxation * solution_new[c];
      solution_old[c] = solution_new[c];
    }

    T_physics += dT;
    itrs++;
  }

  ti_specs.num_itrs = itrs;
  return 0;
}


/* ******************************************************************
* Calculates steady-state solution using the Picard Newton method.                                             
****************************************************************** */
int Richards_PK::AdvanceToSteadyState_PicardNewton(TI_Specs& ti_specs)
{
  Epetra_Vector  solution_old(*solution);
  Epetra_Vector& solution_new = *solution;
  Epetra_Vector  residual(*solution);

  Epetra_Vector* solution_old_cells = FS->CreateCellView(solution_old);
  Epetra_Vector* solution_new_cells = FS->CreateCellView(solution_new);

  Epetra_Vector& flux = FS->ref_darcy_flux();
  Epetra_Vector  flux_new(flux);

  // update steady state boundary conditions
  double time = T_physics;
  bc_pressure->Compute(time);
  bc_flux->Compute(time);
  if (shift_water_table_.getRawPtr() == NULL)
    bc_head->Compute(time);
  else
    bc_head->ComputeShift(time, shift_water_table_->Values());

  // convergence criteria
  int max_itrs_nonlinear = ti_specs.max_itrs;
  double residual_tol_nonlinear = ti_specs.residual_tol;
  LinearSolver_Specs& ls_specs = ti_specs.ls_specs;

  int itrs = 0;
  double L2norm, L2error = 1.0;

  while (L2error > residual_tol_nonlinear && itrs < max_itrs_nonlinear) {
    // update dynamic boundary conditions
    bc_seepage->Compute(time);
    ProcessBoundaryConditions(
        bc_pressure, bc_head, bc_flux, bc_seepage,
        *solution_cells, *solution_faces, atm_pressure,
        rainfall_factor, bc_submodel, bc_model, bc_values);

    // update permeabilities
    rel_perm->Compute(*solution, bc_model, bc_values);

    // create algebraic problem
    rhs = matrix_->rhs();  // export RHS from the matrix class
    matrix_->CreateMFDstiffnessMatrices(*rel_perm);
    matrix_->CreateMFDrhsVectors();
    AddGravityFluxes_MFD(K, &*matrix_, *rel_perm);
    // AddNewtonFluxes_MFD(*rel_perm, *solution_cells, flux, *rhs,
    //                     static_cast<Matrix_MFD_PLambda*>(matrix_));
    matrix_->ApplyBoundaryConditions(bc_model, bc_values);
    matrix_->AssembleGlobalMatrices();

    // create preconditioner
    preconditioner_->CreateMFDstiffnessMatrices(*rel_perm);
    preconditioner_->CreateMFDrhsVectors();
    // AddNewtonFluxes_MFD(*rel_perm, *solution_cells, flux, residual,
    //                     static_cast<Matrix_MFD_PLambda*>(preconditioner_));
    preconditioner_->ApplyBoundaryConditions(bc_model, bc_values);
    preconditioner_->AssembleSchurComplement(bc_model, bc_values);
    preconditioner_->UpdatePreconditioner();

    // check convergence of non-linear residual
    L2error = matrix_->ComputeResidual(solution_new, residual);
    residual.Norm2(&L2error);
    rhs->Norm2(&L2norm);
    L2error /= L2norm;

    // solve linear problem
    AmanziSolvers::LinearOperatorFactory<Matrix_MFD, Epetra_Vector, Epetra_Map> factory;
    Teuchos::RCP<AmanziSolvers::LinearOperator<Matrix_MFD, Epetra_Vector, Epetra_Map> >
       solver = factory.Create(ls_specs.solver_name, solver_list_, matrix_, preconditioner_);

    solver->ApplyInverse(*rhs, *solution);

    int num_itrs_linear = solver->num_itrs();
    double linear_residual = solver->residual();

    // update relaxation
    double relaxation;
    relaxation = CalculateRelaxationFactor(*solution_old_cells, *solution_new_cells);

    if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH) {
      Teuchos::OSTab tab = vo_->getOSTab();
      *(vo_->os()) << "Picard:" << itrs << " ||r||=" << L2error << " relax=" << relaxation 
                   << " solver(" << linear_residual << ", " << num_itrs_linear << ")" << endl;
    }

    matrix_->DeriveDarcyMassFlux(*solution, *face_importer_, flux_new);
    AddGravityFluxes_DarcyFlux(K, flux_new, *rel_perm);

    int ndof = ncells_owned + nfaces_owned;
    for (int c = 0; c < ndof; c++) {
      solution_new[c] = (1.0 - relaxation) * solution_old[c] + relaxation * solution_new[c];
      solution_old[c] = solution_new[c];
    }

    for (int f = 0; f < nfaces_owned; f++) {
      flux[f] = (1.0 - relaxation) * flux[f] + relaxation * flux_new[f];
    }

    // DEBUG
    CommitState(FS); WriteGMVfile(FS);

    T_physics += dT;
    itrs++;
  }

  ti_specs.num_itrs = itrs;
  return 0;
}


/* ******************************************************************
* Calculate relazation factor.                                                       
****************************************************************** */
double Richards_PK::CalculateRelaxationFactor(const Epetra_Vector& uold, const Epetra_Vector& unew)
{ 
  double relaxation = 1.0;

  if (error_control_ & FLOW_TI_ERROR_CONTROL_SATURATION) {
    Epetra_Vector dSdP(mesh_->cell_map(false));
    rel_perm->DerivedSdP(uold, dSdP);

    for (int c = 0; c < ncells_owned; c++) {
      double diff = dSdP[c] * fabs(unew[c] - uold[c]);
      if (diff > 3e-2) relaxation = std::min(relaxation, 3e-2 / diff);
    }
  }

  if (error_control_ & FLOW_TI_ERROR_CONTROL_PRESSURE) {
    for (int c = 0; c < ncells_owned; c++) {
      double diff = fabs(unew[c] - uold[c]);
      double umax = std::max(fabs(unew[c]), fabs(uold[c]));
      if (diff > 1e-2 * umax) relaxation = std::min(relaxation, 1e-2 * umax / diff);
    }
  }

#ifdef HAVE_MPI
  double relaxation_tmp = relaxation;
  mesh_->get_comm()->MinAll(&relaxation_tmp, &relaxation, 1);  // find the global minimum
#endif

  return relaxation;
}

}  // namespace AmanziFlow
}  // namespace Amanzi

