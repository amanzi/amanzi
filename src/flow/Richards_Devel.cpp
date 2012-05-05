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

#include "Richards_PK.hpp"

namespace Amanzi {
namespace AmanziFlow {

/* ******************************************************************
* Makes one Picard step during transient time integration.
* This is the experimental method.                                                 
****************************************************************** */
int Richards_PK::PicardStep(double Tp, double dTp, double& dTnext)
{
  Epetra_Vector solution_old(*solution);
  Epetra_Vector solution_new(*solution);

  Teuchos::RCP<Epetra_Vector> solution_old_cells = Teuchos::rcp(FS->createCellView(solution_old));
  Teuchos::RCP<Epetra_Vector> solution_new_cells = Teuchos::rcp(FS->createCellView(solution_new));

  if (!is_matrix_symmetric) solver->SetAztecOption(AZ_solver, AZ_gmres);
  solver->SetAztecOption(AZ_output, AZ_none);

  int itrs;
  for (itrs = 0; itrs < 20; itrs++) {
    SetAbsolutePermeabilityTensor(K);

    if (!is_matrix_symmetric) {  // Define K and Krel_faces
      CalculateRelativePermeabilityFace(*solution_old_cells);
      for (int c = 0; c < K.size(); c++) K[c] *= rho / mu;
    } else {  // Define K and Krel_cells, Krel_faces is always one
      CalculateRelativePermeabilityCell(*solution_old_cells);
      for (int c = 0; c < K.size(); c++) K[c] *= (*Krel_cells)[c] * rho / mu;
    }

    // update boundary conditions
    double time = Tp + dTp;
    bc_pressure->Compute(time);
    bc_flux->Compute(time);
    bc_head->Compute(time);
    bc_seepage->Compute(time);
    UpdateBoundaryConditions(
        bc_pressure, bc_head, bc_flux, bc_seepage,
        *solution_old_cells, atm_pressure,
        bc_markers, bc_values);

    // create algebraic problem
    matrix->createMFDstiffnessMatrices(mfd3d_method, K, *Krel_faces);
    matrix->createMFDrhsVectors();
    addGravityFluxes_MFD(K, *Krel_faces, matrix);
    AddTimeDerivative_MFDpicard(*solution_cells, *solution_old_cells, dTp, matrix);
    matrix->applyBoundaryConditions(bc_markers, bc_values);
    matrix->assembleGlobalMatrices();
    rhs = matrix->get_rhs();

    // create preconditioner
    int disc_method = AmanziFlow::FLOW_MFD3D_TWO_POINT_FLUX;
    preconditioner->createMFDstiffnessMatrices(disc_method, K, *Krel_faces);
    preconditioner->createMFDrhsVectors();
    addGravityFluxes_MFD(K, *Krel_faces, preconditioner);
    preconditioner->applyBoundaryConditions(bc_markers, bc_values);
    preconditioner->assembleGlobalMatrices();
    preconditioner->computeSchurComplement(bc_markers, bc_values);
    preconditioner->update_ML_preconditioner();

    // call AztecOO solver
    solver->SetRHS(&*rhs);  // AztecOO modifies the right-hand-side.
    solution_new = solution_old;  // initial solution guess
    solver->SetLHS(&solution_new);

    solver->Iterate(max_itrs, convergence_tol);
    int num_itrs = solver->NumIters();
    double linear_residual = solver->TrueResidual();

    int ndof = ncells_owned + nfaces_owned;
    for (int c = 0; c < ndof; c++) {
      solution_new[c] = (solution_old[c] + solution_new[c]) / 2;
    }

    // error analysis
    double error = ErrorNorm(solution_old, solution_new);

    if (MyPID == 0 && verbosity >= FLOW_VERBOSITY_HIGH) {
      std::printf("Picard:%4d   Pressure(error=%9.4e)  solver(%8.3e, %4d)\n",
          itrs, error, linear_residual, num_itrs);
    }

    if (error < 10.0) 
      break;
    else 
      solution_old = solution_new;
  }

  if (itrs < 10) {
    *solution = solution_new;
    dTnext = dT * 1.25;
    return itrs;
  } else if (itrs < 15) {
    *solution = solution_new;
    dTnext = dT;
    return itrs;
  } else if (itrs < 19) {
    *solution = solution_new;
    dTnext = dT * 0.8;
    *solution = solution_new;
  }

  throw itrs;
}


/* ******************************************************************
* Calculates steady-state solution assuming that abosolute and relative
* permeabilities do not depend explicitly on time.
* WARNING: temporary replacement for missing BDF1 time integrator.                                                    
****************************************************************** */
int Richards_PK::AdvanceSteadyState_BackwardEuler()
{
  Epetra_Vector  solution_old(*solution);
  Epetra_Vector& solution_new = *solution;

  if (is_matrix_symmetric) solver->SetAztecOption(AZ_solver, AZ_cgs);
  solver->SetAztecOption(AZ_output, AZ_none);

  T_internal = T0_sss;
  dT = dT0_sss;

  int itrs = 0, ifail = 0;
  double L2error = 1.0;
  while (L2error > convergence_tol_sss && itrs < max_itrs_sss) {
    SetAbsolutePermeabilityTensor(K);

    if (!is_matrix_symmetric) {  // Define K and Krel_faces
      CalculateRelativePermeabilityFace(*solution_cells);
      for (int c = 0; c < K.size(); c++) K[c] *= rho / mu;
    } else {  // Define K and Krel_cells, Krel_faces is always one
      CalculateRelativePermeabilityCell(*solution_cells);
      for (int c = 0; c < K.size(); c++) K[c] *= (*Krel_cells)[c] * rho / mu;
    }

    // update boundary conditions
    double time = T_internal + dT;
    bc_pressure->Compute(time);
    bc_flux->Compute(time);
    bc_head->Compute(time);
    bc_seepage->Compute(time);
    UpdateBoundaryConditions(
        bc_pressure, bc_head, bc_flux, bc_seepage,
        *solution_cells, atm_pressure,
        bc_markers, bc_values);

    // create algebraic problem
    matrix->createMFDstiffnessMatrices(mfd3d_method, K, *Krel_faces);
    matrix->createMFDrhsVectors();
    addGravityFluxes_MFD(K, *Krel_faces, matrix);
    AddTimeDerivative_MFDfake(*solution_cells, dT, matrix);
    matrix->applyBoundaryConditions(bc_markers, bc_values);
    matrix->assembleGlobalMatrices();

    // create preconditioner
    int disc_method = AmanziFlow::FLOW_MFD3D_TWO_POINT_FLUX;
    preconditioner->createMFDstiffnessMatrices(disc_method, K, *Krel_faces);
    preconditioner->createMFDrhsVectors();
    addGravityFluxes_MFD(K, *Krel_faces, preconditioner);
    preconditioner->applyBoundaryConditions(bc_markers, bc_values);
    preconditioner->assembleGlobalMatrices();
    preconditioner->computeSchurComplement(bc_markers, bc_values);
    preconditioner->update_ML_preconditioner();

    // call AztecOO solver
    rhs = matrix->get_rhs();
    Epetra_Vector b(*rhs);
    solver->SetRHS(&b);  // AztecOO modifies the right-hand-side.
    solver->SetLHS(&*solution);  // initial solution guess

    solver->Iterate(max_itrs, convergence_tol);
    int num_itrs = solver->NumIters();
    double residual = solver->TrueResidual();

    // error estimates
    double sol_norm = FS->normLpCell(solution_new, 2.0);
    L2error = ErrorNorm(solution_old, solution_new);

    if (L2error > 1.0 && itrs && ifail < 5) {  // itrs=0 allows to avoid bad initial guess.
      dT /= 10;
      solution_new = solution_old;
      if (MyPID == 0 && verbosity >= FLOW_VERBOSITY_HIGH) {
        std::printf("Fail:%4d  Pressure(diff=%9.4e, sol=%9.4e)  solver(%8.3e,%3d), T=%9.3e dT=%7.2e\n",
            itrs, L2error, sol_norm, residual, num_itrs, T_internal, dT);
      }
      ifail++;
    } else {
      T_internal += dT;
      solution_old = solution_new;

      if (MyPID == 0 && verbosity >= FLOW_VERBOSITY_HIGH) {
        std::printf("Step:%4d  Pressure(diff=%9.4e, sol=%9.4e)  solver(%8.3e,%3d), T=%9.3e dT=%7.2e\n",
            itrs, L2error, sol_norm, residual, num_itrs, T_internal, dT);
      }

      ifail = 0;
      dT = std::min(dT*1.25, dTmax_sss);
      itrs++;
    }

    if (T_internal > T1_sss) break;
  }

  num_nonlinear_steps = itrs;
  return 0;
}


/* ******************************************************************
* Adds time derivative to the cell-based part of MFD algebraic system.                                               
****************************************************************** */
void Richards_PK::AddTimeDerivative_MFDpicard(
    Epetra_Vector& pressure_cells, Epetra_Vector& pressure_cells_dSdP, 
    double dT_prec, Matrix_MFD* matrix)
{
  Epetra_Vector dSdP(mesh_->cell_map(false));
  DerivedSdP(pressure_cells_dSdP, dSdP);

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
* Adds time derivative to cell-based part of MFD algebraic system.                                               
****************************************************************** */
void Richards_PK::AddTimeDerivative_MFDfake(
    Epetra_Vector& pressure_cells, double dT_prec, Matrix_MFD* matrix)
{
  std::vector<double>& Acc_cells = matrix->get_Acc_cells();
  std::vector<double>& Fc_cells = matrix->get_Fc_cells();

  for (int c = 0; c < ncells_owned; c++) {
    double volume = mesh_->cell_volume(c);
    double factor = volume / dT_prec;
    Acc_cells[c] += factor;
    Fc_cells[c] += factor * pressure_cells[c];
  }
}


/* ******************************************************************
* Check difference between solutions at times T and T+dT.                                                 
****************************************************************** */
double Richards_PK::ErrorNorm(const Epetra_Vector& uold, const Epetra_Vector& unew)
{
  double error_norm = 0.0;
  for (int n = 0; n < ncells_owned; n++) {
    double tmp = abs(uold[n] - unew[n]) / (absolute_tol + relative_tol * abs(uold[n]));
    error_norm = std::max<double>(error_norm, tmp);
  }

#ifdef HAVE_MPI
  double buf = error_norm;
  uold.Comm().MaxAll(&buf, &error_norm, 1);  // find the global maximum
#endif
  return error_norm;
}

}  // namespace AmanziFlow
}  // namespace Amanzi

