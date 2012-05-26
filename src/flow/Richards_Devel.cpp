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

  Teuchos::RCP<Epetra_Vector> solution_old_cells = Teuchos::rcp(FS->CreateCellView(solution_old));
  Teuchos::RCP<Epetra_Vector> solution_old_faces = Teuchos::rcp(FS->CreateFaceView(solution_old));
  Teuchos::RCP<Epetra_Vector> solution_new_cells = Teuchos::rcp(FS->CreateCellView(solution_new));

  if (!is_matrix_symmetric) solver->SetAztecOption(AZ_solver, AZ_gmres);
  solver->SetAztecOption(AZ_output, AZ_none);
  solver->SetAztecOption(AZ_conv, AZ_rhs);

  int itrs;
  for (itrs = 0; itrs < 20; itrs++) {
    if (!is_matrix_symmetric) {
      CalculateRelativePermeabilityFace(*solution_old_cells);
      Krel_cells->PutScalar(1.0);
    } else {
      CalculateRelativePermeabilityCell(*solution_old_cells);
      Krel_faces->PutScalar(1.0);
    }

    // update boundary conditions
    double time = Tp + dTp;
    bc_pressure->Compute(time);
    bc_flux->Compute(time);
    bc_head->Compute(time);
    bc_seepage->Compute(time);
    UpdateBoundaryConditions(
        bc_pressure, bc_head, bc_flux, bc_seepage,
        *solution_old_faces, atm_pressure,
        bc_markers, bc_values);

    // create algebraic problem
    matrix->CreateMFDstiffnessMatrices(*Krel_cells, *Krel_faces);
    matrix->CreateMFDrhsVectors();
    AddGravityFluxes_MFD(K, *Krel_cells, *Krel_faces, matrix);
    AddTimeDerivative_MFDpicard(*solution_cells, *solution_old_cells, dTp, matrix);
    matrix->ApplyBoundaryConditions(bc_markers, bc_values);
    matrix->AssembleGlobalMatrices();
    rhs = matrix->rhs();

    // create preconditioner
    preconditioner->CreateMFDstiffnessMatrices(*Krel_cells, *Krel_faces);
    preconditioner->CreateMFDrhsVectors();
    AddGravityFluxes_MFD(K, *Krel_cells, *Krel_faces, preconditioner);
    AddTimeDerivative_MFD(*solution_old_cells, dTp, preconditioner);
    preconditioner->ApplyBoundaryConditions(bc_markers, bc_values);
    preconditioner->AssembleGlobalMatrices();
    preconditioner->ComputeSchurComplement(bc_markers, bc_values);
    preconditioner->UpdateML_Preconditioner();

    // call AztecOO solver
    solver->SetRHS(&*rhs);  // AztecOO modifies the right-hand-side.
    solution_new = solution_old;  // initial solution guess
    solver->SetLHS(&solution_new);

    solver->Iterate(max_itrs, convergence_tol);
    int num_itrs = solver->NumIters();
    double linear_residual = solver->TrueResidual();

    // error analysis
    Epetra_Vector solution_diff(*solution_old_cells);
    solution_diff.Update(-1.0, *solution_new_cells, 1.0);
    double error = ErrorNormSTOMP(*solution_old_cells, solution_diff);

    if (MyPID == 0 && verbosity >= FLOW_VERBOSITY_HIGH) {
      std::printf("Picard:%4d   Pressure(error=%9.4e)  solver(%8.3e, %4d)\n",
          itrs, error, linear_residual, num_itrs);
    }

    if (error < 1e-4) 
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

  Teuchos::RCP<Epetra_Vector> solution_old_cells = Teuchos::rcp(FS->CreateCellView(solution_old));
  Teuchos::RCP<Epetra_Vector> solution_new_cells = Teuchos::rcp(FS->CreateCellView(solution_new));

  if (! is_matrix_symmetric) solver->SetAztecOption(AZ_solver, AZ_gmres);
  solver->SetAztecOption(AZ_output, AZ_none);

  T_internal = T0_sss;
  dT = dT0_sss;

  int itrs = 0, ifail = 0;
  double L2error = 1.0;
  while (L2error > convergence_tol_sss && itrs < max_itrs_sss) {
    if (!is_matrix_symmetric) {  // Define K and Krel_faces
      CalculateRelativePermeabilityFace(*solution_cells);
      Krel_cells->PutScalar(1.0);
    } else {  // Define K and Krel_cells, Krel_faces is always one
      CalculateRelativePermeabilityCell(*solution_cells);
      Krel_faces->PutScalar(1.0);
    }

    // update boundary conditions
    double time = T_internal + dT;
    bc_pressure->Compute(time);
    bc_flux->Compute(time);
    bc_head->Compute(time);
    bc_seepage->Compute(time);
    UpdateBoundaryConditions(
        bc_pressure, bc_head, bc_flux, bc_seepage,
        *solution_faces, atm_pressure,
        bc_markers, bc_values);

    // create algebraic problem
    matrix->CreateMFDstiffnessMatrices(*Krel_cells, *Krel_faces);
    matrix->CreateMFDrhsVectors();
    AddGravityFluxes_MFD(K, *Krel_cells, *Krel_faces, matrix);
    AddTimeDerivative_MFDfake(*solution_cells, dT, matrix);
    matrix->ApplyBoundaryConditions(bc_markers, bc_values);
    matrix->AssembleGlobalMatrices();

    // create preconditioner
    preconditioner->CreateMFDstiffnessMatrices(*Krel_cells, *Krel_faces);
    preconditioner->CreateMFDrhsVectors();
    AddGravityFluxes_MFD(K, *Krel_cells, *Krel_faces, preconditioner);
    AddTimeDerivative_MFDfake(*solution_cells, dT, preconditioner);
    preconditioner->ApplyBoundaryConditions(bc_markers, bc_values);
    preconditioner->AssembleGlobalMatrices();
    preconditioner->ComputeSchurComplement(bc_markers, bc_values);
    preconditioner->UpdateML_Preconditioner();

    // call AztecOO solver
    rhs = matrix->rhs();
    Epetra_Vector b(*rhs);
    solver->SetRHS(&b);  // AztecOO modifies the right-hand-side.
    solver->SetLHS(&*solution);  // initial solution guess

    solver->Iterate(max_itrs, convergence_tol);
    int num_itrs = solver->NumIters();
    double residual = solver->TrueResidual();

    // error estimates
    double sol_norm = FS->normLpCell(solution_new, 2.0);
    Epetra_Vector solution_diff(*solution_old_cells);
    solution_diff.Update(-1.0, *solution_new_cells, 1.0);
    L2error = ErrorNormPicardExperimental(*solution_old_cells, solution_diff);

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
* Adds time derivative to cell-based part of MFD algebraic system.                                               
****************************************************************** */
void Richards_PK::AddTimeDerivative_MFDfake(
    Epetra_Vector& pressure_cells, double dT_prec, Matrix_MFD* matrix)
{
  std::vector<double>& Acc_cells = matrix->Acc_cells();
  std::vector<double>& Fc_cells = matrix->Fc_cells();

  for (int c = 0; c < ncells_owned; c++) {
    double volume = mesh_->cell_volume(c);
    double factor = volume / dT_prec;
    Acc_cells[c] += factor;
    Fc_cells[c] += factor * pressure_cells[c];
  }
}


/* ******************************************************************
* Check difference du between the predicted and converged solutions.                                            
****************************************************************** */
double Richards_PK::ErrorNormPicardExperimental(const Epetra_Vector& u, const Epetra_Vector& du)
{
  double error_p = 0.0;
  for (int c = 0; c < ncells_owned; c++) {
    double tmp = fabs(du[c]) / (fabs(u[c] - atm_pressure) + atm_pressure);
    error_p = std::max<double>(error_p, tmp);
  }

#ifdef HAVE_MPI
  double buf = error_p;
  u.Comm().MaxAll(&buf, &error_p, 1);  // find the global maximum
#endif
  return error_p;
}

}  // namespace AmanziFlow
}  // namespace Amanzi

