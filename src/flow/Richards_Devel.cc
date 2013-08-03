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

#include "Richards_PK.hh"
#include "TI_Specs.hh"


namespace Amanzi {
namespace AmanziFlow {

/* ******************************************************************
* Calculates steady-state solution assuming that abosolute and relative
* permeabilities do not depend explicitly on time.
* WARNING: temporary replacement for missing BDF1 time integrator.                                                    
****************************************************************** */
int Richards_PK::AdvanceToSteadyState_BackwardEuler(TI_Specs& ti_specs)
{
  Epetra_Vector  solution_old(*solution);
  Epetra_Vector& solution_new = *solution;

  Teuchos::RCP<Epetra_Vector> solution_old_cells = Teuchos::rcp(FS->CreateCellView(solution_old));
  Teuchos::RCP<Epetra_Vector> solution_new_cells = Teuchos::rcp(FS->CreateCellView(solution_new));

  if (is_matrix_symmetric) 
      solver->SetAztecOption(AZ_solver, AZ_cg);
  else 
      solver->SetAztecOption(AZ_solver, AZ_gmres);
  solver->SetAztecOption(AZ_output, verbosity_AztecOO);

  int max_itrs_nonlinear = ti_specs_sss_.max_itrs;
  double T1 = ti_specs_sss_.T1;
  double dTmax = ti_specs_sss_.dTmax;
  double residual_tol_nonlinear = ti_specs_sss_.residual_tol;

  int itrs = 0, ifail = 0;
  double L2error = 1.0;
  while (L2error > residual_tol_nonlinear && itrs < max_itrs_nonlinear) {
    // update permeabilities
    CalculateRelativePermeability(*solution);
 
    // update boundary conditions
    double time = T_physics + dT;
    bc_pressure->Compute(time);
    bc_flux->Compute(time);
    bc_seepage->Compute(time);
    if (shift_water_table_.getRawPtr() == NULL)
      bc_head->Compute(time);
    else
      bc_head->ComputeShift(time, shift_water_table_->Values());

    ProcessBoundaryConditions(
        bc_pressure, bc_head, bc_flux, bc_seepage,
        *solution_cells, *solution_faces, atm_pressure,
        rainfall_factor, bc_submodel, bc_model, bc_values);

    // create algebraic problem
    matrix_->CreateMFDstiffnessMatrices(*Krel_cells, *Krel_faces, Krel_method);
    matrix_->CreateMFDrhsVectors();
    AddGravityFluxes_MFD(K, *Krel_cells, *Krel_faces, Krel_method, matrix_);
    AddTimeDerivative_MFDfake(*solution_cells, dT, matrix_);
    matrix_->ApplyBoundaryConditions(bc_model, bc_values);
    matrix_->AssembleGlobalMatrices();

    // create preconditioner
    preconditioner_->CreateMFDstiffnessMatrices(*Krel_cells, *Krel_faces, Krel_method);
    preconditioner_->CreateMFDrhsVectors();
    AddTimeDerivative_MFDfake(*solution_cells, dT, preconditioner_);
    preconditioner_->ApplyBoundaryConditions(bc_model, bc_values);
    preconditioner_->AssembleSchurComplement(bc_model, bc_values);
    preconditioner_->UpdatePreconditioner();

    // call AztecOO solver
    rhs = matrix_->rhs();
    Epetra_Vector b(*rhs);
    solver->SetRHS(&b);  // AztecOO modifies the right-hand-side.
    solver->SetLHS(&*solution);  // initial solution guess

    solver->Iterate(max_itrs_linear, convergence_tol_linear);
    int num_itrs = solver->NumIters();

    // error estimates
    double sol_norm = FS->normLpCell(solution_new, 2.0);
    Epetra_Vector solution_diff(*solution_old_cells);
    solution_diff.Update(-1.0, *solution_new_cells, 1.0);
    L2error = ErrorNormPicardExperimental(*solution_old_cells, solution_diff);

    if (L2error > 1.0 && itrs && ifail < 5) {  // itrs=0 allows to avoid bad initial guess.
      dT /= 10;
      solution_new = solution_old;
      if (MyPID == 0 && verbosity >= FLOW_VERBOSITY_HIGH) {
        double linear_residual = solver->ScaledResidual();
        std::printf("Fail:%4d  Pressure(diff=%9.4e, sol=%9.4e)  solver(%8.3e,%3d), T=%9.3e dT=%7.2e\n",
            itrs, L2error, sol_norm, linear_residual, num_itrs, T_physics, dT);
      }
      ifail++;
    } else {
      T_physics += dT;
      solution_old = solution_new;

      if (MyPID == 0 && verbosity >= FLOW_VERBOSITY_HIGH) {
        double linear_residual = solver->ScaledResidual();
        std::printf("Step:%4d  Pressure(diff=%9.4e, sol=%9.4e)  solver(%8.3e,%3d), T=%9.3e dT=%7.2e\n",
            itrs, L2error, sol_norm, linear_residual, num_itrs, T_physics, dT);
      }

      ifail = 0;
      dT = std::min(dT*1.25, dTmax);
      itrs++;
    }

    if (T_physics > T1) break;
  }

  ti_specs_sss_.num_itrs = itrs;
  return 0;
}


/* ******************************************************************
* Adds time derivative to cell-based part of MFD algebraic system.                                               
****************************************************************** */
void Richards_PK::AddTimeDerivative_MFDfake(
    Epetra_Vector& pressure_cells, double dT_prec, Matrix_MFD* matrix_operator)
{
  std::vector<double>& Acc_cells = matrix_operator->Acc_cells();
  std::vector<double>& Fc_cells = matrix_operator->Fc_cells();

  for (int c = 0; c < ncells_owned; c++) {
    double volume = mesh_->cell_volume(c);
    double factor = volume / dT_prec;
    Acc_cells[c] += factor;
    Fc_cells[c] += factor * pressure_cells[c];
  }
}


}  // namespace AmanziFlow
}  // namespace Amanzi

