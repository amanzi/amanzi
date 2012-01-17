/*
This is the flow component of the Amanzi code. 
License: BSD
Author:  Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include "Richards_PK.hpp"

namespace Amanzi {
namespace AmanziFlow {

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

  T_internal = T0_sss;
  dT = dT0_sss;

  int itrs = 0, ifail = 0;
  double L2error = 1.0;
  while (L2error > convergence_tol_sss && itrs < max_itrs_sss) {
    calculateRelativePermeability(*solution_cells);
    setAbsolutePermeabilityTensor(K);
    if (flag_upwind) {  // Define K and Krel_faces
      calculateRelativePermeabilityUpwindGravity(*solution_cells);
      for (int c=0; c<K.size(); c++) K[c] *= rho / mu;
    } else {  // Define K and Krel_cells, Krel_faces is always one
      for (int c=0; c<K.size(); c++) K[c] *= (*Krel_cells)[c] * rho / mu;  
    }

    // update boundary conditions
    double time = T_internal + dT;
    bc_pressure->Compute(time);
    bc_flux->Compute(time);
    bc_head->Compute(time);
    updateBoundaryConditions(bc_pressure, bc_head, bc_flux, bc_markers, bc_values);

    // create algebraic problem (matrix = preconditioner)
    matrix->createMFDstiffnessMatrices(K, *Krel_faces);
    matrix->createMFDrhsVectors();
    addGravityFluxes_MFD(K, *Krel_faces, matrix);
    addTimeDerivative_MFDfake(*solution_cells, dT, matrix);
    matrix->applyBoundaryConditions(bc_markers, bc_values);
    matrix->assembleGlobalMatrices();
    matrix->computeSchurComplement(bc_markers, bc_values);
    matrix->update_ML_preconditioner();

    // call AztecOO solver
    rhs = matrix->get_rhs();
    Epetra_Vector b(*rhs);
    solver->SetRHS(&b);  // AztecOO modifies the right-hand-side.
    solver->SetLHS(&*solution);  // initial solution guess 

    solver->Iterate(max_itrs, convergence_tol);
    int num_itrs = solver->NumIters();
    double residual = solver->TrueResidual();

    // error estimates
    double sol_norm = FS->normL2cell(solution_new);
    L2error = errorSolutionDiff(solution_old, solution_new);

    if (L2error > 10.0 && itrs && ifail < 5) {  // itrs=0 allows to avoid bad initial guess. 
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
  while (L2error > convergence_tol_sss) {
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
    if (MyPID == 0 && verbosity >= FLOW_VERBOSITY_HIGH) {
      std::printf("Time step:%6d   Pressure(diff=%9.4e, sol=%9.4e, res=%9.4e)  time=%8.3e\n", 
          itrs, sol_error, sol_norm, L2error, T_internal);
    }

    T_internal += dT;
    itrs++;
  }
  return 0;
}


/* ******************************************************************
* Adds time derivative to cell-based part of MFD algebraic system.                                               
****************************************************************** */
void Richards_PK::addTimeDerivative_MFDfake(
   Epetra_Vector& pressure_cells, double dT_prec, Matrix_MFD* matrix)
{
  std::vector<double>& Acc_cells = matrix->get_Acc_cells();
  std::vector<double>& Fc_cells = matrix->get_Fc_cells();

  int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);

  for (int c=0; c<ncells; c++) {
    double volume = mesh_->cell_volume(c);
    double factor = volume / dT_prec;
    Acc_cells[c] += factor;
    Fc_cells[c] += factor * pressure_cells[c];
  }
}


/* ******************************************************************
* Check difference between solutions at times T and T+dT.                                                 
****************************************************************** */
double Richards_PK::errorSolutionDiff(const Epetra_Vector& uold, const Epetra_Vector& unew)
{
  double error_norm = 0.0; 
  for (int n=0; n<uold.MyLength(); n++) {
    double tmp = abs(uold[n] - unew[n]) / (absolute_tol_sss + relative_tol_sss * abs(uold[n]));
    error_norm = std::max<double>(error_norm, tmp);
  }

  // find the global maximum
#ifdef HAVE_MPI
  double buf = error_norm;
  MPI_Allreduce(&buf, &error_norm, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#endif
  return  error_norm;
}

}  // namespace AmanziFlow
}  // namespace Amanzi

