/*
This is the flow component of the Amanzi code. 

Copyright 2010-2012 held jointly by LANS/LANL, LBNL, and PNNL. 
Amanzi is released under the three-clause BSD License. 
The terms of use and "as is" disclaimer for this license are 
provided in the top-level COPYRIGHT file.

Author:  Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include "Richards_PK.hh"


namespace Amanzi {
namespace AmanziFlow {

/* ******************************************************************
* Makes one Picard step during transient time integration.
* This is the experimental method.                                                 
****************************************************************** */
int Richards_PK::PicardTimeStep(double Tp, double dTp, double& dTnext)
{
  // p^{k-1} = solution_old, p^k = solution_new
  Epetra_Vector solution_old(*solution);
  Epetra_Vector solution_new(*solution);

  Epetra_Vector* solution_old_cells = FS->CreateCellView(solution_old);
  Epetra_Vector* solution_old_faces = FS->CreateFaceView(solution_old);
  Epetra_Vector* solution_new_cells = FS->CreateCellView(solution_new);

  if (is_matrix_symmetric) 
      solver->SetAztecOption(AZ_solver, AZ_cg);
  else 
      solver->SetAztecOption(AZ_solver, AZ_gmres);
  solver->SetAztecOption(AZ_output, AZ_none);
  solver->SetAztecOption(AZ_conv, AZ_rhs);

  double error0;
  int itrs;
  for (itrs = 0; itrs < 20; itrs++) {
    CalculateRelativePermeability(solution_old);

    double time = Tp + dTp;
    UpdateSourceBoundaryData(time, *solution_old_faces);

    // create algebraic problem
    matrix_->CreateMFDstiffnessMatrices(*Krel_cells, *Krel_faces, Krel_method);
    matrix_->CreateMFDrhsVectors();
    AddGravityFluxes_MFD(K, *Krel_cells, *Krel_faces, Krel_method, matrix_);
    AddTimeDerivative_MFDpicard(*solution, *solution_cells, dTp, matrix_);
    matrix_->ApplyBoundaryConditions(bc_model, bc_values);
    matrix_->AssembleGlobalMatrices();
    rhs = matrix_->rhs();

    // create preconditioner
    preconditioner_->CreateMFDstiffnessMatrices(*Krel_cells, *Krel_faces, Krel_method);
    preconditioner_->CreateMFDrhsVectors();
    AddGravityFluxes_MFD(K, *Krel_cells, *Krel_faces, Krel_method, preconditioner_);
    AddTimeDerivative_MFDpicard(*solution, *solution_cells, dTp, preconditioner_);
    preconditioner_->ApplyBoundaryConditions(bc_model, bc_values);
    preconditioner_->AssembleGlobalMatrices();
    preconditioner_->AssembleSchurComplement(bc_model, bc_values);
    preconditioner_->UpdatePreconditioner();

    // call AztecOO solver
    solver->SetRHS(&*rhs);
    solution_new = solution_old;  // initial solution guess
    solver->SetLHS(&solution_new);

    solver->Iterate(max_itrs_linear, convergence_tol_linear);
    int num_itrs = solver->NumIters();

    // error analysis
    Epetra_Vector solution_diff(*solution_old_cells);
    solution_diff.Update(-1.0, *solution_new_cells, 1.0);
    double error = ErrorNormPicardExperimental(*solution_old_cells, solution_diff);

    if (MyPID == 0 && verbosity >= FLOW_VERBOSITY_HIGH) {
      double linear_residual = solver->ScaledResidual();
      std::printf("Picard:%4d   Pressure(error=%9.4e)  solver(%8.3e, %4d)\n",
          itrs, error, linear_residual, num_itrs);
    }

    if (error < 1e-5 && itrs > 0) 
      break;
    else 
      solution_old = solution_new;
  }

  if (itrs < 10) {
    *solution = solution_new;
    dTnext = dT * 1.2;
  } else if (itrs < 15) {
    *solution = solution_new;
    dTnext = dT;
  } else if (itrs < 19) {
    *solution = solution_new;
    dTnext = dT * 0.5;
    throw itrs;
  }

  return itrs;
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

