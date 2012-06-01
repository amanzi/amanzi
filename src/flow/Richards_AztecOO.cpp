/*
This is the flow component of the Amanzi code. 
Solvers based on AztecOO are collected here.

Copyright 2010-2012 held jointly by LANS/LANL, LBNL, and PNNL. 
Amanzi is released under the three-clause BSD License. 
The terms of use and "as is" disclaimer for this license are 
provided Reconstruction.cppin the top-level COPYRIGHT file.

Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include "Matrix_MFD.hpp"
#include "Richards_PK.hpp"

namespace Amanzi {
namespace AmanziFlow {

/* ******************************************************************
* Calculate saturated pressure solution using boundary conditions 
* at time Tp.
* WARNING: data in vectors Krel and rhs are destroyed.
****************************************************************** */
void Richards_PK::SolveFullySaturatedProblem(double Tp, Epetra_Vector& u)
{
  Epetra_Vector* u_faces = FS->CreateFaceView(u);

  // update boundary conditions
  bc_pressure->Compute(Tp);
  bc_head->Compute(Tp);
  bc_flux->Compute(Tp);
  bc_seepage->Compute(Tp);
  ProcessBoundaryConditions(
      bc_pressure, bc_head, bc_flux, bc_seepage,
      *u_faces, atm_pressure,
      bc_markers, bc_values);

  // set fully saturated media
  Krel_cells->PutScalar(1.0);
  Krel_faces->PutScalar(1.0);

  // calculate and assemble elemental stiffness matrices
  AssembleSteadyStateProblem_MFD(matrix, false);
  AssembleSteadyStateProblem_MFD(preconditioner, true);
  preconditioner->UpdateML_Preconditioner();

  // solve symmetric problem
  AztecOO* solver_tmp = new AztecOO;

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
* Calculate transient pressure solution using boundary conditions 
* at time Tp.
* WARNING: solution u is the input/output parameter.
****************************************************************** */
void Richards_PK::SolveTransientProblem(double Tp, double dTp, Epetra_Vector& u)
{
  Epetra_Vector* u_faces = FS->CreateFaceView(u);

  // update boundary conditions
  bc_pressure->Compute(Tp);
  bc_head->Compute(Tp);
  bc_flux->Compute(Tp);
  bc_seepage->Compute(Tp);
  ProcessBoundaryConditions(
      bc_pressure, bc_head, bc_flux, bc_seepage,
      *u_faces, atm_pressure,
      bc_markers, bc_values);

  // calculate and assemble elemental stiffness matrices
  AssembleTransientProblem_MFD(matrix, dTp, u, false);
  AssembleTransientProblem_MFD(preconditioner, dTp, u, true);
  preconditioner->UpdateML_Preconditioner();

  // solve symmetric problem
  AztecOO* solver_tmp = new AztecOO;

  solver_tmp->SetUserOperator(matrix);
  solver_tmp->SetPrecOperator(preconditioner);
  solver_tmp->SetAztecOption(AZ_solver, AZ_gmres);
  solver_tmp->SetAztecOption(AZ_output, AZ_none);
  solver_tmp->SetAztecOption(AZ_conv, AZ_rhs);

  Epetra_Vector b(*(matrix->rhs()));
  solver_tmp->SetRHS(&b);

  solver_tmp->SetLHS(&u);
  solver_tmp->Iterate(max_itrs, convergence_tol);

  if (MyPID == 0 && verbosity >= FLOW_VERBOSITY_HIGH) {
    int num_itrs = solver_tmp->NumIters();
    double linear_residual = solver_tmp->TrueResidual();
    std::printf("Richards PK: transient pressure solver(%8.3e, %4d)\n", linear_residual, num_itrs);
  }

  delete solver_tmp;
}

}  // namespace AmanziFlow
}  // namespace Amanzi

