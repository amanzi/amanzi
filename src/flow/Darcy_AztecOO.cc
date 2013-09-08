/*
This is the flow component of the Amanzi code.

Copyright 2010-2012 held jointly by LANS/LANL, LBNL, and PNNL. 
Amanzi is released under the three-clause BSD License. 
The terms of use and "as is" disclaimer for this license are 
provided in the top-level COPYRIGHT file.

Authors: Konstantin Lipnikov (version 2) (lipnikov@lanl.gov)
*/

#include <vector>

#include "Darcy_PK.hh"

namespace Amanzi {
namespace AmanziFlow {

/* ******************************************************************
* Gathers together routines to compute steady-state MFD matrices.                            
****************************************************************** */
void Darcy_PK::AssembleMatrixMFD()
{
  matrix_->CreateMFDstiffnessMatrices();
  matrix_->CreateMFDrhsVectors();
  AddGravityFluxes_MFD(K, matrix_);
  matrix_->ApplyBoundaryConditions(bc_model, bc_values);
  matrix_->AssembleGlobalMatrices();
  matrix_->AssembleSchurComplement(bc_model, bc_values);
  matrix_->UpdatePreconditioner();
}


/* ******************************************************************
* Calculates steady-state solution assuming that absolute permeability 
* does not depend on time. The boundary conditions are calculated
* only once, during the initialization step.                                                
****************************************************************** */
void Darcy_PK::SolveFullySaturatedProblem(double Tp, Epetra_Vector& u)
{
  LinearSolver_Specs& ls_specs = ti_specs->ls_specs;

  AztecOO* solver = new AztecOO;
  solver->SetUserOperator(matrix_);
  solver->SetPrecOperator(preconditioner_);

  solver->SetAztecOption(AZ_solver, ls_specs.method);  // Must be an AZ_xxx method.
  solver->SetAztecOption(AZ_output, verbosity_AztecOO);
  solver->SetAztecOption(AZ_conv, AZ_rhs);

  // calculate and assemble elemental stifness matrices
  matrix_->CreateMFDstiffnessMatrices();
  matrix_->CreateMFDrhsVectors();
  AddGravityFluxes_MFD(K, matrix_);
  matrix_->ApplyBoundaryConditions(bc_model, bc_values);
  matrix_->AssembleGlobalMatrices();
  matrix_->AssembleSchurComplement(bc_model, bc_values);
  matrix_->UpdatePreconditioner();

  rhs = matrix_->rhs();
  Epetra_Vector b(*rhs);
  solver->SetRHS(&b);  // Aztec00 modifies the right-hand-side.
  solver->SetLHS(&u);  // initial solution guess

  solver->Iterate(ls_specs.max_itrs, ls_specs.convergence_tol);

  int num_itrs = solver->NumIters();
  double linear_residual = solver->ScaledResidual();

  if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH) {
    Teuchos::OSTab tab = vo_->getOSTab();
    *(vo_->os()) << "pressure solver: ||r||=" << linear_residual << " itr=" << num_itrs << endl;
  }

  delete solver;
}


/* ******************************************************************
* Calculates steady-state solution using a user-given rhs vector. 
* The matrix has to be assembled before this call.
****************************************************************** */
void Darcy_PK::SolveFullySaturatedProblem(double Tp, const Epetra_Vector& rhs, Epetra_Vector& u)
{
  LinearSolver_Specs& ls_specs = ti_specs->ls_specs;

  AztecOO* solver = new AztecOO;
  solver->SetUserOperator(matrix_);
  solver->SetPrecOperator(preconditioner_);

  solver->SetAztecOption(AZ_solver, ls_specs.method);  // Must be an AZ_xxx method.
  solver->SetAztecOption(AZ_output, verbosity_AztecOO);
  solver->SetAztecOption(AZ_conv, AZ_rhs);

  Epetra_Vector b(rhs);
  solver->SetRHS(&b);  // Aztec00 modifies the right-hand-side.
  solver->SetLHS(&u);  // initial solution guess

  solver->Iterate(ls_specs.max_itrs, ls_specs.convergence_tol);

  if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH) {
    int num_itrs = solver->NumIters();
    double linear_residual = solver->ScaledResidual();

    Teuchos::OSTab tab = vo_->getOSTab();
    *(vo_->os()) << "pressure solver: ||r||=" << linear_residual << " itr=" << num_itrs << endl;
  }

  delete solver;
}

}  // namespace AmanziFlow
}  // namespace Amanzi


