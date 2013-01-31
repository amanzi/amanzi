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
  matrix_->CreateMFDstiffnessMatrices(*Krel_cells, *Krel_faces, FLOW_RELATIVE_PERM_NONE);
  matrix_->CreateMFDrhsVectors();
  AddGravityFluxes_MFD(K, *Krel_cells, *Krel_faces, FLOW_RELATIVE_PERM_NONE, matrix_);
  matrix_->ApplyBoundaryConditions(bc_model, bc_values);
  matrix_->AssembleGlobalMatrices();
  matrix_->ComputeSchurComplement(bc_model, bc_values);
  matrix_->UpdatePreconditioner();
}


/* ******************************************************************
* Calculates steady-state solution assuming that absolute permeability 
* does not depend on time. The boundary conditions are calculated
* only once, during the initialization step.                                                
****************************************************************** */
void Darcy_PK::SolveFullySaturatedProblem(double Tp, Epetra_Vector& u)
{
  solver->SetAztecOption(AZ_output, verbosity_AztecOO);
  solver->SetAztecOption(AZ_conv, AZ_rhs);

  // calculate and assemble elemental stifness matrices
  matrix_->CreateMFDstiffnessMatrices(*Krel_cells, *Krel_faces, FLOW_RELATIVE_PERM_NONE);
  matrix_->CreateMFDrhsVectors();
  AddGravityFluxes_MFD(K, *Krel_cells, *Krel_faces, FLOW_RELATIVE_PERM_NONE, matrix_);
  matrix_->ApplyBoundaryConditions(bc_model, bc_values);
  matrix_->AssembleGlobalMatrices();
  matrix_->ComputeSchurComplement(bc_model, bc_values);
  matrix_->UpdatePreconditioner();

  rhs = matrix_->rhs();
  Epetra_Vector b(*rhs);
  solver->SetRHS(&b);  // Aztec00 modifies the right-hand-side.
  solver->SetLHS(&u);  // initial solution guess

  int max_itrs = ti_specs_sss.ls_specs.max_itrs;
  double convergence_tol = ti_specs_sss.ls_specs.convergence_tol;

  solver->Iterate(max_itrs, convergence_tol);

  int num_itrs = solver->NumIters();
  double linear_residual = solver->ScaledResidual();

  if (verbosity >= FLOW_VERBOSITY_HIGH && MyPID == 0) {
    std::printf("Flow PK: pressure solver: ||r||=%8.3e itr=%d\n", linear_residual, num_itrs);
  }
}


/* ******************************************************************
* Calculates steady-state solution using a user-given rhs vector. 
* The matrix has to be assembled before this call.
****************************************************************** */
void Darcy_PK::SolveFullySaturatedProblem(double Tp, const Epetra_Vector& rhs, Epetra_Vector& u)
{
  solver->SetAztecOption(AZ_output, verbosity_AztecOO);
  solver->SetAztecOption(AZ_conv, AZ_rhs);

  Epetra_Vector b(rhs);
  solver->SetRHS(&b);  // Aztec00 modifies the right-hand-side.
  solver->SetLHS(&u);  // initial solution guess

  int max_itrs = ti_specs_sss.ls_specs.max_itrs;
  double convergence_tol = ti_specs_sss.ls_specs.convergence_tol;

  solver->Iterate(max_itrs, convergence_tol);

  int num_itrs = solver->NumIters();
  double linear_residual = solver->ScaledResidual();

  if (verbosity >= FLOW_VERBOSITY_HIGH && MyPID == 0) {
    std::printf("Flow PK: pressure solver: ||r||=%8.3e itr=%d\n", linear_residual, num_itrs);
  }
}

}  // namespace AmanziFlow
}  // namespace Amanzi


