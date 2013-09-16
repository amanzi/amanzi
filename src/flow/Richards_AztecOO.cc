/*
This is the flow component of the Amanzi code. 
Solvers based on AztecOO are collected here.

Copyright 2010-2012 held jointly by LANS/LANL, LBNL, and PNNL. 
Amanzi is released under the three-clause BSD License. 
The terms of use and "as is" disclaimer for this license are 
provided Reconstruction.cppin the top-level COPYRIGHT file.

Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include "Matrix_MFD.hh"
#include "Matrix_Audit.hh"
#include "Richards_PK.hh"

namespace Amanzi {
namespace AmanziFlow {

/* ******************************************************************
* Calculate saturated pressure solution using boundary conditions 
* at time Tp.
* WARNING: data in vectors Krel and rhs are destroyed.
****************************************************************** */
void Richards_PK::SolveFullySaturatedProblem(double Tp, Epetra_Vector& u, LinearSolver_Specs& ls_specs)
{
  Epetra_Vector* u_cells = FS->CreateCellView(u);
  Epetra_Vector* u_faces = FS->CreateFaceView(u);
  UpdateSourceBoundaryData(Tp, *u_cells, *u_faces);

  // set fully saturated media
  rel_perm->SetFullySaturated();

  // calculate and assemble elemental stiffness matrices
  AssembleSteadyStateMatrix_MFD(&*matrix_);
  AssembleSteadyStatePreconditioner_MFD(&*preconditioner_);
  preconditioner_->UpdatePreconditioner();

  // solve symmetric problem
  AztecOO* solver_tmp = new AztecOO;

  solver_tmp->SetUserOperator(&*matrix_);
  solver_tmp->SetPrecOperator(&*preconditioner_);
  solver_tmp->SetAztecOption(AZ_solver, ls_specs.method);  // Must be AZ_xxx method.
  solver_tmp->SetAztecOption(AZ_output, verbosity_AztecOO);
  solver_tmp->SetAztecOption(AZ_conv, AZ_rhs);

  Epetra_Vector b(*(matrix_->rhs()));
  solver_tmp->SetRHS(&b);
  solver_tmp->SetLHS(&u);

  solver_tmp->Iterate(ls_specs.max_itrs, ls_specs.convergence_tol);

  // Matrix_Audit audit(mesh_, matrix_);
  // audit.InitAudit();
  // audit.RunAudit();

  if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH) {
    int num_itrs = solver_tmp->NumIters();
    double linear_residual = solver_tmp->ScaledResidual();

    Teuchos::OSTab tab = vo_->getOSTab();
    *(vo_->os()) << "saturated solver: ||r||=" << linear_residual 
                 << " itr=" << num_itrs << endl;
  }

  delete solver_tmp;
}


/* ******************************************************************
* Enforce constraints at time Tp by solving diagonalized MFD problem.
* Algorithm is based on de-coupling pressure-lambda system.
****************************************************************** */
void Richards_PK::EnforceConstraints_MFD(double Tp, Epetra_Vector& u)
{
  Epetra_Vector utmp(u);
  Epetra_Vector* u_cells = FS->CreateCellView(u);
  Epetra_Vector* u_faces = FS->CreateFaceView(u);
  Epetra_Vector* utmp_faces = FS->CreateFaceView(utmp);

  UpdateSourceBoundaryData(Tp, *u_cells, *u_faces);

  // calculate and assemble elemental stiffness matrix
  rel_perm->Compute(u, bc_model, bc_values);
  AssembleSteadyStateMatrix_MFD(&*matrix_);
  matrix_->ReduceGlobalSystem2LambdaSystem(u);

  // copy stiffness matrix to preconditioner (raw-data)
  preconditioner_->PopulatePreconditioner(*matrix_);
  preconditioner_->UpdatePreconditioner();

  LinearSolver_Specs& ls_specs = ti_specs->ls_specs_constraints;

  // solve non-symmetric problem
  AztecOO* solver_tmp = new AztecOO;

  solver_tmp->SetUserOperator(&*matrix_);
  solver_tmp->SetPrecOperator(&*preconditioner_);
  solver_tmp->SetAztecOption(AZ_solver, ls_specs.method);  // Must be AZ_xxx method. 
  solver_tmp->SetAztecOption(AZ_output, verbosity_AztecOO);
  solver_tmp->SetAztecOption(AZ_conv, AZ_rhs);

  Epetra_Vector b(*(matrix_->rhs()));
  solver_tmp->SetRHS(&b);

  solver_tmp->SetLHS(&utmp);
  solver_tmp->Iterate(ls_specs.max_itrs, ls_specs.convergence_tol);
  *u_faces = *utmp_faces;

  if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH) {
    int num_itrs = solver_tmp->NumIters();
    double linear_residual = solver_tmp->ScaledResidual();

    Teuchos::OSTab tab = vo_->getOSTab();
    *(vo_->os()) << "constraints solver: ||r||=" << linear_residual
                 << " itr=" << num_itrs << endl;
  }

  delete solver_tmp;
}

}  // namespace AmanziFlow
}  // namespace Amanzi

