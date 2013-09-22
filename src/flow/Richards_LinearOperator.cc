/*
This is the flow component of the Amanzi code. 

Copyright 2010-2012 held jointly by LANS/LANL, LBNL, and PNNL. 
Amanzi is released under the three-clause BSD License. 
The terms of use and "as is" disclaimer for this license are 
provided Reconstruction.cppin the top-level COPYRIGHT file.

Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include "Matrix_MFD.hh"
#include "Matrix_Audit.hh"
#include "LinearOperatorFactory.hh"
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
  const Epetra_Vector& rhs = *(matrix_->rhs());

  AssembleSteadyStatePreconditioner_MFD(&*preconditioner_);
  preconditioner_->UpdatePreconditioner();

  // solve linear problem
  AmanziSolvers::LinearOperatorFactory<Matrix_MFD, Epetra_Vector, Epetra_Map> factory;

  Teuchos::RCP<AmanziSolvers::LinearOperator<Matrix_MFD, Epetra_Vector, Epetra_Map> >
     solver = factory.Create(ls_specs.solver_name, solver_list_, matrix_, preconditioner_);

  solver->ApplyInverse(rhs, u);

  if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH) {
    int num_itrs = solver->num_itrs();
    double residual = solver->residual();

    Teuchos::OSTab tab = vo_->getOSTab();
    *(vo_->os()) << "saturated solver(" << solver->name() 
                 << "): ||r||=" << residual << " itr=" << num_itrs << endl;
  }
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

  // solve non-symmetric problem
  LinearSolver_Specs& ls_specs = ti_specs->ls_specs_constraints;

  AmanziSolvers::LinearOperatorFactory<Matrix_MFD, Epetra_Vector, Epetra_Map> factory;
  Teuchos::RCP<AmanziSolvers::LinearOperator<Matrix_MFD, Epetra_Vector, Epetra_Map> >
     solver = factory.Create(ls_specs.solver_name, solver_list_, matrix_, preconditioner_);

  Epetra_Vector& rhs = *(matrix_->rhs());
  solver->ApplyInverse(rhs, utmp);

  *u_faces = *utmp_faces;

  if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH) {
    int num_itrs = solver->num_itrs();
    double residual = solver->residual();

    Teuchos::OSTab tab = vo_->getOSTab();
    *(vo_->os()) << "constraints solver(" << solver->name() 
                 << "): ||r||=" << residual << " itr=" << num_itrs << endl;
  }
}

}  // namespace AmanziFlow
}  // namespace Amanzi

