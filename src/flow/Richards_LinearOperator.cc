/*
This is the flow component of the Amanzi code. 

Copyright 2010-2012 held jointly by LANS/LANL, LBNL, and PNNL. 
Amanzi is released under the three-clause BSD License. 
The terms of use and "as is" disclaimer for this license are 
provided Reconstruction.cppin the top-level COPYRIGHT file.

Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include "Matrix_MFD.hh"
#include "LinearOperatorFactory.hh"
#include "Richards_PK.hh"

namespace Amanzi {
namespace AmanziFlow {

/* ******************************************************************
* Calculate saturated pressure solution using boundary conditions 
* at time Tp.
* WARNING: data in vectors Krel and rhs are destroyed.
****************************************************************** */
void Richards_PK::SolveFullySaturatedProblem(
    double Tp, CompositeVector& u, LinearSolver_Specs& ls_specs)
{
  UpdateSourceBoundaryData(Tp, u);

  // set fully saturated media
  rel_perm->SetFullySaturated();

  // calculate and assemble elemental stiffness matrices
  AssembleSteadyStateMatrix(&*matrix_);
  AssembleSteadyStatePreconditioner(&*preconditioner_);
  preconditioner_->UpdatePreconditioner();

  // solve linear problem
  AmanziSolvers::LinearOperatorFactory<FlowMatrix, CompositeVector, CompositeVectorSpace> factory;
  Teuchos::RCP<AmanziSolvers::LinearOperator<FlowMatrix, CompositeVector, CompositeVectorSpace> >
     solver = factory.Create(ls_specs.solver_name, linear_operator_list_, matrix_, preconditioner_);

  const CompositeVector& rhs = *matrix_->rhs();
  int ierr = solver->ApplyInverse(rhs, u);

  if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH) {
    int num_itrs = solver->num_itrs();
    double residual = solver->residual();
    int code = solver->returned_code();

    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "saturated solver (" << solver->name() 
               << "): ||r||=" << residual << " itr=" << num_itrs 
               << " code=" << code << std::endl;
  }
  if (ierr != 0) {
    Errors::Message msg;
    msg << "\nLinear solver returned an unrecoverable error code.\n";
    Exceptions::amanzi_throw(msg);
  }
}


/* ******************************************************************
* Enforce constraints at time Tp by solving diagonalized MFD problem.
* Algorithm is based on de-coupling pressure-lambda system.
****************************************************************** */
void Richards_PK::EnforceConstraints(double Tp, CompositeVector& u)
{
  UpdateSourceBoundaryData(Tp, u);

  CompositeVector utmp(u);
  Epetra_MultiVector& utmp_faces = *utmp.ViewComponent("face");
  Epetra_MultiVector& u_faces = *u.ViewComponent("face");

  // calculate and assemble elemental stiffness matrix
  darcy_flux->ScatterMasterToGhosted("face");
  rel_perm->Compute(u, *darcy_flux, bc_model, bc_values);
  AssembleSteadyStateMatrix(&*matrix_);
  Matrix_MFD* matrix_tmp = dynamic_cast<Matrix_MFD*>(&*matrix_);
  matrix_tmp->ReduceGlobalSystem2LambdaSystem(u);

  // copy stiffness matrix to preconditioner (raw-data)
  Matrix_MFD* preconditioner_tmp = dynamic_cast<Matrix_MFD*>(&*preconditioner_);
  preconditioner_tmp->PopulatePreconditioner(*matrix_tmp);
  preconditioner_tmp->UpdatePreconditioner();

  // solve non-symmetric problem
  LinearSolver_Specs& ls_specs = ti_specs->ls_specs_constraints;

  AmanziSolvers::LinearOperatorFactory<FlowMatrix, CompositeVector, CompositeVectorSpace> factory;
  Teuchos::RCP<AmanziSolvers::LinearOperator<FlowMatrix, CompositeVector, CompositeVectorSpace> >
     solver = factory.Create(ls_specs.solver_name, linear_operator_list_, matrix_, preconditioner_);

  CompositeVector& rhs = *matrix_->rhs();
  int ierr = solver->ApplyInverse(rhs, utmp);

  u_faces = utmp_faces;

  if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH) {
    int num_itrs = solver->num_itrs();
    double residual = solver->residual();
    int code = solver->returned_code();

    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "constraints solver (" << solver->name() 
               << "): ||r||=" << residual << " itr=" << num_itrs
               << " code=" << code << std::endl;
  }
  if (ierr != 0) {
    Errors::Message msg;
    msg << "\nLinear solver returned an unrecoverable error code.\n";
    Exceptions::amanzi_throw(msg);
  }
}

}  // namespace AmanziFlow
}  // namespace Amanzi

