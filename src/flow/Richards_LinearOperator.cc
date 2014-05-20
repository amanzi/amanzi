/*
This is the flow component of the Amanzi code. 

Copyright 2010-2012 held jointly by LANS/LANL, LBNL, and PNNL. 
Amanzi is released under the three-clause BSD License. 
The terms of use and "as is" disclaimer for this license are 
provided Reconstruction.cppin the top-level COPYRIGHT file.

Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include "OperatorDefs.hh"
#include "OperatorDiffusion.hh"
#include "OperatorDiffusionFactory.hh"
#include "LinearOperatorFactory.hh"

#include "Matrix_MFD.hh"
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

  Teuchos::ParameterList op_list;
  Teuchos::ParameterList& tmp_list = op_list.sublist("diffusion operator");
  tmp_list.set<std::string>("discretization primary", "optimized mfd scaled");
  tmp_list.set<std::string>("discretization secondary", "optimized mfd scaled");
  Teuchos::Array<std::string> stensil(2);
  stensil[0] = "face";
  stensil[1] = "cell";
  tmp_list.set<Teuchos::Array<std::string> >("schema", stensil);
  stensil.remove(1);
  tmp_list.set<Teuchos::Array<std::string> >("preconditioner schema", stensil);
  tmp_list.set<bool>("gravity", true);

  Operators::OperatorDiffusionFactory opfactory;
  Teuchos::RCP<Operators::OperatorDiffusion> op = opfactory.Create(mesh_, op_list, gravity_);
  op->InitOperator(K, Teuchos::null, Teuchos::null, rho_, mu_);
  op->UpdateMatrices(Teuchos::null);

  int schema_prec_dofs = op->schema_prec_dofs();
  op->SymbolicAssembleMatrix(schema_prec_dofs);

  // calculate and assemble elemental stifness matrices
  int n = bc_model.size();
  std::vector<double> bc_values_copy(n);
  for (int i = 0; i < n; i++) bc_values_copy[i] = bc_values[i][0];

  // add diffusion operator
  op->ApplyBCs(bc_model, bc_values_copy);
  op->AssembleMatrix(schema_prec_dofs);
  op->InitPreconditioner(ti_specs->preconditioner_name, preconditioner_list_, bc_model, bc_values_copy);

  AmanziSolvers::LinearOperatorFactory<Operators::OperatorDiffusion, CompositeVector, CompositeVectorSpace> sfactory;
  Teuchos::RCP<AmanziSolvers::LinearOperator<Operators::OperatorDiffusion, CompositeVector, CompositeVectorSpace> >
     solver = sfactory.Create(ls_specs.solver_name, linear_operator_list_, op);

  solver->add_criteria(AmanziSolvers::LIN_SOLVER_MAKE_ONE_ITERATION);  // Make at least one iteration

  CompositeVector& rhs = *op->rhs();
  int ierr = solver->ApplyInverse(rhs, u);

  if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH) {
    int num_itrs = solver->num_itrs();
    int code = solver->returned_code();
    double pnorm;
    u.Norm2(&pnorm);

    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "saturated solver (" << solver->name() 
               << "): ||p,lambda||=" << pnorm << " itr=" << num_itrs 
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
  darcy_flux_copy->ScatterMasterToGhosted("face");
  rel_perm->Compute(u, *darcy_flux_copy, bc_model, bc_values);
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

