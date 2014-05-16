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
#include "LinearOperatorFactory.hh"

#include "OperatorDefs.hh"
#include "OperatorDiffusionFactory.hh"
#include "OperatorDiffusion.hh"
#include "OperatorGravity.hh"

namespace Amanzi {
namespace AmanziFlow {

/* ******************************************************************
* Calculates steady-state solution assuming that absolute permeability 
* does not depend on time. The boundary conditions are calculated
* only once, during the initialization step.                                                
****************************************************************** */
void Darcy_PK::SolveFullySaturatedProblem(double Tp, CompositeVector& u)
{
  // calculate and assemble elemental stifness matrices
#if 0
  Teuchos::ParameterList op_list;
  Teuchos::ParameterList& tmp_list = op_list.sublist("diffusion operator");
  tmp_list.set<std::string>("discretization primary", "optimized mfd scaled");
  tmp_list.set<std::string>("discretization secondary", "optimized mfd scaled");
  Teuchos::Array<std::string> stensil(2);
  stensil[0] = "cell";
  stensil[1] = "face";
  tmp_list.set<Teuchos::Array<std::string> >("schema", stensil);

  // add diffusion operator
  Operators::OperatorDiffusionFactory opfactory;
  Teuchos::RCP<Operators::OperatorDiffusion> op1 = opfactory.Create(mesh_, op_list);
  int schema_dofs = op1->schema_dofs();

  op1->InitOperator(K, Teuchos::null);
  op1->UpdateMatrices(Teuchos::null);

  // add gravity operator
  AmanziGeometry::Point rho_g(gravity_);
  rho_g *= rho_;

  int n = bc_model.size();
  std::vector<double> bc_values_copy(n);
  for (int i = 0; i < n ; i++) bc_values_copy[i] = bc_values[i][0];

  Teuchos::RCP<Operators::OperatorGravity> op2 = Teuchos::rcp(new Operators::OperatorGravity(*op1));
  op2->UpdateMatrices(K, rho_g);
  op2->SymbolicAssembleMatrix(schema_dofs);
  op2->AssembleMatrix(schema_dofs);
  op2->ApplyBCs(bc_model, bc_values_copy);
  op2->InitPreconditioner(ti_specs->preconditioner_name, preconditioner_list_);

  AmanziSolvers::LinearOperatorFactory<Operators::OperatorGravity, CompositeVector, CompositeVectorSpace> sfactory;
  Teuchos::RCP<AmanziSolvers::LinearOperator<Operators::OperatorGravity, CompositeVector, CompositeVectorSpace> >
     solver = sfactory.Create(ti_specs->ls_specs.solver_name, linear_operator_list_, op2);

  solver->add_criteria(AmanziSolvers::LIN_SOLVER_MAKE_ONE_ITERATION);  // Make at least one iteration

  CompositeVector& rhs = *op2->rhs();
  int ierr = solver->ApplyInverse(rhs, *solution);
#endif

#if 1
  matrix_->CreateStiffnessMatricesDarcy(mfd3d_method_);
  matrix_->CreateRHSVectors();
  matrix_->AddGravityFluxesDarcy(rho_, gravity_);
  matrix_->ApplyBoundaryConditions(bc_model, bc_values);
  matrix_->Assemble();
  matrix_->AssembleDerivatives(u, bc_model, bc_values);
  matrix_->UpdatePreconditioner();

  // create linear solver
  LinearSolver_Specs& ls_specs = ti_specs->ls_specs;

  AmanziSolvers::LinearOperatorFactory<FlowMatrix, CompositeVector, CompositeVectorSpace> factory;
  Teuchos::RCP<AmanziSolvers::LinearOperator<FlowMatrix, CompositeVector, CompositeVectorSpace> >
     solver = factory.Create(ls_specs.solver_name, linear_operator_list_, matrix_);

  CompositeVector& rhs = *matrix_->rhs();
  int ierr = solver->ApplyInverse(rhs, *solution);
#endif

  if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH) {
    int num_itrs = solver->num_itrs();
    double residual = solver->residual();
    int code = solver->returned_code();

    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "pressure solver (" << solver->name() 
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


