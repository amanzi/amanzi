/*
  This is the flow component of the Amanzi code.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <vector>

#include "OperatorDefs.hh"
#include "OperatorDiffusion.hh"

#include "Darcy_PK.hh"
#include "LinearOperatorFactory.hh"

namespace Amanzi {
namespace Flow {

/* ******************************************************************
* Calculates steady-state solution assuming that absolute permeability 
* does not depend on time. The boundary conditions are calculated
* only once, during the initialization step.                                                
****************************************************************** */
void Darcy_PK::SolveFullySaturatedProblem(CompositeVector& u)
{
  // add diffusion operator
  op_->RestoreCheckPoint();
  op_diff_->ApplyBCs(true, true);
  op_->AssembleMatrix();
  op_->InitPreconditioner(preconditioner_name_, *preconditioner_list_);

  AmanziSolvers::LinearOperatorFactory<Operators::Operator, CompositeVector, CompositeVectorSpace> sfactory;
  Teuchos::RCP<AmanziSolvers::LinearOperator<Operators::Operator, CompositeVector, CompositeVectorSpace> >
     solver = sfactory.Create(solver_name_, *linear_operator_list_, op_);

  solver->add_criteria(AmanziSolvers::LIN_SOLVER_MAKE_ONE_ITERATION);  // Make at least one iteration

  CompositeVector& rhs = *op_->rhs();
  int ierr = solver->ApplyInverse(rhs, *solution);

  if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH) {
    int num_itrs = solver->num_itrs();
    double residual = solver->residual();
    int code = solver->returned_code();

    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "pressure solver (" << solver->name() 
               << "): ||r||=" << residual << " itr=" << num_itrs
               << " code=" << code << std::endl;
  }

  // if (ierr != 0) {
  //   Errors::Message msg;
  //   msg << "\nLinear solver returned an unrecoverable error code.\n";
  //   Exceptions::amanzi_throw(msg);
  // }
  if (ierr < 0) {
    Errors::Message msg;
    switch(ierr){
    case  Amanzi::AmanziSolvers::LIN_SOLVER_NON_SPD_APPLY:
      msg << "Linear system is not SPD.\n";
    case  Amanzi::AmanziSolvers::LIN_SOLVER_NON_SPD_APPLY_INVERSE:
      msg << "Linear system is not SPD.\n";
    case  Amanzi::AmanziSolvers::LIN_SOLVER_MAX_ITERATIONS:
      msg << "Maximum iterations are reached in solution of linear system.\n";
    case  Amanzi::AmanziSolvers::LIN_SOLVER_RESIDUAL_OVERFLOW:
      msg << "Residual overflow in solution of linear system.\n";
    default:
      msg << "\nLinear solver returned an unrecoverable error code: "<<ierr<<".\n";
    }
    Exceptions::amanzi_throw(msg);
  }
}

}  // namespace Flow
}  // namespace Amanzi


