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
#include "OperatorDiffusion.hh"

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
  int n = bc_model.size();
  std::vector<double> bc_values_copy(n);
  for (int i = 0; i < n; i++) bc_values_copy[i] = bc_values[i][0];

  // add diffusion operator
  int schema_dofs = op->schema_dofs();
  int schema_prec_dofs = op->schema_prec_dofs();

  op->RestoreCheckPoint();
  op->ApplyBCs(bc_model, bc_values_copy);
  op->AssembleMatrix(schema_prec_dofs);
  op->InitPreconditioner(ti_specs->preconditioner_name, preconditioner_list_, bc_model, bc_values_copy);

  AmanziSolvers::LinearOperatorFactory<Operators::OperatorDiffusion, CompositeVector, CompositeVectorSpace> sfactory;
  Teuchos::RCP<AmanziSolvers::LinearOperator<Operators::OperatorDiffusion, CompositeVector, CompositeVectorSpace> >
     solver = sfactory.Create(ti_specs->ls_specs.solver_name, linear_operator_list_, op);

  solver->add_criteria(AmanziSolvers::LIN_SOLVER_MAKE_ONE_ITERATION);  // Make at least one iteration

  CompositeVector& rhs = *op->rhs();
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
  if (ierr != 0) {
    Errors::Message msg;
    msg << "\nLinear solver returned an unrecoverable error code.\n";
    Exceptions::amanzi_throw(msg);
  }
}

}  // namespace AmanziFlow
}  // namespace Amanzi


