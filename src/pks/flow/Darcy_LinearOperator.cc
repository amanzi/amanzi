/*
  Flow PK

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <vector>

#include "OperatorDefs.hh"
#include "PDE_Diffusion.hh"

#include "Darcy_PK.hh"

namespace Amanzi {
namespace Flow {

/* ******************************************************************
* Calculates steady-state solution assuming that absolute permeability 
* does not depend on time. The boundary conditions are calculated
* only once, during the initialization step.                                                
****************************************************************** */
void Darcy_PK::SolveFullySaturatedProblem(CompositeVector& u, bool wells_on)
{
  // add diffusion operator
  op_->RestoreCheckPoint();
 
  if (wells_on && S_->HasField("well_index")) {
    const CompositeVector& wi = *S_->GetFieldData("well_index");
    op_acc_->AddAccumulationTerm(wi, "cell");
  }

  op_diff_->ApplyBCs(true, true, true);
  CompositeVector& rhs = *op_->rhs();
  if (wells_on) AddSourceTerms(rhs);

  op_->AssembleMatrix();
  op_->UpdatePreconditioner();

     solver = sfactory.Create(solver_name_, *linear_operator_list_, op_);

  solver->add_criteria(AmanziSolvers::LIN_SOLVER_MAKE_ONE_ITERATION);  // Make at least one iteration

  int ierr = solver->ApplyInverse(rhs, *solution);

  if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH) {
    int num_itrs = solver->num_itrs();
    double residual = solver->residual();
    int code = solver->returned_code();

    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "pressure solver (" << solver->name() 
               << "): ||r||_H=" << residual << " itr=" << num_itrs
               << " code=" << code << std::endl;

    residual = solver->TrueResidual(rhs, *solution);
    *vo_->os() << "true l2 residual: ||r||=" << residual << std::endl;
  }

  // catastrophic failure.
  if (ierr < 0) {
    Errors::Message msg;
    msg = solver->DecodeErrorCode(ierr);
    Exceptions::amanzi_throw(msg);
  }
}

}  // namespace Flow
}  // namespace Amanzi


