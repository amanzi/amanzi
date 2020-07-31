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

  op_->ComputeInverse();
  int ierr = op_->ApplyInverse(rhs, *solution);

  if (vo_->os_OK(Teuchos::VERB_HIGH)) {
    int num_itrs = op_->num_itrs();
    double residual = op_->residual();
    int code = op_->returned_code();

    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "pressure solver (" << solver_name_
               << "): ||r||_H=" << residual << " itr=" << num_itrs
               << " code=" << code << std::endl;

    // Do you really want this just for reporting?  Seems extremely inefficient --etc
    //    residual = op_->TrueResidual(rhs, *solution);
    //    *vo_->os() << "true l2 residual: ||r||=" << residual << std::endl;
  }

  // catastrophic failure.
  if (ierr < 0) {
    Errors::Message msg("Transport_PK solver failed with message: \"");
    msg << op_->returned_code_string() << "\"";
    Exceptions::amanzi_throw(msg);
  }
}

}  // namespace Flow
}  // namespace Amanzi


