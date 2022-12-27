/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Flow PK

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
void
Darcy_PK::SolveFullySaturatedProblem(CompositeVector& u, bool wells_on)
{
  // add diffusion operator
  op_->RestoreCheckPoint();

  if (wells_on && S_->HasRecord("well_index")) {
    const auto& wi = S_->Get<CompositeVector>("well_index");
    op_acc_->AddAccumulationTerm(wi, "cell");
  }

  op_diff_->ApplyBCs(true, true, true);
  CompositeVector& rhs = *op_->rhs();
  if (wells_on) AddSourceTerms(rhs, 1.0);

  op_->ComputeInverse();
  int ierr = op_->ApplyInverse(rhs, *solution);

  if (vo_->os_OK(Teuchos::VERB_HIGH)) {
    int num_itrs = op_->num_itrs();

    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "pressure solver (" << solver_name_ << "): ||r||_H=" << op_->residual()
               << " itr=" << num_itrs << " code=" << op_->returned_code() << std::endl;

    // verify true resdual if the convergence was too slow
    if (num_itrs > 10) {
      CompositeVector r(*solution);
      op_->ComputeResidual(*solution, r);

      double true_residual;
      r.Norm2(&true_residual);
      *vo_->os() << "true l2 residual: ||r||=" << true_residual << std::endl;
    }
  }

  // catastrophic failure.
  if (ierr < 0) {
    Errors::Message msg("Transport_PK solver failed with message: \"");
    msg << op_->returned_code_string() << "\"";
    Exceptions::amanzi_throw(msg);
  }
}

} // namespace Flow
} // namespace Amanzi
