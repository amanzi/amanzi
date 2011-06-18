#include "dbc.hh"
#include "RichardsNoxInterface.hpp"

namespace Amanzi
{

bool RichardsNoxInterface::computeF(const Epetra_Vector &x, Epetra_Vector &f, FillType flag)
{
  (*problem_).ComputeF(x, f);
  return true;
}


bool RichardsNoxInterface::computeJacobian(const Epetra_Vector &x, Epetra_Operator &J)
{
  // Shouldn't be called -- not required for JFNK.
  ASSERT(false);
}


bool RichardsNoxInterface::computePreconditioner(const Epetra_Vector &x, Epetra_Operator &M, Teuchos::ParameterList *params)
{
  // We assume the input operator is the same one we handed to NOX.
  ASSERT(&M == &(problem_->Precon()));

  if (lag_count_ == 0) (*problem_).ComputePrecon(x);
  lag_count_++;
  lag_count_ %= lag_prec_;

  return true;
}


} // close namespace Amanzi
