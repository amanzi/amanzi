#include "DarcyMatvec.hpp"

DarcyMatvec::DarcyMatvec(DarcyProblem *problem) : problem_(problem)
{
  label_ = strdup("DarcyMatvec");
}


int DarcyMatvec::Apply(const Epetra_MultiVector &X, Epetra_MultiVector &Y) const
{
  const Epetra_Vector *xx = dynamic_cast<const Epetra_Vector*>(&X);
  Epetra_Vector *yy = dynamic_cast<Epetra_Vector*>(&Y);
  if (!xx || !yy) {
    std::cout << "DarcyMatvec::Apply() -- ERROR: the arguments are not Epetra_Vectors" << std::endl;
    throw;
  }
  problem_->ComputeF(*xx, *yy);
  yy->Update(1.0, problem_->RHS(), 1.0);
  return 0;
}
