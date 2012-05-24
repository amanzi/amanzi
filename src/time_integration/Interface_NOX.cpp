/*
This is the flow component of the Amanzi code. 
License: BSD
Authors: Neil Carlson (version 1) 
         Konstantin Lipnikov (version 2) (lipnikov@lanl.gov)
*/

// #include "dbc.hh"

// #include "Flow_PK.hpp"
#include "Interface_NOX.hpp"

namespace Amanzi {
namespace AmanziFlow {

/* ******************************************************************
* NOT USED
****************************************************************** */
bool Interface_NOX::computeF(const Epetra_Vector& x, Epetra_Vector& f, FillType flag)
{
  for (int i=0; i<x.MyLength(); i++) f[i] = x[i]*x[i] - 9;
  return true;
}


/* ******************************************************************
* NOT USED
****************************************************************** */
bool Interface_NOX::computePreconditioner(
    const Epetra_Vector& x, Epetra_Operator& M, Teuchos::ParameterList* params)
{
  // We assume the input operator is the same one we handed to NOX.
  //ASSERT(&M == &(FPK_->NOX_Precon())); XXX

//   if (lag_count_ == 0) FPK_->NOX_ComputePrecon(x);
  lag_count_++;
  lag_count_ %= lag_prec_;

  return true;
}

}  // namespace AmanziFlow
}  // namespace Amanzi

