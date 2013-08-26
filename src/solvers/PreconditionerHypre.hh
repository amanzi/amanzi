/*
This is the Linear Solver component of the Amanzi code.
 
License: BSD
Authors: Konstantin Lipnikov (lipnikov@lanl.gov)

Conjugate gradient method.
Usage: 
*/

#ifndef __PRECONDITIONER_HYPRE_HH__
#define __PRECONDITIONER_HYPRE_HH__

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Epetra_Vector.h"
#include "Epetra_FECrsMatrix.h"
#include "Ifpack.h" 

#include "exceptions.hh"
#include "Preconditioner.hh"
 
namespace Amanzi {
namespace AmanziPreconditioners {

class PreconditionerHypre : public Preconditioner {
 public:
  PreconditionerHypre() {};
  ~PreconditionerHypre() {};

  void Init(const Teuchos::ParameterList& list);
  void Update(Teuchos::RCP<Epetra_FECrsMatrix> A);
  void Destroy() {};

  void ApplyInverse(const Epetra_Vector& v, Epetra_Vector& hv);

 private:
  Teuchos::ParameterList list_;

  Teuchos::RCP<Ifpack_Hypre> IfpHypre_;
  double tol, strong_threshold;
  int nsmooth, ncycles, verbosity;
};

}  // namespace AmanziPreconditioners
}  // namespace Amanzi



#endif


