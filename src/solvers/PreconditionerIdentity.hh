/*
This is the Linear Solver component of the Amanzi code.
 
License: BSD
Authors: Konstantin Lipnikov (lipnikov@lanl.gov)

Conjugate gradient method.
Usage: 
*/

#ifndef __PRECONDITIONER_IDENTITY_HH__
#define __PRECONDITIONER_IDENTITY_HH__

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Epetra_Vector.h"
#include "Epetra_FECrsMatrix.h"

#include "exceptions.hh"
#include "Preconditioner.hh"
 
namespace Amanzi {
namespace AmanziPreconditioners {

class PreconditionerIdentity : public Preconditioner {
 public:
  PreconditionerIdentity() {};
  ~PreconditionerIdentity() {};

  void Init(const Teuchos::ParameterList& list) {};
  void Update(Teuchos::RCP<Epetra_FECrsMatrix> A) {};
  void Destroy() {};

  void ApplyInverse(const Epetra_Vector& v, Epetra_Vector& hv) { hv = v; }
};

}  // namespace AmanziPreconditioners
}  // namespace Amanzi



#endif


