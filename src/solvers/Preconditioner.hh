/*
This is the Linear Solver component of the Amanzi code.
 
License: BSD
Authors: Konstantin Lipnikov (lipnikov@lanl.gov)

Conjugate gradient method.
Usage: 
*/

#ifndef __PRECONDITIONER_HH__
#define __PRECONDITIONER_HH__

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Epetra_Vector.h"
#include "Epetra_FECrsMatrix.h"

#include "exceptions.hh"
 
namespace Amanzi {
namespace AmanziPreconditioners {

class Preconditioner {
 public:
  Preconditioner() {};
  ~Preconditioner() {};

  virtual void Init(const Teuchos::ParameterList& list) = 0;
  virtual void Update(Teuchos::RCP<Epetra_FECrsMatrix> A) = 0;
  virtual void Destroy() = 0;

  virtual void ApplyInverse(const Epetra_Vector& v, Epetra_Vector& hv) = 0;
};

}  // namespace AmanziPreconditioners
}  // namespace Amanzi


#endif


