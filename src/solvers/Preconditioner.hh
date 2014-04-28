/*
This is the Linear Solver component of the Amanzi code.

License: BSD
Authors: Konstantin Lipnikov (lipnikov@lanl.gov)

Base class for preconditioners.

Usage:
*/

#ifndef AMANZI_PRECONDITIONER_HH_
#define AMANZI_PRECONDITIONER_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Epetra_MultiVector.h"
#include "Epetra_RowMatrix.h"

#include "exceptions.hh"

namespace Amanzi {
namespace AmanziPreconditioners {

class Preconditioner {
 public:
  Preconditioner() {};
  ~Preconditioner() {};

  // Initializes the solver with provided parameters.
  // This need not be called by preconditioners created using the factory.
  virtual void Init(const std::string& name,
                    const Teuchos::ParameterList& list) = 0;

  // Rebuild the preconditioner using the given matrix A.
  virtual void Update(const Teuchos::RCP<Epetra_RowMatrix>& A) = 0;

  // Destroy the preconditioner and auxiliary data structures.
  virtual void Destroy() = 0;

  // Apply the preconditioner.
  virtual int ApplyInverse(const Epetra_MultiVector& v,
			   Epetra_MultiVector& hv) = 0;

  virtual int returned_code() = 0;
};

}  // namespace AmanziPreconditioners
}  // namespace Amanzi


#endif


