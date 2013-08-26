/*
This is the Linear Solver component of the Amanzi code.
 
License: BSD
Authors: Konstantin Lipnikov (lipnikov@lanl.gov)

Conjugate gradient method.
Usage: 
*/

#include "Teuchos_RCP.hpp"

#include "Preconditioner.hh"
#include "PreconditionerFactory.hh"
#include "PreconditionerIdentity.hh"
#include "PreconditionerHypre.hh"
 
namespace Amanzi {
namespace AmanziPreconditioners {

/* ******************************************************************
* Initialization of the preconditioner                                                 
****************************************************************** */
Teuchos::RCP<Preconditioner> PreconditionerFactory::Create(const string& name)
{
  if (name == "hypre") {
    Teuchos::RCP<PreconditionerHypre> prec = Teuchos::rcp(new PreconditionerHypre());
    return prec;
  }
  Teuchos::RCP<PreconditionerIdentity> prec = Teuchos::rcp(new PreconditionerIdentity());
  return prec;
}

}  // namespace AmanziPreconditioners
}  // namespace Amanzi



