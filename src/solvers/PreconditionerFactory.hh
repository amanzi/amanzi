/*
This is the Linear Solver component of the Amanzi code.
 
License: BSD
Authors: Konstantin Lipnikov (lipnikov@lanl.gov)

Conjugate gradient method.
Usage: 
*/

#ifndef __PRECONDITIONER_FACTORY_HH__
#define __PRECONDITIONER_FACTORY_HH__

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Epetra_Vector.h"

#include "Preconditioner.hh"
#include "PreconditionerFactory.hh"
#include "PreconditionerHypre.hh"
 
namespace Amanzi {
namespace AmanziPreconditioners {

class PreconditionerFactory {
 public:
  PreconditionerFactory() {};
  ~PreconditionerFactory() {};

  Teuchos::RCP<Preconditioner> Create(const string& name);
};

}  // namespace AmanziPreconditioners
}  // namespace Amanzi

#endif


