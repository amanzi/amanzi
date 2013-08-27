/*
This is the Linear Solver component of the Amanzi code.
 
License: BSD
Authors: Konstantin Lipnikov (lipnikov@lanl.gov)

Conjugate gradient method.
Usage: 
*/

#include "Teuchos_RCP.hpp"

#include "exceptions.hh"
#include "PreconditionerML.hh"
 
namespace Amanzi {
namespace AmanziPreconditioners {

/* ******************************************************************
* Apply the preconditioner.                                                 
****************************************************************** */
void PreconditionerML::ApplyInverse(const Epetra_Vector& v, Epetra_Vector& hv)
{
  ML_->ApplyInverse(v, hv);
}


/* ******************************************************************
* Initialize the preconditioner.                                                 
****************************************************************** */
void PreconditionerML::Init(const std::string& name, const Teuchos::ParameterList& list)
{
  list_ = list;
  initialized = false;
}


/* ******************************************************************
* Rebuild the preconditioner suing the given matrix A.                                                
****************************************************************** */
void PreconditionerML::Update(Teuchos::RCP<Epetra_FECrsMatrix> A)
{
  if (initialized) ML_->DestroyPreconditioner();
  ML_ = Teuchos::rcp(new ML_Epetra::MultiLevelPreconditioner(*A, list_, false));
  ML_->ComputePreconditioner();
  initialized = true;
}


/* ******************************************************************
* Destroy the preconditioner and auxiliary data structures.
****************************************************************** */
void PreconditionerML::Destroy()
{
  if (initialized) ML_->DestroyPreconditioner();
  initialized = false;
}

}  // namespace AmanziPreconditioners
}  // namespace Amanzi



