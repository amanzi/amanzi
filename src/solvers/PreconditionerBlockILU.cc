/*
This is the Linear Solver component of the Amanzi code.
 
License: BSD
Authors: Konstantin Lipnikov (lipnikov@lanl.gov)

Incomplete LU preconditioner.
Usage: 
*/

#include "Teuchos_RCP.hpp"
#include "Ifpack.h"

#include "exceptions.hh"
#include "PreconditionerBlockILU.hh"
 
namespace Amanzi {
namespace AmanziPreconditioners {

/* ******************************************************************
* Apply the preconditioner.                                                 
****************************************************************** */
void PreconditionerBlockILU::ApplyInverse(const Epetra_Vector& v, Epetra_Vector& hv)
{
  IfpILU_->ApplyInverse(v, hv);
}


/* ******************************************************************
* Initialize the preconditioner.                                                 
****************************************************************** */
void PreconditionerBlockILU::Init(const std::string& name, const Teuchos::ParameterList& list)
{
  list_ = list;
  initialized = false;
}


/* ******************************************************************
* Rebuild the preconditioner suing the given matrix A.                                                
****************************************************************** */
void PreconditionerBlockILU::Update(Teuchos::RCP<Epetra_FECrsMatrix> A)
{
  Ifpack factory;
  std::string type("ILU");

  int overlap = list_.get<int>("overlap", 0);
  list_.set<std::string>("schwarz: combine mode", "Add");

  IfpILU_ = Teuchos::rcp(factory.Create(type, &*A, overlap));

  IfpILU_->SetParameters(list_);
  IfpILU_->Initialize();
  IfpILU_->Compute();

  initialized = true;
}


/* ******************************************************************
* Destroy the preconditioner and auxiliary data structures.
****************************************************************** */
void PreconditionerBlockILU::Destroy()
{
}

}  // namespace AmanziPreconditioners
}  // namespace Amanzi



