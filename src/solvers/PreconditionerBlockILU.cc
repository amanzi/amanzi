/*
  Solvers

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Incomplete LU preconditioner.
*/

#include "Teuchos_RCP.hpp"
#include "Ifpack.h"

#include "exceptions.hh"
#include "PreconditionerBlockILU.hh"

namespace Amanzi {
namespace AmanziPreconditioners {

/* ******************************************************************
* Apply the preconditioner.
* According to IfPack documentation, the error code is set to 0 if 
* the inversion was successful. 
****************************************************************** */
int PreconditionerBlockILU::ApplyInverse(const Epetra_MultiVector& v, Epetra_MultiVector& hv)
{
  returned_code_ = IfpILU_->ApplyInverse(v, hv);
  return returned_code_;
}


/* ******************************************************************
 * Initialize the preconditioner.
 ****************************************************************** */
void PreconditionerBlockILU::Init(const std::string& name, const Teuchos::ParameterList& list)
{
  list_ = list;
  initialized_ = false;
}


/* ******************************************************************
 * Rebuild the preconditioner suing the given matrix A.
 ****************************************************************** */
void PreconditionerBlockILU::Update(const Teuchos::RCP<Epetra_RowMatrix>& A)
{
  Ifpack factory;
  std::string type("ILU");

  int overlap = list_.get<int>("overlap", 0);
  list_.set<std::string>("schwarz: combine mode", "Add");

  IfpILU_ = Teuchos::rcp(factory.Create(type, &*A, overlap));

  IfpILU_->SetParameters(list_);
  IfpILU_->Initialize();
  IfpILU_->Compute();

  initialized_ = true;
}


/* ******************************************************************
 * Destroy the preconditioner and auxiliary data structures.
 ****************************************************************** */
void PreconditionerBlockILU::Destroy()
{
}

}  // namespace AmanziPreconditioners
}  // namespace Amanzi
