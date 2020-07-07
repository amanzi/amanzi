/*
  Solvers

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Trilinos' multigrid solver.
*/

#include "Teuchos_RCP.hpp"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"

#include "exceptions.hh"
#include "PreconditionerML.hh"

namespace Amanzi {
namespace AmanziSolvers {

/* ******************************************************************
 * Apply the preconditioner. 
 * ML's return code is set to 0 if successful, see Trilinos webpages. 
 ****************************************************************** */
int PreconditionerML::ApplyInverse(const Epetra_Vector& v, Epetra_Vector& hv) const
{
  returned_code_ = ML_->ApplyInverse(v, hv);
  AMANZI_ASSERT(returned_code_ == 0);
  returned_code_ = 1;
  return returned_code_;
}


/* ******************************************************************
 * Initialize the preconditioner.
 ****************************************************************** */
void PreconditionerML::InitInverse(Teuchos::ParameterList& list)
{
  list_ = list;
  list_.remove("verbose object", false); // note, ML validates parameter lists...
  list_.remove("method", false); // note, ML validates parameter lists...
  initialized_ = true;
}


/* ******************************************************************
 * Rebuild the preconditioner suing the given matrix A.
 ****************************************************************** */
void PreconditionerML::UpdateInverse()
{
  AMANZI_ASSERT(initialized_);
  AMANZI_ASSERT(h_.get());
  ML_ = Teuchos::rcp(new ML_Epetra::MultiLevelPreconditioner(*h_, list_, false));
}

void PreconditionerML::ComputeInverse()
{
  AMANZI_ASSERT(ML_.get());
  ML_->ComputePreconditioner();
}


}  // namespace AmanziSolvers
}  // namespace Amanzi
