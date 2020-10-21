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

// setting to null allows us to destroy the ML_ object
void
PreconditionerML::set_matrices(const Teuchos::RCP<Epetra_CrsMatrix>& m,
			       const Teuchos::RCP<Epetra_CrsMatrix>& h)
{
  if (h == Teuchos::null && ML_.get()) {
    ML_->DestroyPreconditioner();
  }
  Preconditioner::set_matrices(m,h);
}

/* ******************************************************************
 * Apply the preconditioner.
 * ML's return code is set to 0 if successful, see Trilinos webpages.
 ****************************************************************** */
int PreconditionerML::ApplyInverse(const Epetra_Vector& v, Epetra_Vector& hv) const
{
  returned_code_ = ML_->ApplyInverse(v, hv);
  return returned_code_;
}


/* ******************************************************************
 * Initialize the preconditioner.
 ****************************************************************** */
void PreconditionerML::set_inverse_parameters(Teuchos::ParameterList& list)
{
  list_ = list;
  list_.remove("verbose object", false); // note, ML validates parameter lists...
  list_.remove("method", false); // note, ML validates parameter lists...
  if (list_.isParameter("max levels") && list_.get<int>("max levels") > 10)
    list_.set("max levels", 10);
  initialized_ = true;
}


/* ******************************************************************
 * Rebuild the preconditioner suing the given matrix A.
 ****************************************************************** */
void PreconditionerML::InitializeInverse()
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
