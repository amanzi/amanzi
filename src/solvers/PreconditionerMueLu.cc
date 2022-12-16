/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "PreconditionerMueLu.hh"

namespace Amanzi {
namespace AmanziSolvers {

/* ******************************************************************
* TBW
****************************************************************** */
void
PreconditionerMueLu::set_matrices(const Teuchos::RCP<Epetra_CrsMatrix>& m,
                                  const Teuchos::RCP<Epetra_CrsMatrix>& h)
{
  Preconditioner::set_matrices(m, h);
}


/* ******************************************************************
* TBW
****************************************************************** */
void
PreconditionerMueLu::InitializeInverse()
{
  plist_.remove("method");
  plist_.remove("verbose object");
}


void
PreconditionerMueLu::ComputeInverse()
{
#if defined(HAVE_MUELU_EPETRA)
  MueLu_ = MueLu::CreateEpetraPreconditioner(h_, plist_);
#endif
}


/* ******************************************************************
* TBW
****************************************************************** */
int
PreconditionerMueLu::ApplyInverse(const Vector_t& v, Vector_t& hv) const
{
#if defined(HAVE_MUELU_EPETRA)
  AMANZI_ASSERT(MueLu_.get());
  MueLu_->ApplyInverse(v, hv);
#endif
  returned_code_ = 0;
  return returned_code_;
}

} // namespace AmanziSolvers
} // namespace Amanzi
