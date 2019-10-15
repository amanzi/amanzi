/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Konstantin Lipnikov (lipnikov@lanl.gov)
*/

//! <MISSING_ONELINE_DOCSTRING>

#include "Teuchos_RCP.hpp"
#include "Ifpack_Hypre.h"

#include "exceptions.hh"
#include "PreconditionerEuclid.hh"

namespace Amanzi {
namespace AmanziPreconditioners {

/* ******************************************************************
 * Apply the preconditioner.
 ****************************************************************** */
int
PreconditionerEuclid::ApplyInverse(const Epetra_MultiVector& v,
                                   Epetra_MultiVector& hv) const
{
  returned_code_ = IfpHypre_->ApplyInverse(v, hv);
  return returned_code_;
}


/* ******************************************************************
 * Initialize the preconditioner.
 ****************************************************************** */
void
PreconditionerEuclid::Init(const std::string& name,
                           const Teuchos::ParameterList& list)
{
  plist_ = list;
  funcs_.push_back(Teuchos::rcp(new FunctionParameter(
    (Hypre_Chooser)1, &HYPRE_EuclidSetStats, plist_.get<int>("verbosity", 0))));

  if (plist_.isParameter("ilu(k) fill level"))
    funcs_.push_back(Teuchos::rcp(
      new FunctionParameter((Hypre_Chooser)1,
                            &HYPRE_EuclidSetLevel,
                            plist_.get<int>("ilu(k) fill level"))));

  if (plist_.isParameter("rescale rows")) {
    bool rescale_rows = plist_.get<bool>("rescale rows");
    funcs_.push_back(Teuchos::rcp(new FunctionParameter(
      (Hypre_Chooser)1, &HYPRE_EuclidSetRowScale, rescale_rows ? 1 : 0)));
  }

  if (plist_.isParameter("ilut drop tolerance"))
    funcs_.push_back(Teuchos::rcp(
      new FunctionParameter((Hypre_Chooser)1,
                            &HYPRE_EuclidSetILUT,
                            plist_.get<double>("ilut drop tolerance"))));
}


/* ******************************************************************
 * Rebuild the preconditioner suing the given matrix A.
 ****************************************************************** */
void
PreconditionerEuclid::Update(const Teuchos::RCP<const Epetra_RowMatrix>& A)
{
  auto A_nc = Teuchos::rcp_const_cast<Epetra_RowMatrix>(A);
  IfpHypre_ = Teuchos::rcp(new Ifpack_Hypre(&*A_nc));

  Teuchos::ParameterList hypre_list("Preconditioner List");
  hypre_list.set("Preconditioner", Euclid);
  hypre_list.set("SolveOrPrecondition", (Hypre_Chooser)1);
  hypre_list.set("SetPreconditioner", true);
  hypre_list.set("NumFunctions", (int)funcs_.size());
  hypre_list.set<Teuchos::RCP<FunctionParameter>*>("Functions", &funcs_[0]);

  IfpHypre_->SetParameters(hypre_list);
  IfpHypre_->Initialize();
  IfpHypre_->Compute();
}

} // namespace AmanziPreconditioners
} // namespace Amanzi
