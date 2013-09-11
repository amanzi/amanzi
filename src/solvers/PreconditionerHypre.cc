/*
  This is the Linear Solver component of the Amanzi code.

  License: BSD
  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)

  Conjugate gradient method.
  Usage:
*/

#include "Teuchos_RCP.hpp"
#include "Ifpack.h"

#include "exceptions.hh"
#include "PreconditionerHypre.hh"

namespace Amanzi {
namespace AmanziPreconditioners {

/* ******************************************************************
 * Apply the preconditioner.
 ****************************************************************** */
void PreconditionerHypre::ApplyInverse(const Epetra_MultiVector& v, Epetra_MultiVector& hv)
{
  IfpHypre_->ApplyInverse(v, hv);
}


/* ******************************************************************
 * Initialize the preconditioner.
 ****************************************************************** */
void PreconditionerHypre::Init(const std::string& name, const Teuchos::ParameterList& list)
{
  list_ = list;
#ifdef HAVE_HYPRE
  ncycles_ = list_.get<int>("cycle applications", 5);  // Boomer AMG parameters
  nsmooth_ = list_.get<int>("smoother sweeps", 3);
  tol_ = list_.get<double>("tolerance", 0.0);
  strong_threshold_ = list_.get<double>("strong threshold", 0.0);
  verbosity_ = list_.get<int>("verbosity", 0);
#else
  Errors::Message msg("Hypre is not available in this installation of Amanzi.  To use Hypre, please reconfigure.");
  Exceptions::amanzi_throw(msg);
#endif
}


/* ******************************************************************
 * Rebuild the preconditioner suing the given matrix A.
 ****************************************************************** */
void PreconditionerHypre::Update(const Teuchos::RCP<Epetra_RowMatrix>& A)
{
#ifdef HAVE_HYPRE
  IfpHypre_ = Teuchos::rcp(new Ifpack_Hypre(&*A));

  Teuchos::RCP<FunctionParameter> functs[8];
  functs[0] = Teuchos::rcp(new FunctionParameter((Hypre_Chooser)1, &HYPRE_BoomerAMGSetCoarsenType, 0));
  functs[1] = Teuchos::rcp(new FunctionParameter((Hypre_Chooser)1, &HYPRE_BoomerAMGSetPrintLevel, verbosity_));
  functs[2] = Teuchos::rcp(new FunctionParameter((Hypre_Chooser)1, &HYPRE_BoomerAMGSetNumSweeps, nsmooth_));
  functs[3] = Teuchos::rcp(new FunctionParameter((Hypre_Chooser)1, &HYPRE_BoomerAMGSetMaxIter, ncycles_));
  functs[4] = Teuchos::rcp(new FunctionParameter((Hypre_Chooser)1, &HYPRE_BoomerAMGSetRelaxType, 6));
  functs[5] = Teuchos::rcp(new FunctionParameter((Hypre_Chooser)1, &HYPRE_BoomerAMGSetStrongThreshold, strong_threshold_));
  functs[6] = Teuchos::rcp(new FunctionParameter((Hypre_Chooser)1, &HYPRE_BoomerAMGSetTol, tol_));
  functs[7] = Teuchos::rcp(new FunctionParameter((Hypre_Chooser)1, &HYPRE_BoomerAMGSetCycleType, 1));

  Teuchos::ParameterList hypre_list("Preconditioner List");
  hypre_list.set("Preconditioner", BoomerAMG);
  hypre_list.set("SolveOrPrecondition", (Hypre_Chooser)1);
  hypre_list.set("SetPreconditioner", true);
  hypre_list.set("NumFunctions", 8);
  hypre_list.set<Teuchos::RCP<FunctionParameter>*>("Functions", functs);

  IfpHypre_->SetParameters(hypre_list);
  IfpHypre_->Initialize();
  IfpHypre_->Compute();
#endif
}

}  // namespace AmanziPreconditioners
}  // namespace Amanzi
