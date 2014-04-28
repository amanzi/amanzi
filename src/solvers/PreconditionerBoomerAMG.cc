/*
  This is the Linear Solver component of the Amanzi code.

  License: BSD
  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)

  Conjugate gradient method.
  Usage:
*/

#include "Teuchos_RCP.hpp"
#include "Ifpack_Hypre.h"

#include "exceptions.hh"
#include "PreconditionerBoomerAMG.hh"

namespace Amanzi {
namespace AmanziPreconditioners {

/* ******************************************************************
 * Apply the preconditioner.
 ****************************************************************** */
int PreconditionerBoomerAMG::ApplyInverse(const Epetra_MultiVector& v, Epetra_MultiVector& hv)
{
  returned_code_ = IfpHypre_->ApplyInverse(v, hv);
  return (returned_code_ == 0) ? 0 : 1;
}


/* ******************************************************************
 * Initialize the preconditioner.
 ****************************************************************** */
void PreconditionerBoomerAMG::Init(const std::string& name, const Teuchos::ParameterList& list)
{
  plist_ = list;
#ifdef HAVE_HYPRE
  funcs_.push_back(Teuchos::rcp(new FunctionParameter((Hypre_Chooser)1, &HYPRE_BoomerAMGSetTol,
          plist_.get<double>("tolerance", 0.0))));

  funcs_.push_back(Teuchos::rcp(new FunctionParameter((Hypre_Chooser)1, &HYPRE_BoomerAMGSetPrintLevel,
          plist_.get<int>("verbosity", 0))));

  if (plist_.isParameter("cycle applications"))
    funcs_.push_back(Teuchos::rcp(new FunctionParameter((Hypre_Chooser)1, &HYPRE_BoomerAMGSetMaxIter,
            plist_.get<int>("cycle applications"))));

  if (plist_.isParameter("max multigrid levels"))
    funcs_.push_back(Teuchos::rcp(new FunctionParameter((Hypre_Chooser)1, &HYPRE_BoomerAMGSetMaxLevels,
            plist_.get<int>("max multigrid levels"))));
  
  if (plist_.isParameter("coarsen type"))
    funcs_.push_back(Teuchos::rcp(new FunctionParameter((Hypre_Chooser)1, &HYPRE_BoomerAMGSetCoarsenType,
            plist_.get<int>("coarsen type", 0))));

  if (plist_.isParameter("max coarse size"))
    funcs_.push_back(Teuchos::rcp(new FunctionParameter((Hypre_Chooser)1, &HYPRE_BoomerAMGSetMaxCoarseSize,
            plist_.get<int>("max coarse size"))));

  if (plist_.isParameter("strong threshold"))
    funcs_.push_back(Teuchos::rcp(new FunctionParameter((Hypre_Chooser)1, &HYPRE_BoomerAMGSetStrongThreshold,
            plist_.get<double>("strong threshold"))));

  if (plist_.isParameter("cycle type"))
    funcs_.push_back(Teuchos::rcp(new FunctionParameter((Hypre_Chooser)1, &HYPRE_BoomerAMGSetCycleType,
            plist_.get<int>("cycle type"))));

  if (plist_.isParameter("smoother sweeps"))
    funcs_.push_back(Teuchos::rcp(new FunctionParameter((Hypre_Chooser)1, &HYPRE_BoomerAMGSetNumSweeps,
            plist_.get<int>("smoother sweeps"))));

  if (plist_.isParameter("relaxation type"))
    funcs_.push_back(Teuchos::rcp(new FunctionParameter((Hypre_Chooser)1, &HYPRE_BoomerAMGSetNumSweeps,
            plist_.get<int>("relaxation type"))));

#else HAVE_HYPRE
  Errors::Message msg("Hypre (BoomerAMG) is not available in this installation of Amanzi.  To use Hypre, please reconfigure.");
  Exceptions::amanzi_throw(msg);
#endif
}


/* ******************************************************************
 * Rebuild the preconditioner using the given matrix A.
 ****************************************************************** */
void PreconditionerBoomerAMG::Update(const Teuchos::RCP<Epetra_RowMatrix>& A)
{
#ifdef HAVE_HYPRE
  IfpHypre_ = Teuchos::rcp(new Ifpack_Hypre(&*A));

  Teuchos::ParameterList hypre_list("Preconditioner List");
  hypre_list.set("Preconditioner", BoomerAMG);
  hypre_list.set("SolveOrPrecondition", (Hypre_Chooser)1);
  hypre_list.set("SetPreconditioner", true);
  hypre_list.set("NumFunctions", funcs_.size());
  hypre_list.set<Teuchos::RCP<FunctionParameter>*>("Functions", &funcs_[0]);

  IfpHypre_->SetParameters(hypre_list);
  IfpHypre_->Initialize();
  IfpHypre_->Compute();
#endif
}

}  // namespace AmanziPreconditioners
}  // namespace Amanzi
