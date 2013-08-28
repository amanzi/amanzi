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
void PreconditionerHypre::ApplyInverse(const Epetra_Vector& v, Epetra_Vector& hv)
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
  ncycles = list_.get<int>("cycle applications", 5);  // Boomer AMG parameters
  nsmooth = list_.get<int>("smoother sweeps", 3);
  tol = list_.get<double>("tolerance", 0.0);
  strong_threshold = list_.get<double>("strong threshold", 0.0);
  verbosity = list_.get<int>("verbosity", 0);
#endif
}


/* ******************************************************************
* Rebuild the preconditioner suing the given matrix A.                                                
****************************************************************** */
void PreconditionerHypre::Update(Teuchos::RCP<Epetra_FECrsMatrix> A)
{
#ifdef HAVE_HYPRE
  IfpHypre_ = Teuchos::rcp(new Ifpack_Hypre(&*A));

  Teuchos::RCP<FunctionParameter> functs[8];
  functs[0] = Teuchos::rcp(new FunctionParameter((Hypre_Chooser)1, &HYPRE_BoomerAMGSetCoarsenType, 0));
  functs[1] = Teuchos::rcp(new FunctionParameter((Hypre_Chooser)1, &HYPRE_BoomerAMGSetPrintLevel, verbosity)); 
  functs[2] = Teuchos::rcp(new FunctionParameter((Hypre_Chooser)1, &HYPRE_BoomerAMGSetNumSweeps, nsmooth));
  functs[3] = Teuchos::rcp(new FunctionParameter((Hypre_Chooser)1, &HYPRE_BoomerAMGSetMaxIter, ncycles));
  functs[4] = Teuchos::rcp(new FunctionParameter((Hypre_Chooser)1, &HYPRE_BoomerAMGSetRelaxType, 6)); 
  functs[5] = Teuchos::rcp(new FunctionParameter((Hypre_Chooser)1, &HYPRE_BoomerAMGSetStrongThreshold, strong_threshold)); 
  functs[6] = Teuchos::rcp(new FunctionParameter((Hypre_Chooser)1, &HYPRE_BoomerAMGSetTol, tol)); 
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



