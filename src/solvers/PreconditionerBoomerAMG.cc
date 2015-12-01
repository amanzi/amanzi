/*
  This is the Linear Solver component of the Amanzi code.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)

  Conjugate gradient method.
*/

#include "Teuchos_RCP.hpp"
#include "Ifpack_Hypre.h"

#include "errors.hh"
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
  // returned_code_ = 0 means success. This is the only code returned by IfPack.
  return returned_code_;
}


/* ******************************************************************
 * Initialize the preconditioner.
 ****************************************************************** */
void PreconditionerBoomerAMG::Init(const std::string& name, const Teuchos::ParameterList& list)
{
  plist_ = list;
#ifdef HAVE_HYPRE
  // check for old input spec and error
  if (plist_.isParameter("number of cycles")) {
    Errors::Message msg("\"boomer amg\" ParameterList uses old style, \"number of cycles\".  Please update to the new style using \"cycle applications\" and \"smoother sweeps\"");
    Exceptions::amanzi_throw(msg);
  }
    
  funcs_.push_back(Teuchos::rcp(new FunctionParameter((Hypre_Chooser)1, &HYPRE_BoomerAMGSetTol,
						      plist_.get<double>("tolerance", 0.0))));
  funcs_.push_back(Teuchos::rcp(new FunctionParameter((Hypre_Chooser)1, &HYPRE_BoomerAMGSetPrintLevel,
						      plist_.get<int>("verbosity", 0))));
  funcs_.push_back(Teuchos::rcp(new FunctionParameter((Hypre_Chooser)1, &HYPRE_BoomerAMGSetMaxIter,
						      plist_.get<int>("cycle applications", 5))));
  funcs_.push_back(Teuchos::rcp(new FunctionParameter((Hypre_Chooser)1, &HYPRE_BoomerAMGSetCoarsenType,
						      plist_.get<int>("coarsen type", 0))));
  funcs_.push_back(Teuchos::rcp(new FunctionParameter((Hypre_Chooser)1, &HYPRE_BoomerAMGSetStrongThreshold,
						      plist_.get<double>("strong threshold", 0.5))));
  funcs_.push_back(Teuchos::rcp(new FunctionParameter((Hypre_Chooser)1, &HYPRE_BoomerAMGSetCycleType,
						      plist_.get<int>("cycle type", 1))));
  funcs_.push_back(Teuchos::rcp(new FunctionParameter((Hypre_Chooser)1, &HYPRE_BoomerAMGSetNumSweeps,
						      plist_.get<int>("smoother sweeps", 3))));
  if (plist_.isParameter("relaxation type down") && plist_.isParameter("relaxation type up")) {
    funcs_.push_back(Teuchos::rcp(new FunctionParameter((Hypre_Chooser)1, &HYPRE_BoomerAMGSetCycleRelaxType,
                                                        plist_.get<int>("relaxation type down"), 1)));
    funcs_.push_back(Teuchos::rcp(new FunctionParameter((Hypre_Chooser)1, &HYPRE_BoomerAMGSetCycleRelaxType,
                                                        plist_.get<int>("relaxation type up"), 2)));
  } else {
    funcs_.push_back(Teuchos::rcp(new FunctionParameter((Hypre_Chooser)1, &HYPRE_BoomerAMGSetRelaxType,
                                                        plist_.get<int>("relaxation type", 3))));
  }

  if (plist_.isParameter("max multigrid levels"))
    funcs_.push_back(Teuchos::rcp(new FunctionParameter((Hypre_Chooser)1, &HYPRE_BoomerAMGSetMaxLevels,
							plist_.get<int>("max multigrid levels"))));
  if (plist_.isParameter("max coarse size"))
    funcs_.push_back(Teuchos::rcp(new FunctionParameter((Hypre_Chooser)1, &HYPRE_BoomerAMGSetMaxCoarseSize,
							plist_.get<int>("max coarse size"))));

  if (plist_.isParameter("number of functions")) {
    // num_funcs > 1 tells BoomerAMG to use the "systems of PDEs" code.  Note
    // that, to use this approach, unknowns must be ordered with DoF fastest
    // varying (i.e. not the native Epetra_MultiVector order).  By default, it
    // uses the "unknown" approach in which each equation is coarsened and
    // interpolated independently.  Comments below are taken from Allison
    // Baker's email to the PETSc mailing list, 25 Apr 2007, as these features
    // of BoomerAMG are not documented very well.  Here we ignore her option
    // 2, as she warns it is inefficient and likely not useful.
    // http://lists.mcs.anl.gov/pipermail/petsc-users/2007-April/001487.html

    int num_funcs = plist_.get<int>("number of functions");
    funcs_.push_back(Teuchos::rcp(new FunctionParameter((Hypre_Chooser)1,
            &HYPRE_BoomerAMGSetNumFunctions, num_funcs)));

    if (plist_.isParameter("integer blocking")) {
      int* blocking = plist_.get<int*>("integer blocking");
      funcs_.push_back(Teuchos::rcp(new FunctionParameter((Hypre_Chooser)1, &HYPRE_BoomerAMGSetDofFunc, blocking)));
    }


    // additional options
    if (num_funcs > 1) {
      // HYPRE_BOOMERAMGSetNodal(solver, int nodal ) tells AMG to coarsen such
      // that each variable has the same coarse grid - sometimes this is more
      // "physical" for a particular problem. The value chosen here for nodal
      // determines how strength of connection is determined between the
      // coupled system.  I suggest setting nodal = 1, which uses a Frobenius
      // norm.  This does NOT tell AMG to use nodal relaxation.
      if (plist_.isParameter("nodal strength of connection norm")) {
        int nodal = plist_.get<int>("nodal strength of connection norm", 0);
        funcs_.push_back(Teuchos::rcp(new FunctionParameter((Hypre_Chooser)1,
                &HYPRE_BoomerAMGSetNodal, nodal)));
      }

      // You can additionally do nodal relaxation via the schwarz
      // smoother option in hypre.
      //   HYPRE_BoomerAMGSetSmoothType(solver, 6);
      //   HYPRE_BoomerAMGSetDomainType(solver, 1);
      //   HYPRE_BoomerAMGSetOverlap(solver, 0);
      //   HYPRE_BoomerAMGSetSmoothNumLevels(solver, num_levels);
      // Set num_levels to number of levels on which you want nodal
      // smoothing, i.e. 1=just the fine grid, 2= fine grid and the grid
      // below, etc.  I find that doing nodal relaxation on just the finest
      // level is generally sufficient.)  Note that the interpolation scheme
      // used will be the same as in the unknown approach - so this is what
      // we call a hybrid systems method.
      // 
      if (plist_.isParameter("nodal relaxation levels")) {
	int num_levels = plist_.get<int>("nodal relaxation levels");

	// I believe this works, but needs testing -- we do not pop previous
	// settings, and instead just call the function twice. --ETC
	funcs_.push_back(Teuchos::rcp(new FunctionParameter((Hypre_Chooser)1,
			 &HYPRE_BoomerAMGSetSmoothType, 6)));
	funcs_.push_back(Teuchos::rcp(new FunctionParameter((Hypre_Chooser)1,
                         &HYPRE_BoomerAMGSetDomainType, 1)));
	funcs_.push_back(Teuchos::rcp(new FunctionParameter((Hypre_Chooser)1,
                         &HYPRE_BoomerAMGSetOverlap, 0)));
	funcs_.push_back(Teuchos::rcp(new FunctionParameter((Hypre_Chooser)1,
                         &HYPRE_BoomerAMGSetSmoothNumLevels, num_levels)));
	funcs_.push_back(Teuchos::rcp(new FunctionParameter((Hypre_Chooser)1,
		         &HYPRE_BoomerAMGSetSchwarzUseNonSymm, 1))); // should provide an option for non-sym

	// Note that if num_levels > 1, you MUST also do nodal coarsening (to maintain the nodes on coarser grids).
	if (num_levels > 1) {
	  int nodal = plist_.get<int>("nodal strength of connection norm", 1);
	  funcs_.push_back(Teuchos::rcp(new FunctionParameter((Hypre_Chooser)1,
                           &HYPRE_BoomerAMGSetNodal, nodal)));
	}
      }
    }
  }

#else
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
  hypre_list.set("NumFunctions", (int)funcs_.size());
  hypre_list.set<Teuchos::RCP<FunctionParameter>*>("Functions", &funcs_[0]);

  IfpHypre_->SetParameters(hypre_list);
  IfpHypre_->Initialize();
  IfpHypre_->Compute();
#endif
}

}  // namespace AmanziPreconditioners
}  // namespace Amanzi
