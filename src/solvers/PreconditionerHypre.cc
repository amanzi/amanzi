/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
      Ethan Coon (coonet@ornl.gov)
*/

//! Hypre based preconditioners include Algebraic MultiGrid and global ILU
#include "HYPRE_parcsr_ls.h"
#include "Ifpack_Hypre.h"
#include "Teuchos_RCP.hpp"

#include "dbc.hh"
#include "errors.hh"
#include "exceptions.hh"
#include "VerboseObject.hh"
#include "PreconditionerHypre.hh"

namespace Amanzi {
namespace AmanziSolvers {

/* ******************************************************************
* Apply the preconditioner.
****************************************************************** */
int
PreconditionerHypre::ApplyInverse(const Epetra_Vector& v, Epetra_Vector& hv) const
{
  returned_code_ = IfpHypre_->ApplyInverse(v, hv);
  AMANZI_ASSERT(returned_code_ == 0);
  return returned_code_;
}


/* ******************************************************************
* Initialize the preconditioner.
****************************************************************** */
void
PreconditionerHypre::set_inverse_parameters(Teuchos::ParameterList& list)
{
  plist_ = list;

  std::string vo_name = this->name() + " (" + plist_.get<std::string>("method") + ")";
  vo_ = Teuchos::rcp(new VerboseObject(vo_name, plist_));
}


/* ******************************************************************
* Initialize AMG preconditioner.
****************************************************************** */
void
PreconditionerHypre::InitBoomer_()
{
#ifdef HAVE_HYPRE
  method_ = BoomerAMG;

  // check for old input spec and error
  if (plist_.isParameter("number of cycles")) {
    Errors::Message msg(
      "\"boomer amg\" ParameterList uses old style, \"number of cycles\".  Please update to the "
      "new style using \"cycle applications\" and \"smoother sweeps\"");
    Exceptions::amanzi_throw(msg);
  }

  // verbosity
  int vlevel_int = 0;
  std::string vlevel =
    plist_.sublist("verbose object").get<std::string>("verbosity level", "medium");
  if (vlevel == "high") {
    vlevel_int = 1;
  } else if (vlevel == "extreme") {
    vlevel_int = 3;
  }

  if (plist_.isParameter("verbosity")) vlevel_int = plist_.get<int>("verbosity");
  IfpHypre_->SetParameter((Hypre_Chooser)1, &HYPRE_BoomerAMGSetPrintLevel, vlevel_int);
  IfpHypre_->SetParameter(
    (Hypre_Chooser)1, &HYPRE_BoomerAMGSetTol, plist_.get<double>("tolerance", 0.0));
  IfpHypre_->SetParameter(
    (Hypre_Chooser)1, &HYPRE_BoomerAMGSetMaxIter, plist_.get<int>("cycle applications", 5));
  IfpHypre_->SetParameter(
    (Hypre_Chooser)1, &HYPRE_BoomerAMGSetCoarsenType, plist_.get<int>("coarsen type", 0));
  IfpHypre_->SetParameter((Hypre_Chooser)1,
                          &HYPRE_BoomerAMGSetStrongThreshold,
                          plist_.get<double>("strong threshold", 0.5));
  IfpHypre_->SetParameter(
    (Hypre_Chooser)1, &HYPRE_BoomerAMGSetCycleType, plist_.get<int>("cycle type", 1));
  IfpHypre_->SetParameter(
    (Hypre_Chooser)1, &HYPRE_BoomerAMGSetNumSweeps, plist_.get<int>("smoother sweeps", 3));

  if (plist_.isParameter("relaxation type down") && plist_.isParameter("relaxation type up")) {
    IfpHypre_->SetParameter((Hypre_Chooser)1,
                            &HYPRE_BoomerAMGSetCycleRelaxType,
                            plist_.get<int>("relaxation type down"),
                            1);
    IfpHypre_->SetParameter((Hypre_Chooser)1,
                            &HYPRE_BoomerAMGSetCycleRelaxType,
                            plist_.get<int>("relaxation type up"),
                            2);
  } else if (plist_.isParameter("relaxation type")) {
    IfpHypre_->SetParameter(
      (Hypre_Chooser)1, &HYPRE_BoomerAMGSetRelaxType, plist_.get<int>("relaxation type"));
  } else {
    // use Hypre's defaults
  }

  if (plist_.isParameter("max multigrid levels"))
    IfpHypre_->SetParameter(
      (Hypre_Chooser)1, &HYPRE_BoomerAMGSetMaxLevels, plist_.get<int>("max multigrid levels"));
  if (plist_.isParameter("max coarse size"))
    IfpHypre_->SetParameter(
      (Hypre_Chooser)1, &HYPRE_BoomerAMGSetMaxCoarseSize, plist_.get<int>("max coarse size"));

  if (plist_.get<bool>("use block indices", false)) {
    num_blocks_ = plist_.get<int>("number of unique block indices");
    IfpHypre_->SetParameter((Hypre_Chooser)1, &HYPRE_BoomerAMGSetNumFunctions, num_blocks_);

    // Block indices is an array of integers, indicating what unknowns are
    // coarsened as a system.
    block_indices_ = plist_.get<Teuchos::RCP<std::vector<int>>>("block indices");
  }

  if (plist_.isParameter("number of functions")) {
    if (num_blocks_ > 0) {
      Errors::Message msg(
        "Hypre (BoomerAMG) cannot be given both \"use block indices\" and \"number of functions\" "
        "options as these are two ways of specifying the same thing.");
      Exceptions::amanzi_throw(msg);
    }

    // num_funcs > 1 tells BoomerAMG to use the automatic "systems of
    // PDEs" code.  Note that, to use this approach, unknowns must be
    // ordered with DoF fastest varying (i.e. not the native
    // Epetra_MultiVector order).  By default, it uses the "unknown"
    // approach in which each equation is coarsened and interpolated
    // independently.  Comments below are taken from Allison Baker's
    // email to the PETSc mailing list, 25 Apr 2007, as these features
    // of BoomerAMG are not documented very well.  Here we ignore her
    // option 2, as she warns it is inefficient and likely not useful.
    // http://lists.mcs.anl.gov/pipermail/petsc-users/2007-April/001487.html

    int num_funcs = plist_.get<int>("number of functions");
    IfpHypre_->SetParameter((Hypre_Chooser)1, &HYPRE_BoomerAMGSetNumFunctions, num_funcs);

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
        IfpHypre_->SetParameter((Hypre_Chooser)1, &HYPRE_BoomerAMGSetNodal, nodal);
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
        IfpHypre_->SetParameter((Hypre_Chooser)1, &HYPRE_BoomerAMGSetSmoothType, 6);
        IfpHypre_->SetParameter((Hypre_Chooser)1, &HYPRE_BoomerAMGSetDomainType, 1);
        IfpHypre_->SetParameter((Hypre_Chooser)1, &HYPRE_BoomerAMGSetOverlap, 0);
        IfpHypre_->SetParameter((Hypre_Chooser)1, &HYPRE_BoomerAMGSetSmoothNumLevels, num_levels);
        IfpHypre_->SetParameter((Hypre_Chooser)1,
                                &HYPRE_BoomerAMGSetSchwarzUseNonSymm,
                                1); // should provide an option for non-sym

        // Note that if num_levels > 1, you MUST also do nodal coarsening (to maintain the nodes on coarser grids).
        if (num_levels > 1) {
          int nodal = plist_.get<int>("nodal strength of connection norm", 1);
          IfpHypre_->SetParameter((Hypre_Chooser)1, &HYPRE_BoomerAMGSetNodal, nodal);
        }
      }
    }
  }

#else
  Errors::Message msg("Hypre (BoomerAMG) is not available in this installation of Amanzi.  To use "
                      "Hypre, please reconfigure.");
  Exceptions::amanzi_throw(msg);
#endif
}


/* ******************************************************************
* Initialize modified ILU preconditioner.
****************************************************************** */
void
PreconditionerHypre::InitEuclid_()
{
#ifdef HAVE_HYPRE
  method_ = Euclid;

  if (plist_.isParameter("verbosity"))
    IfpHypre_->SetParameter((Hypre_Chooser)1, &HYPRE_EuclidSetStats, plist_.get<int>("verbosity"));

  if (plist_.isParameter("ilu(k) fill level"))
    IfpHypre_->SetParameter(
      (Hypre_Chooser)1, &HYPRE_EuclidSetLevel, plist_.get<int>("ilu(k) fill level"));

  if (plist_.isParameter("rescale rows")) {
    bool rescale_rows = plist_.get<bool>("rescale rows");
    IfpHypre_->SetParameter((Hypre_Chooser)1, &HYPRE_EuclidSetRowScale, rescale_rows ? 1 : 0);
  }

  if (plist_.isParameter("ilut drop tolerance"))
    IfpHypre_->SetParameter(
      (Hypre_Chooser)1, &HYPRE_EuclidSetILUT, plist_.get<double>("ilut drop tolerance"));
#else
  Errors::Message msg("Hypre (Euclid) is not available in this installation of Amanzi.  To use "
                      "Hypre, please reconfigure.");
  Exceptions::amanzi_throw(msg);
#endif
}


/* ******************************************************************
* Initialize specialized AMS preconditioner.
****************************************************************** */
void
PreconditionerHypre::InitAMS_()
{
#ifdef HAVE_HYPRE
  method_ = AMS;

  if (plist_.isParameter("verbosity"))
    IfpHypre_->SetParameter(
      (Hypre_Chooser)1, &HYPRE_AMSSetPrintLevel, plist_.get<int>("verbosity"));

  // additional data: parallel vector with coordinates and gradient matrix
  // G_ * xyz_ should be the valid operation
  Teuchos::RCP<Epetra_MultiVector> xyz;
  Teuchos::RCP<Epetra_CrsMatrix> G;

  if (plist_.isParameter("graph coordinates") &&
      plist_.isType<Teuchos::RCP<Epetra_MultiVector>>("graph coordinates")) {
    xyz = plist_.get<Teuchos::RCP<Epetra_MultiVector>>("graph coordinates");
  }

  if (plist_.isParameter("discrete gradient operator") &&
      plist_.isType<Teuchos::RCP<Epetra_CrsMatrix>>("discrete gradient operator")) {
    G = plist_.get<Teuchos::RCP<Epetra_CrsMatrix>>("discrete gradient operator");
  }

  if (!xyz.get() || !G.get()) {
    Errors::Message msg("Hypre (AMS) needs additional data: gradient matrix and graph coordinates");
    Exceptions::amanzi_throw(msg);
  }

  // PList must go first
  Teuchos::ParameterList tmp;
  tmp.set<std::string>("hypre: Preconditioner", "AMS")
    .set<std::string>("hypre: SolveOrPrecondition", "Preconditioner");
  tmp.sublist("Coordinates").set<Teuchos::RCP<Epetra_MultiVector>>("Coordinates", xyz);
  tmp.sublist("Operators").set<Teuchos::RCP<const Epetra_CrsMatrix>>("G", G);
  tmp.sublist("HYPRE_AMSSetAlphaAMGOptions")
    .set<int>("arg 0", 6) // 10
    .set<int>("arg 1", 0) // 1
    .set<int>("arg 2", 6)
    .set<double>("arg 3", 0.25)
    .set<int>("arg 4", 0)
    .set<int>("arg 5", 0);
  tmp.sublist("HYPRE_AMSSetBetaAMGOptions")
    .set<int>("arg 0", 6) // 10
    .set<int>("arg 1", 0) // 1
    .set<int>("arg 2", 6)
    .set<double>("arg 3", 0.25)
    .set<int>("arg 4", 0)
    .set<int>("arg 5", 0);
  IfpHypre_->SetParameters(tmp);

  IfpHypre_->SetParameter((Hypre_Chooser)1, &HYPRE_AMSSetTol, plist_.get<double>("tolerance", 0.0));
  IfpHypre_->SetParameter(
    (Hypre_Chooser)1, &HYPRE_AMSSetPrintLevel, plist_.get<int>("verbosity", 0));
  IfpHypre_->SetParameter(
    (Hypre_Chooser)1, &HYPRE_AMSSetMaxIter, plist_.get<int>("cycle applications", 5));
  IfpHypre_->SetParameter(
    (Hypre_Chooser)1, &HYPRE_AMSSetCycleType, plist_.get<int>("cycle type", 1));
  IfpHypre_->SetParameter((Hypre_Chooser)1,
                          &HYPRE_AMSSetAlphaAMGCoarseRelaxType,
                          plist_.get<int>("coarse level relaxation", 9));
  IfpHypre_->SetParameter((Hypre_Chooser)1,
                          &HYPRE_AMSSetBetaAMGCoarseRelaxType,
                          plist_.get<int>("coarse level relaxation", 9));

#else
  Errors::Message msg("Hypre (AMS) is not available in this installation of Amanzi. To use Hypre, "
                      "please reconfigure.");
  Exceptions::amanzi_throw(msg);
#endif
}


/* ******************************************************************
* Amanzi preconditioner initialization
****************************************************************** */
void
PreconditionerHypre::InitializeInverse()
{
  IfpHypre_ = Teuchos::rcp(new Ifpack_Hypre(&*h_));

  // clear the list of functions
  // Teuchos::ParameterList plist;
  // IfpHypre_->SetParameters(plist);

  IfpHypre_->SetParameter((Hypre_Chooser)1);
  IfpHypre_->SetParameter(true);

  std::string method_name = plist_.get<std::string>("method");
  if (method_name == "boomer amg") {
    InitBoomer_();
  } else if (method_name == "euclid") {
    InitEuclid_();
  } else if (method_name == "ams") {
    InitAMS_();
  } else {
    Errors::Message msg;
    msg << "PreconditionerHypre: unknown method name \"" << method_name << "\"";
    Exceptions::amanzi_throw(msg);
  }

  // must reset the parameters every time to reset the block index
  IfpHypre_->SetParameter((Hypre_Chooser)1, method_);

  if (block_indices_.get()) {
    // must NEW the index array EVERY time, as it gets freed every time by
    // Hypre on cleanup.  This is pretty stupid, but we can't reuse it.  --etc
    //
    // Note this is not a memory leak -- this gets freed the second time
    // IfpHypre_::Compute() gets called (for every call but the last) and when
    // IfpHypre_ gets destroyed (for the last call).
    int* indices = new int[block_indices_->size()];
    for (int i = 0; i != block_indices_->size(); ++i) { indices[i] = (*block_indices_)[i]; }
    IfpHypre_->SetParameter((Hypre_Chooser)1, &HYPRE_BoomerAMGSetDofFunc, indices);
  }

  IfpHypre_->Initialize();
}


/* ******************************************************************
* Rebuild the preconditioner using the given matrix A.
****************************************************************** */
void
PreconditionerHypre::ComputeInverse()
{
#ifdef HAVE_HYPRE
  AMANZI_ASSERT(IfpHypre_.get());
  IfpHypre_->Compute();
#endif
}

} // namespace AmanziSolvers
} // namespace Amanzi
