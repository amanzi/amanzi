/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -----------------------------------------------------------------------------
This is the overland flow component of ATS.
License: BSD
Authors: Ethan Coon (ecoon@lanl.gov)
----------------------------------------------------------------------------- */

#include "Epetra_FECrsMatrix.h"
#include "EpetraExt_RowMatrixOut.h"
#include "boost/math/special_functions/fpclassify.hpp"
#include "MatrixMFD_TPFA.hh"

#include "snow_distribution.hh"

#define DEBUG_FLAG 1

namespace Amanzi {
namespace Flow {

// Overland is a BDFFnBase
// -----------------------------------------------------------------------------
// computes the non-linear functional g = g(t,u,udot)
// -----------------------------------------------------------------------------
void SnowDistribution::Functional( double t_old,
                        double t_new,
                        Teuchos::RCP<TreeVector> u_old,
                        Teuchos::RCP<TreeVector> u_new,
                        Teuchos::RCP<TreeVector> g ) {
  // VerboseObject stuff.
  Teuchos::OSTab tab = vo_->getOSTab();

  // bookkeeping
  double h = t_new - t_old;
  ASSERT(std::abs(S_inter_->time() - t_old) < 1.e-4*h);
  ASSERT(std::abs(S_next_->time() - t_new) < 1.e-4*h);

  Teuchos::RCP<CompositeVector> u = u_new->Data();

  // zero out residual
  Teuchos::RCP<CompositeVector> res = g->Data();
  res->PutScalar(0.0);

#if DEBUG_FLAG
  if (vo_->os_OK(Teuchos::VERB_HIGH))
    *vo_->os() << "----------------------------------------------------------------" << std::endl
               << "Residual calculation: t0 = " << t_old
               << " t1 = " << t_new << " h = " << h << std::endl;
#endif

  // unnecessary here if not debeugging, but doesn't hurt either
  S_next_->GetFieldEvaluator("snow_skin_potential")->HasFieldChanged(S_next_.ptr(), name_);

#if DEBUG_FLAG
  // dump u_old, u_new
  db_->WriteCellInfo(true);
  std::vector<std::string> vnames;
  vnames.push_back("precip_new");
  vnames.push_back("potential");

  std::vector< Teuchos::Ptr<const CompositeVector> > vecs;
  vecs.push_back(u.ptr());
  vecs.push_back(S_next_->GetFieldData("snow_skin_potential").ptr());

  db_->WriteVectors(vnames, vecs, true);
#endif

  // pointer-copy temperature into state and update any auxilary data
  solution_to_state(u_new, S_next_);

  // diffusion term, treated implicitly
  ApplyDiffusion_(S_next_.ptr(), res.ptr());

#if DEBUG_FLAG
  db_->WriteVector("k_s", S_next_->GetFieldData("upwind_snow_conductivity").ptr(), true);
  db_->WriteVector("res (post diffusion)", res.ptr(), true);
#endif

  // accumulation term
  AddAccumulation_(res.ptr());
#if DEBUG_FLAG
  if (vo_->os_OK(Teuchos::VERB_HIGH)) {
    std::vector<double> time(1,S_next_->time());
    double precip = (*precip_func_)(time);
    *vo_->os() << " precip = " << precip << std::endl;
  }
  db_->WriteVector("res (post accumulation)", res.ptr(), true);
#endif
};


// -----------------------------------------------------------------------------
// Apply the preconditioner to u and return the result in Pu.
// -----------------------------------------------------------------------------
void SnowDistribution::ApplyPreconditioner(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu) {
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_HIGH))
    *vo_->os() << "Precon application:" << std::endl;

#if DEBUG_FLAG
  db_->WriteVector("h_res", u->Data().ptr(), true);
#endif

  // apply the preconditioner
  mfd_preconditioner_->ApplyInverse(*u->Data(), *Pu->Data());

#if DEBUG_FLAG
  db_->WriteVector("PC*h_res", Pu->Data().ptr(), true);
#endif
};


// -----------------------------------------------------------------------------
// Update the preconditioner at time t and u = up
// -----------------------------------------------------------------------------
void SnowDistribution::UpdatePreconditioner(double t, Teuchos::RCP<const TreeVector> up, double h) {
  // VerboseObject stuff.
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_EXTREME))
    *vo_->os() << "Precon update at t = " << t << std::endl;

  // update state with the solution up.
  ASSERT(std::abs(S_next_->time() - t) <= 1.e-4*t);
  PKDefaultBase::solution_to_state(up, S_next_);

  // update the rel perm according to the scheme of choice
  UpdatePermeabilityData_(S_next_.ptr());

  Teuchos::RCP<const CompositeVector> cond =
    S_next_->GetFieldData("upwind_snow_conductivity");

  // calculating the operator is done in 3 steps:
  // 1. Create all local matrices.
  mfd_preconditioner_->CreateMFDstiffnessMatrices(cond.ptr());
  mfd_preconditioner_->CreateMFDrhsVectors();

  // 2.b Update local matrices diagonal with the accumulation terms.
  // -- update the accumulation derivatives
  const Epetra_MultiVector& cell_volume = *S_next_->GetFieldData("surface_cell_volume")
      ->ViewComponent("cell",false);

  std::vector<double>& Acc_cells = mfd_preconditioner_->Acc_cells();
  int ncells = cell_volume.MyLength();
  for (int c=0; c!=ncells; ++c) {
    // accumulation term
    Acc_cells[c] += cell_volume[0][c] * 10.;
  }
};

double SnowDistribution::ErrorNorm(Teuchos::RCP<const TreeVector> u,
                       Teuchos::RCP<const TreeVector> du) {
  Teuchos::OSTab tab = vo_->getOSTab();

  Teuchos::RCP<const CompositeVector> res = du->Data();
  const Epetra_MultiVector& res_c = *res->ViewComponent("cell",false);
  const Epetra_MultiVector& precip_c = *u->Data()->ViewComponent("cell",false);
  const Epetra_MultiVector& cv = *S_next_->GetFieldData("surface_cell_volume")
      ->ViewComponent("cell",false);
  double dt = S_next_->time() - S_inter_->time();
  std::vector<double> time(1, S_next_->time());
  
  // Cell error is based upon error in mass conservation
  double Qe = (*precip_func_)(time);
  double enorm_cell(0.);
  int bad_cell = -1;
  unsigned int ncells = res_c.MyLength();
  for (unsigned int c=0; c!=ncells; ++c) {
    double tmp = std::abs(res_c[0][c]*dt)
        / (atol_ * .1 * cv[0][c] + rtol_ * 10 * dt * precip_c[0][c]);
    if (tmp > enorm_cell) {
      enorm_cell = tmp;
      bad_cell = c;
    }
  }

  // Write out Inf norms too.
  if (vo_->os_OK(Teuchos::VERB_MEDIUM)) {
    double infnorm_c(0.);
    res_c.NormInf(&infnorm_c);

    ENorm_t err_c;
#ifdef HAVE_MPI
    ENorm_t l_err_c;
    l_err_c.value = enorm_cell;
    l_err_c.gid = res_c.Map().GID(bad_cell);

    MPI_Allreduce(&l_err_c, &err_c, 1, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD);
#else
    err_c.value = enorm_cell;
    err_c.gid = bad_cell;
#endif

    *vo_->os() << "ENorm (cells) = " << err_c.value << "[" << err_c.gid << "] (" << infnorm_c << ")" << std::endl;
  }

#ifdef HAVE_MPI
  double buf = enorm_cell;
  MPI_Allreduce(&buf, &enorm_cell, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#endif
  return enorm_cell;
};

}  // namespace Flow
}  // namespace Amanzi
