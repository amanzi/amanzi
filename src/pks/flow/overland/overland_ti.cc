/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -----------------------------------------------------------------------------
This is the overland flow component of ATS.
License: BSD
Authors: Ethan Coon (ecoon@lanl.gov)
----------------------------------------------------------------------------- */

#include "Epetra_FECrsMatrix.h"
#include "EpetraExt_RowMatrixOut.h"
#include "boost/math/special_functions/fpclassify.hpp"

#include "overland.hh"
#include "Op.hh"

namespace Amanzi {
namespace Flow {

#define DEBUG_FLAG 1
#define DEBUG_ICE_FLAG 0
#define DEBUG_RES_FLAG 0


// Overland is a BDFFnBase
// -----------------------------------------------------------------------------
// computes the non-linear functional g = g(t,u,udot)
// -----------------------------------------------------------------------------
void OverlandFlow::Functional( double t_old,
                        double t_new,
                        Teuchos::RCP<TreeVector> u_old,
                        Teuchos::RCP<TreeVector> u_new,
                        Teuchos::RCP<TreeVector> g ) {
  // VerboseObject stuff.
  Teuchos::OSTab tab = vo_->getOSTab();
  niter_++;

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
  S_next_->GetFieldEvaluator("pres_elev")->HasFieldChanged(S_next_.ptr(), name_);

#if DEBUG_FLAG
  // dump u_old, u_new
  db_->WriteCellInfo(true);
  std::vector<std::string> vnames;
  vnames.push_back("z");
  vnames.push_back("h_old");
  vnames.push_back("h_new");
  vnames.push_back("h+z");

  std::vector< Teuchos::Ptr<const CompositeVector> > vecs;
  vecs.push_back(S_inter_->GetFieldData("elevation").ptr());
  vecs.push_back(S_inter_->GetFieldData("ponded_depth").ptr());
  vecs.push_back(S_next_->GetFieldData("ponded_depth").ptr());
  vecs.push_back(S_next_->GetFieldData("pres_elev").ptr());

  db_->WriteVectors(vnames, vecs, true);
#endif

  // pointer-copy temperature into state and update any auxilary data
  solution_to_state(*u_new, S_next_);

  // update boundary conditions
  bc_head_->Compute(t_new);
  bc_flux_->Compute(t_new);
  bc_seepage_head_->Compute(t_new);
  UpdateBoundaryConditions_(S_next_.ptr());

  // diffusion term, treated implicitly
  ApplyDiffusion_(S_next_.ptr(), res.ptr());

#if DEBUG_FLAG
  db_->WriteVector("k_s", S_next_->GetFieldData("upwind_overland_conductivity").ptr(), true);
  db_->WriteVector("res (diff)", res.ptr(), true);
#endif

  // accumulation term
  AddAccumulation_(res.ptr());
#if DEBUG_FLAG
  db_->WriteVector("res (acc)", res.ptr(), true);
#endif

  // add rhs load value
  AddSourceTerms_(res.ptr());
#if DEBUG_FLAG
  db_->WriteVector("res (src)", res.ptr(), true);
#endif

};


// -----------------------------------------------------------------------------
// Apply the preconditioner to u and return the result in Pu.
// -----------------------------------------------------------------------------
void OverlandFlow::ApplyPreconditioner(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu) {
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_HIGH))
    *vo_->os() << "Precon application:" << std::endl;

#if DEBUG_FLAG
  db_->WriteVector("h_res", u->Data().ptr(), true);
#endif

  // apply the preconditioner
  preconditioner_->ApplyInverse(*u->Data(), *Pu->Data());

#if DEBUG_FLAG
  db_->WriteVector("PC*h_res (h-coords)", Pu->Data().ptr(), true);
#endif

};


// -----------------------------------------------------------------------------
// Update the preconditioner at time t and u = up
// -----------------------------------------------------------------------------
void OverlandFlow::UpdatePreconditioner(double t, Teuchos::RCP<const TreeVector> up, double h) {
  // VerboseObject stuff.
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_EXTREME))
    *vo_->os() << "Precon update at t = " << t << std::endl;

  // update state with the solution up.
  ASSERT(std::abs(S_next_->time() - t) <= 1.e-4*t);
  PKDefaultBase::solution_to_state(*up, S_next_);

  // calculating the operator is done in 3 steps:
  // 1. Diffusion components

  // 1.a: Pre-assembly updates.
  // -- update boundary condition markers, which set the BC type
  bc_head_->Compute(t);
  bc_flux_->Compute(t);
  bc_seepage_head_->Compute(t);
  UpdateBoundaryConditions_(S_next_.ptr());

  // -- update the rel perm according to the boundary info and upwinding
  // -- scheme of choice
  UpdatePermeabilityData_(S_next_.ptr());
  UpdatePermeabilityDerivativeData_(S_next_.ptr());

  Teuchos::RCP<const CompositeVector> cond =
    S_next_->GetFieldData("upwind_overland_conductivity");
  Teuchos::RCP<const CompositeVector> dcond =
    S_next_->GetFieldData("dupwind_overland_conductivity_dponded_depth");

  // 1.b: Create all local matrices.
  preconditioner_->Init();
  preconditioner_diff_->Setup(cond, Teuchos::null);
  Teuchos::RCP<const CompositeVector> pres_elev = S_next_->GetFieldData("pres_elev");
  preconditioner_diff_->UpdateMatrices(Teuchos::null, pres_elev.ptr());


  // 2. Accumulation shift
  //    The desire is to keep this matrix invertible for pressures less than
  //    atmospheric.  To do that, we keep the accumulation derivative
  //    non-zero, calculating dWC_bar / dh_bar, where bar indicates (p -
  //    p_atm), not max(p - p_atm,0).  Note that this operator is in h
  //    coordinates, not p coordinates, as the diffusion operator is applied
  //    to h.
  //

  // -- update the accumulation derivatives
  const Epetra_MultiVector& cv =
      *S_next_->GetFieldData("surface_cell_volume")->ViewComponent("cell",false);

  std::vector<double>& Acc_cells = preconditioner_acc_->local_matrices()->vals;
  unsigned int ncells = Acc_cells.size();
  for (unsigned int c=0; c!=ncells; ++c) {
    Acc_cells[c] += cv[0][c] / h;
  }

  preconditioner_diff_->ApplyBCs(true, true);
  preconditioner_->AssembleMatrix();
  preconditioner_->InitPreconditioner(plist_->sublist("preconditioner"));
};

}  // namespace Flow
}  // namespace Amanzi
