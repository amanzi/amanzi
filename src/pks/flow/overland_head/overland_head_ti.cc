/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -----------------------------------------------------------------------------
This is the overland flow component of ATS.
License: BSD
Authors: Ethan Coon (ecoon@lanl.gov)
----------------------------------------------------------------------------- */

#include "Epetra_FECrsMatrix.h"
#include "EpetraExt_RowMatrixOut.h"
#include "boost/math/special_functions/fpclassify.hpp"

#include "overland_head.hh"
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
void OverlandHeadFlow::Functional( double t_old,
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
  vnames.push_back("p_old");
  vnames.push_back("p_new");
  vnames.push_back("z");
  vnames.push_back("h_old");
  vnames.push_back("h_new");
  vnames.push_back("h+z");

  std::vector< Teuchos::Ptr<const CompositeVector> > vecs;
  vecs.push_back(S_inter_->GetFieldData(key_).ptr());
  vecs.push_back(u.ptr());
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
  bc_seepage_pressure_->Compute(t_new);
  UpdateBoundaryConditions_(S_next_.ptr());

  // diffusion term, treated implicitly
  ApplyDiffusion_(S_next_.ptr(), res.ptr());

#if DEBUG_FLAG
  if (S_next_->HasField("unfrozen_fraction")) {
    vnames.resize(2);
    vecs.resize(2);
    vnames[0] = "uf_frac_old";
    vnames[1] = "uf_frac_new";
    vecs[0] = S_inter_->GetFieldData("unfrozen_fraction").ptr();
    vecs[1] = S_next_->GetFieldData("unfrozen_fraction").ptr();
    db_->WriteVectors(vnames, vecs, false);
  }
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

#if DEBUG_RES_FLAG
  if (niter_ < 23) {
    Teuchos::RCP<const CompositeVector> depth= S_next_->GetFieldData("ponded_depth");

    std::stringstream namestream;
    namestream << "flow_residual_" << niter_;
    *S_next_->GetFieldData(namestream.str(),name_) = *res;

    std::stringstream solnstream;
    solnstream << "flow_solution_" << niter_;
    *S_next_->GetFieldData(solnstream.str(),name_) = *depth;
  }
#endif

};


// -----------------------------------------------------------------------------
// Apply the preconditioner to u and return the result in Pu.
// -----------------------------------------------------------------------------
void OverlandHeadFlow::ApplyPreconditioner(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu) {
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_HIGH))
    *vo_->os() << "Precon application:" << std::endl;

#if DEBUG_FLAG
  db_->WriteVector("h_res", u->Data().ptr(), true);
#endif

  // apply the preconditioner
  lin_solver_->ApplyInverse(*u->Data(), *Pu->Data());

#if DEBUG_FLAG
  db_->WriteVector("PC*h_res (h-coords)", Pu->Data().ptr(), true);
#endif

  // tack on the variable change
  const Epetra_MultiVector& dh_dp =
      *S_next_->GetFieldData("dponded_depth_bar_d"+key_)->ViewComponent("cell",false);
  Epetra_MultiVector& Pu_c = *Pu->Data()->ViewComponent("cell",false);
  unsigned int ncells = Pu_c.MyLength();
  for (unsigned int c=0; c!=ncells; ++c) {
    Pu_c[0][c] /= dh_dp[0][c];
  }
#if DEBUG_FLAG
  db_->WriteVector("PC*h_res (p-coords)", Pu->Data().ptr(), true);
#endif
};


// -----------------------------------------------------------------------------
// Update the preconditioner at time t and u = up
// -----------------------------------------------------------------------------
void OverlandHeadFlow::UpdatePreconditioner(double t, Teuchos::RCP<const TreeVector> up, double h) {
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
  Teuchos::RCP<CompositeVector> flux = S_next_->GetFieldData("surface_flux", name_);
  preconditioner_diff_->UpdateFlux(*up->Data(), *flux);
  preconditioner_diff_->UpdateMatricesNewtonCorrection(flux.ptr(), Teuchos::null);

  // 2. Accumulation shift
  //    The desire is to keep this matrix invertible for pressures less than
  //    atmospheric.  To do that, we keep the accumulation derivative
  //    non-zero, calculating dWC_bar / dh_bar, where bar indicates (p -
  //    p_atm), not max(p - p_atm,0).  Note that this operator is in h
  //    coordinates, not p coordinates, as the diffusion operator is applied
  //    to h.
  //
  // -- update dh_bar / dp
  S_next_->GetFieldEvaluator("ponded_depth_bar")
      ->HasFieldDerivativeChanged(S_next_.ptr(), name_, key_);
  const Epetra_MultiVector& dh_dp =
      *S_next_->GetFieldData("dponded_depth_bar_d"+key_)
      ->ViewComponent("cell",false);
  const double& p_atm = *S_next_->GetScalarData("atmospheric_pressure");

  // -- update the accumulation derivatives
  S_next_->GetFieldEvaluator("surface_water_content_bar")
      ->HasFieldDerivativeChanged(S_next_.ptr(), name_, key_);
  const Epetra_MultiVector& dwc_dp =
      *S_next_->GetFieldData("dsurface_water_content_bar_d"+key_)
      ->ViewComponent("cell",false);

  db_->WriteVector("    dwc_dp", S_next_->GetFieldData("dsurface_water_content_bar_dsurface_pressure").ptr());
  db_->WriteVector("    dh_dp", S_next_->GetFieldData("dponded_depth_bar_dsurface_pressure").ptr());

  // -- pull out other needed data
  std::vector<double>& Acc_cells = preconditioner_acc_->local_matrices()->vals;
  unsigned int ncells = Acc_cells.size();
  for (unsigned int c=0; c!=ncells; ++c) {
    Acc_cells[c] += dwc_dp[0][c] / dh_dp[0][c] / h;
  }

  // // -- update the source term derivatives
  // if (S_next_->GetFieldEvaluator(mass_source_key_)->IsDependency(S_next_.ptr(), key_)) {
  //   S_next_->GetFieldEvaluator(mass_source_key_)
  //       ->HasFieldDerivativeChanged(S_next_.ptr(), name_, key_);
  //   std::string dkey = std::string("d")+mass_source_key_+std::string("_d")+key_;
  //   const Epetra_MultiVector& dq_dp = *S_next_->GetFieldData(dkey)
  //       ->ViewComponent("cell",false);

  //   const Epetra_MultiVector& cv =
  //       *S_next_->GetFieldData("surface_cell_volume")->ViewComponent("cell",false);
    
  //   if (source_in_meters_) {
  //     // External source term is in [m water / s], not in [mols / s], so a
  //     // density is required.  This density should be upwinded.
  //     S_next_->GetFieldEvaluator("surface_molar_density_liquid")
  //         ->HasFieldChanged(S_next_.ptr(), name_);
  //     S_next_->GetFieldEvaluator("surface_source_molar_density")
  //         ->HasFieldChanged(S_next_.ptr(), name_);
  //     const Epetra_MultiVector& nliq1 =
  //         *S_next_->GetFieldData("surface_molar_density_liquid")
  //         ->ViewComponent("cell",false);
  //     const Epetra_MultiVector& nliq1_s =
  //       *S_next_->GetFieldData("surface_source_molar_density")
  //         ->ViewComponent("cell",false);
  //     const Epetra_MultiVector& q = *S_next_->GetFieldData(mass_source_key_)
  //         ->ViewComponent("cell",false);

  //     for (int c=0; c!=cv.MyLength(); ++c) {
  //       double s1 = q[0][c] > 0. ? dq_dp[0][c] * nliq1_s[0][c] : dq_dp[0][c] * nliq1[0][c];
  //       Acc_cells[c] -= cv[0][c] * s1 / dh_dp[0][c];
  //     }
  //   } else {
  //     for (int c=0; c!=cv.MyLength(); ++c) {
  //       Acc_cells[c] -= cv[0][c] * dq_dp[0][c] / dh_dp[0][c];
  //     }
  //   }
  // }      
  

  // 3. Assemble and precompute the Schur complement for inversion.
  // 3.a: Patch up BCs in the case of zero conductivity
  FixBCsForPrecon_(S_next_.ptr());

  // // 3.c: Add in full Jacobian terms
  // if (tpfa_ && full_jacobian_) {
  //   if (vo_->os_OK(Teuchos::VERB_EXTREME))
  //     *vo_->os() << "    including full Jacobian terms" << std::endl;

  //   Teuchos::RCP<const CompositeVector> depth =
  //       S_next_->GetFieldData("ponded_depth");
  //   Teuchos::RCP<const CompositeVector> pres_elev =
  //       S_next_->GetFieldData("pres_elev");

  //   // Update conductivity.  Note the change of variables from pressure to
  //   // height.
  //   // -- Krel_face gets conductivity
  //   S_next_->GetFieldEvaluator("overland_conductivity")
  //       ->HasFieldDerivativeChanged(S_next_.ptr(), name_, key_);
  //   Teuchos::RCP<const CompositeVector> cond =
  //       S_next_->GetFieldData("overland_conductivity");
  //   Teuchos::RCP<const CompositeVector> dcond_dp =
  //       S_next_->GetFieldData("doverland_conductivity_dponded_depth");
  //   CompositeVector dcond_dh(*dcond_dp);
  //   dcond_dh.ViewComponent("cell",false)->ReciprocalMultiply(1., dh_dp,
  //           *dcond_dp->ViewComponent("cell",false), 0.);

  //   // -- Add in the Jacobian
  //   tpfa_preconditioner_->AnalyticJacobian(*upwinding_,
  //                                          S_next_.ptr(), "pres_elev",
  //                                          dcond_dh, bc_markers_,
  //                                          bc_values_);
  // }

  preconditioner_diff_->ApplyBCs(true, true);
  
  // 3.d: Rescale to use as a pressure matrix if used in a coupler
  if (coupled_to_subsurface_via_head_ || coupled_to_subsurface_via_flux_) {
    // Scale Spp by dh/dp (h, NOT h_bar), clobbering rows with p < p_atm
    std::string pd_key = smoothed_ponded_accumulation_ ? "smoothed_ponded_depth" : "ponded_depth";
    std::string pd_deriv_key = "d"+pd_key+"_d"+key_;
    S_next_->GetFieldEvaluator(pd_key)
        ->HasFieldDerivativeChanged(S_next_.ptr(), name_, key_);
    Teuchos::RCP<const CompositeVector> dh0_dp = S_next_->GetFieldData(pd_deriv_key);
    const Epetra_MultiVector& dh0_dp_c = *dh0_dp->ViewComponent("cell",false);
    
    preconditioner_->Rescale(*dh0_dp);
    
    if (vo_->os_OK(Teuchos::VERB_EXTREME))
      *vo_->os() << "  Right scaling TPFA" << std::endl;
    db_->WriteVector("    dh_dp", dh0_dp.ptr());
  }


  if (precon_used_) {
    preconditioner_->AssembleMatrix();
    preconditioner_->InitPreconditioner(plist_->sublist("preconditioner"));
  }      
  
  /*
  // dump the schur complement
  Teuchos::RCP<Epetra_FECrsMatrix> sc = mfd_preconditioner_->Schur();
  std::stringstream filename_s;
  filename_s << "schur_" << S_next_->cycle() << ".txt";
  EpetraExt::RowMatrixToMatlabFile(filename_s.str().c_str(), *sc);
  *vo_->os() << "updated precon " << S_next_->cycle() << std::endl;
  */
};

}  // namespace Flow
}  // namespace Amanzi
