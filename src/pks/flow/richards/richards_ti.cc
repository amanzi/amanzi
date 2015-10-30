/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
A base two-phase, thermal Richard's equation with water vapor.

Authors: Ethan Coon (ATS version) (ecoon@lanl.gov)
*/

#include "Epetra_FECrsMatrix.h"
#include "EpetraExt_RowMatrixOut.h"
#include "boost/math/special_functions/fpclassify.hpp"

#include "Op.hh"
#include "richards.hh"

namespace Amanzi {
namespace Flow {

#define DEBUG_FLAG 1
#define DEBUG_RES_FLAG 0

// Richards is a BDFFnBase
// -----------------------------------------------------------------------------
// computes the non-linear functional g = g(t,u,udot)
// -----------------------------------------------------------------------------
void Richards::Functional(double t_old,
                   double t_new,
                   Teuchos::RCP<TreeVector> u_old,
                   Teuchos::RCP<TreeVector> u_new,
                   Teuchos::RCP<TreeVector> g) {
  // VerboseObject stuff.
  Teuchos::OSTab tab = vo_->getOSTab();

  double h = t_new - t_old;
  ASSERT(std::abs(S_inter_->time() - t_old) < 1.e-4*h);
  ASSERT(std::abs(S_next_->time() - t_new) < 1.e-4*h);

  // pointer-copy temperature into state and update any auxilary data
  solution_to_state(*u_new, S_next_);
  Teuchos::RCP<CompositeVector> u = u_new->Data();

  if (dynamic_mesh_) matrix_diff_->SetTensorCoefficient(K_);

#if DEBUG_FLAG
  if (vo_->os_OK(Teuchos::VERB_HIGH))
    *vo_->os() << "----------------------------------------------------------------" << std::endl
               << "Residual calculation: t0 = " << t_old
               << " t1 = " << t_new << " h = " << h << std::endl;

  // dump u_old, u_new
  db_->WriteCellInfo(true);
  std::vector<std::string> vnames;
  vnames.push_back("p_old"); vnames.push_back("p_new");
  std::vector< Teuchos::Ptr<const CompositeVector> > vecs;
  vecs.push_back(S_inter_->GetFieldData(key_).ptr()); vecs.push_back(u.ptr());
  db_->WriteVectors(vnames, vecs, true);
#endif

  // update boundary conditions
  bc_pressure_->Compute(t_new);
  bc_flux_->Compute(t_new);
  UpdateBoundaryConditions_(S_next_.ptr());

  // zero out residual
  Teuchos::RCP<CompositeVector> res = g->Data();
  res->PutScalar(0.0);

  // diffusion term, treated implicitly
  ApplyDiffusion_(S_next_.ptr(), res.ptr());

  
  // if (vapor_diffusion_) AddVaporDiffusionResidual_(S_next_.ptr(), res.ptr());



#if DEBUG_FLAG
  // dump s_old, s_new
  vnames[0] = "sl_old"; vnames[1] = "sl_new";
  vecs[0] = S_inter_->GetFieldData(getKey(domain_,"saturation_liquid")).ptr();
  vecs[1] = S_next_->GetFieldData(getKey(domain_,"saturation_liquid")).ptr();

  if (S_next_->HasField("saturation_ice")) {
    vnames.push_back("si_old");
    vnames.push_back("si_new");
    vecs.push_back(S_inter_->GetFieldData(getKey(domain_,"saturation_ice")).ptr());
    vecs.push_back(S_next_->GetFieldData(getKey(domain_,"saturation_ice")).ptr());
  }

  vnames.push_back("k_rel");
  vecs.push_back(S_next_->GetFieldData(coef_key_).ptr());
  vnames.push_back("wind");
  vecs.push_back(S_next_->GetFieldData(flux_dir_key_).ptr());
  vnames.push_back("uw_k_rel");
  vecs.push_back(S_next_->GetFieldData(uw_coef_key_).ptr());
  db_->WriteVectors(vnames,vecs,true);

  db_->WriteVector("res (diff)", res.ptr(), true);
#endif

  // accumulation term
  AddAccumulation_(res.ptr());

  // source term
  if (is_source_term_) {
    if (explicit_source_) {
      AddSources_(S_inter_.ptr(), res.ptr());
    } else {
      AddSources_(S_next_.ptr(), res.ptr());
    }
  }
};

// -----------------------------------------------------------------------------
// Apply the preconditioner to u and return the result in Pu.
// -----------------------------------------------------------------------------
int Richards::ApplyPreconditioner(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu) {
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_HIGH))
    *vo_->os() << "Precon application:" << std::endl;

#if DEBUG_FLAG
  db_->WriteVector("p_res", u->Data().ptr(), true);
#endif

  // Apply the preconditioner
  int ierr = lin_solver_->ApplyInverse(*u->Data(), *Pu->Data());

#if DEBUG_FLAG
  db_->WriteVector("PC*p_res", Pu->Data().ptr(), true);
#endif
  
  return (ierr > 0) ? 0 : 1;
};


// -----------------------------------------------------------------------------
// Update the preconditioner at time t and u = up
// -----------------------------------------------------------------------------
void Richards::UpdatePreconditioner(double t, Teuchos::RCP<const TreeVector> up, double h) {
  // VerboseObject stuff.
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_HIGH))
    *vo_->os() << "Precon update at t = " << t << std::endl;

  // Recreate mass matrices
  if (dynamic_mesh_) {
    matrix_diff_->SetTensorCoefficient(K_);
    preconditioner_diff_->SetTensorCoefficient(K_);
  }

  // update state with the solution up.
  ASSERT(std::abs(S_next_->time() - t) <= 1.e-4*t);
  PKDefaultBase::solution_to_state(*up, S_next_);

  // update the rel perm according to the scheme of choice, also upwind derivatives of rel perm
  UpdatePermeabilityData_(S_next_.ptr());
  if (jacobian_) UpdatePermeabilityDerivativeData_(S_next_.ptr());

  // update boundary conditions
  bc_pressure_->Compute(S_next_->time());
  bc_flux_->Compute(S_next_->time());
  UpdateBoundaryConditions_(S_next_.ptr());

  Teuchos::RCP<const CompositeVector> rel_perm =
      S_next_->GetFieldData(uw_coef_key_);

  // Update the preconditioner with darcy and gravity fluxes
  preconditioner_->Init();

  // gravity fluxes
  S_next_->GetFieldEvaluator(mass_dens_key_)->HasFieldChanged(S_next_.ptr(), name_);
  Teuchos::RCP<const CompositeVector> rho = S_next_->GetFieldData(mass_dens_key_);
  preconditioner_diff_->SetDensity(rho);

  // jacobian term
  Teuchos::RCP<const CompositeVector> dkrdp = Teuchos::null;
  if (jacobian_) {
    if (!duw_coef_key_.empty()) {
      dkrdp = S_next_->GetFieldData(duw_coef_key_);
    } else {
      dkrdp = S_next_->GetFieldData(dcoef_key_);
    }
  }

  // create local matrices
  preconditioner_diff_->SetScalarCoefficient(rel_perm, dkrdp);
  preconditioner_diff_->UpdateMatrices(Teuchos::null, up->Data().ptr());

  if (jacobian_ && preconditioner_->RangeMap().HasComponent("face")) {
    Teuchos::RCP<CompositeVector> flux = S_next_->GetFieldData(flux_key_, name_);
    preconditioner_diff_->UpdateFlux(*up->Data(), *flux);
    preconditioner_diff_->UpdateMatricesNewtonCorrection(flux.ptr(), up->Data().ptr());
  }
  
  // if (vapor_diffusion_){
  //   Teuchos::RCP<CompositeVector> vapor_diff_pres = S_next_->GetFieldData("vapor_diffusion_pressure", name_);
  //   ComputeVaporDiffusionCoef(S_next_.ptr(), vapor_diff_pres, "pressure");   

  // // update the stiffness matrix
    // matrix_vapor_->CreateMFDstiffnessMatrices(vapor_diff_pres.ptr());    
    // mfd_preconditioner_->Add2MFDstiffnessMatrices(&matrix_vapor_->Acc_cells(),
    //                                               &matrix_vapor_->Aff_cells(),
    //                                               &matrix_vapor_->Acf_cells(),
    //                                               &matrix_vapor_->Afc_cells());
  // }

  // Update the preconditioner with accumulation terms.
  // -- update the accumulation derivatives
  S_next_->GetFieldEvaluator(conserved_key_)
      ->HasFieldDerivativeChanged(S_next_.ptr(), name_, key_);

  // -- get the accumulation deriv
  Key dwc_dp_key = getDerivKey(conserved_key_, key_);
  const Epetra_MultiVector& dwc_dp =
      *S_next_->GetFieldData(dwc_dp_key)->ViewComponent("cell",false);
  const Epetra_MultiVector& pres =
      *S_next_->GetFieldData(key_)->ViewComponent("cell",false);

#if DEBUG_FLAG
  db_->WriteVector("    dwc_dp", S_next_->GetFieldData(dwc_dp_key).ptr());
#endif

  // -- update the cell-cell block
  std::vector<double>& Acc_cells = preconditioner_acc_->local_matrices()->vals;
  unsigned int ncells = dwc_dp.MyLength();
  for (unsigned int c=0; c!=ncells; ++c) {
    //    ASSERT(dwc_dp[0][c] > 1.e-10);
    Acc_cells[c] += dwc_dp[0][c] / h;
  }

  // -- update preconditioner with source term derivatives if needed
  AddSourcesToPrecon_(S_next_.ptr(), h);
  
  // -- apply BCs
  preconditioner_diff_->ApplyBCs(true, true);

  if (precon_used_) {
    preconditioner_->AssembleMatrix();
    preconditioner_->InitPreconditioner(plist_->sublist("preconditioner"));
  }      
};


}  // namespace Flow
}  // namespace Amanzi



