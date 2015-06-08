/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
A base two-phase, thermal Richard's equation with water vapor.

Authors: Ethan Coon (ATS version) (ecoon@lanl.gov)
*/

#include "Epetra_FECrsMatrix.h"
#include "EpetraExt_RowMatrixOut.h"
#include "boost/math/special_functions/fpclassify.hpp"
#include "Mesh_MSTK.hh"

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
  niter_++;

  double h = t_new - t_old;
  ASSERT(std::abs(S_inter_->time() - t_old) < 1.e-4*h);
  ASSERT(std::abs(S_next_->time() - t_new) < 1.e-4*h);

  // pointer-copy temperature into state and update any auxilary data
  solution_to_state(*u_new, S_next_);
  Teuchos::RCP<CompositeVector> u = u_new->Data();

  if (dynamic_mesh_) matrix_diff_->Setup(K_);

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
  UpdateBoundaryConditions_();

  // zero out residual
  Teuchos::RCP<CompositeVector> res = g->Data();
  res->PutScalar(0.0);

  // diffusion term, treated implicitly
  ApplyDiffusion_(S_next_.ptr(), res.ptr());

  
  // if (vapor_diffusion_) AddVaporDiffusionResidual_(S_next_.ptr(), res.ptr());



#if DEBUG_FLAG
  // dump s_old, s_new
  vnames[0] = "sl_old"; vnames[1] = "sl_new";
  vecs[0] = S_inter_->GetFieldData("saturation_liquid").ptr();
  vecs[1] = S_next_->GetFieldData("saturation_liquid").ptr();

  if (S_next_->HasField("saturation_ice")) {
    vnames.push_back("si_old");
    vnames.push_back("si_new");
    vecs.push_back(S_inter_->GetFieldData("saturation_ice").ptr());
    vecs.push_back(S_next_->GetFieldData("saturation_ice").ptr());
  }

  vnames.push_back("k_rel");
  vecs.push_back(S_next_->GetFieldData("relative_permeability").ptr());
  vnames.push_back("wind");
  vecs.push_back(S_next_->GetFieldData("darcy_flux_direction").ptr());
  vnames.push_back("uw_k_rel");
  vecs.push_back(S_next_->GetFieldData("numerical_rel_perm").ptr());
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
  
#if DEBUG_RES_FLAG
  if (niter_ < 23) {
    std::stringstream namestream;
    namestream << "flow_residual_" << niter_;
    *S_next_->GetFieldData(namestream.str(),name_) = *res;

    std::stringstream solnstream;
    solnstream << "flow_solution_" << niter_;
    *S_next_->GetFieldData(solnstream.str(),name_) = *u;
  }
#endif
};

// -----------------------------------------------------------------------------
// Apply the preconditioner to u and return the result in Pu.
// -----------------------------------------------------------------------------
void Richards::ApplyPreconditioner(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu) {
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_HIGH))
    *vo_->os() << "Precon application:" << std::endl;

#if DEBUG_FLAG
  db_->WriteVector("p_res", u->Data().ptr(), true);
#endif

  // Apply the preconditioner
  preconditioner_->ApplyInverse(*u->Data(), *Pu->Data());

#if DEBUG_FLAG
  db_->WriteVector("PC*p_res", Pu->Data().ptr(), true);
#endif

//   if (precon_wc_) {
//     PreconWC_(u, Pu);
// #if DEBUG_FLAG
//     db_->WriteVector("PC_WC*p_res", Pu->Data().ptr(), true);
// #endif
//   }
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
    matrix_diff_->Setup(K_);
    preconditioner_diff_->Setup(K_);
  }

  // update state with the solution up.
  ASSERT(std::abs(S_next_->time() - t) <= 1.e-4*t);
  PKDefaultBase::solution_to_state(*up, S_next_);

  // update the rel perm according to the scheme of choice, also upwind derivatives of rel perm
  UpdatePermeabilityData_(S_next_.ptr());
  UpdatePermeabilityDerivativeData_(S_next_.ptr());

  // update boundary conditions
  bc_pressure_->Compute(S_next_->time());
  bc_flux_->Compute(S_next_->time());
  UpdateBoundaryConditions_();

  Teuchos::RCP<const CompositeVector> rel_perm =
      S_next_->GetFieldData("numerical_rel_perm");

  if (vo_->os_OK(Teuchos::VERB_HIGH)) {
    const Epetra_MultiVector& kr = *rel_perm->ViewComponent("face",false);
    double min_kr = 1.e6;
    int min_kr_lid = -1;
    for (int f=0; f!=kr.MyLength(); ++f) {
      if (kr[0][f] < min_kr) {
        min_kr = kr[0][f];
        min_kr_lid = f;
      }
    }
    ASSERT(min_kr_lid >= 0);

    ENorm_t global_min_kr;
    ENorm_t local_min_kr;
    local_min_kr.value = min_kr;
    local_min_kr.gid = kr.Map().GID(min_kr_lid);
    ASSERT(local_min_kr.gid >= 0);

    MPI_Allreduce(&local_min_kr, &global_min_kr, 1, MPI_DOUBLE_INT, MPI_MINLOC, MPI_COMM_WORLD);

    *vo_->os() << "Min Kr[face=" << global_min_kr.gid << "] = " << global_min_kr.value << std::endl;
  }

  Teuchos::RCP<CompositeVector> rel_perm_modified =
      Teuchos::rcp(new CompositeVector(*rel_perm));
  *rel_perm_modified = *rel_perm;

  {
    Epetra_MultiVector& rel_perm_mod_f = *rel_perm_modified->ViewComponent("face",false);
    unsigned int nfaces = rel_perm_mod_f.MyLength();
    for (unsigned int f=0; f!=nfaces; ++f) {
      rel_perm_mod_f[0][f] = std::max(rel_perm_mod_f[0][f], 1.e-18);
    }
  }      
  
  // Update the preconditioner with darcy and gravity fluxes
  preconditioner_->Init();

  S_next_->GetFieldEvaluator("mass_density_liquid")->HasFieldChanged(S_next_.ptr(), name_);
  Teuchos::RCP<const CompositeVector> rho = S_next_->GetFieldData("mass_density_liquid");
  preconditioner_diff_->SetDensity(rho);
  
  Teuchos::RCP<const CompositeVector> dkrdp = S_next_->GetFieldData("dnumerical_rel_perm_dpressure");
  preconditioner_diff_->Setup(rel_perm_modified, dkrdp);
  // Teuchos::RCP<const CompositeVector> flux = S_next_->GetFieldData("darcy_flux");
  // preconditioner_diff_->UpdateMatrices(flux.ptr(), up->Data().ptr());
  preconditioner_diff_->UpdateMatrices(Teuchos::null, Teuchos::null);
  Teuchos::RCP<CompositeVector> flux = S_next_->GetFieldData("darcy_flux", name_);
  preconditioner_diff_->UpdateFlux(*up->Data(), *flux);
  preconditioner_diff_->UpdateMatricesNewtonCorrection(flux.ptr(), Teuchos::null);
  

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
  S_next_->GetFieldEvaluator("water_content")
      ->HasFieldDerivativeChanged(S_next_.ptr(), name_, key_);

  // -- get the accumulation deriv
  const Epetra_MultiVector& dwc_dp =
      *S_next_->GetFieldData("dwater_content_d"+key_)->ViewComponent("cell",false);
  const Epetra_MultiVector& pres =
      *S_next_->GetFieldData(key_)->ViewComponent("cell",false);

#if DEBUG_FLAG
  db_->WriteVector("    dwc_dp", S_next_->GetFieldData("dwater_content_d"+key_).ptr());
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



