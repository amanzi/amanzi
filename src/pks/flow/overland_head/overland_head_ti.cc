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

#include "overland_head.hh"

namespace Amanzi {
namespace Flow {

#define DEBUG_FLAG 1
#define DEBUG_ICE_FLAG 0
#define DEBUG_RES_FLAG 0


// Overland is a BDFFnBase
// -----------------------------------------------------------------------------
// computes the non-linear functional g = g(t,u,udot)
// -----------------------------------------------------------------------------
void OverlandHeadFlow::fun( double t_old,
                        double t_new,
                        Teuchos::RCP<TreeVector> u_old,
                        Teuchos::RCP<TreeVector> u_new,
                        Teuchos::RCP<TreeVector> g ) {
  // VerboseObject stuff.
  Teuchos::OSTab tab = getOSTab();
  niter_++;

  // bookkeeping
  double h = t_new - t_old;
  ASSERT(std::abs(S_inter_->time() - t_old) < 1.e-4*h);
  ASSERT(std::abs(S_next_->time() - t_new) < 1.e-4*h);

  Teuchos::RCP<CompositeVector> u = u_new->data();

  // zero out residual
  Teuchos::RCP<CompositeVector> res = g->data();
  res->PutScalar(0.0);

#if DEBUG_FLAG
  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_HIGH, true)) {
    *out_ << "----------------------------------------------------------------" << std::endl;

    *out_ << "Residual calculation: T0 = " << t_old
          << " T1 = " << t_new << " H = " << h << std::endl;
    *out_ << std::setprecision(15);

    S_next_->GetFieldEvaluator("pres_elev")->HasFieldChanged(S_next_.ptr(), name_);
    Teuchos::RCP<const CompositeVector> depth= S_next_->GetFieldData("ponded_depth");
    Teuchos::RCP<const CompositeVector> preselev= S_next_->GetFieldData("pres_elev");

    unsigned int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
    for (std::vector<AmanziMesh::Entity_ID>::const_iterator c0=dc_.begin(); c0!=dc_.end(); ++c0) {
      if (*c0 < ncells) {
        AmanziMesh::Entity_ID_List fnums0;
        std::vector<int> dirs;
        mesh_->cell_get_faces_and_dirs(*c0, &fnums0, &dirs);

        *out_ << "  p(" << *c0 <<"): " << (*u)("cell",0,*c0) << std::endl;
        *out_ << "  h(" << *c0 <<"): " << (*depth)("cell",0,*c0);
        for (unsigned int n=0; n!=fnums0.size(); ++n) *out_ << ",  " << (*depth)("face",fnums0[n]);
        *out_ << std::endl;
        *out_ << "  hz(" << *c0 <<"): " << (*preselev)("cell",0,*c0);
        for (unsigned int n=0; n!=fnums0.size(); ++n) *out_ << ",  " << (*preselev)("face",fnums0[n]);
        *out_ << std::endl;
        *out_ << " --" << std::endl;
      }
    }
  }
#endif

  // pointer-copy temperature into state and update any auxilary data
  solution_to_state(u_new, S_next_);

  // update boundary conditions
  bc_pressure_->Compute(t_new);
  bc_head_->Compute(t_new);
  bc_flux_->Compute(t_new);
  UpdateBoundaryConditions_(S_next_.ptr());

  // diffusion term, treated implicitly
  ApplyDiffusion_(S_next_.ptr(), res.ptr());

#if DEBUG_FLAG
  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_HIGH, true)) {
    Teuchos::RCP<const CompositeVector> cond =
      S_next_->GetFieldData("upwind_overland_conductivity", name_);
    unsigned int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
    for (std::vector<AmanziMesh::Entity_ID>::const_iterator c0=dc_.begin(); c0!=dc_.end(); ++c0) {
      if (*c0 < ncells) {
        AmanziMesh::Entity_ID_List fnums0;
        std::vector<int> dirs;
        mesh_->cell_get_faces_and_dirs(*c0, &fnums0, &dirs);

        *out_ << "  cond(" << *c0 << "): " << (*cond)("cell",*c0);
        for (unsigned int n=0; n!=fnums0.size(); ++n) *out_ << ",  " << (*cond)("face",fnums0[n]);
        *out_ << std::endl;
        *out_ << "  res(" << *c0 << ") (diff): " << (*res)("cell",*c0);
        for (unsigned int n=0; n!=fnums0.size(); ++n) *out_ << ",  " << (*res)("face",fnums0[n]);
        *out_ << std::endl;
      }
    }
  }
#endif

#if DEBUG_ICE_FLAG
  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_HIGH, true)) {
    *out_ << " WATER BUDGET" << std::endl;
    Teuchos::RCP<const CompositeVector> cond0 = S_next_->GetFieldData("overland_conductivity");
    Teuchos::RCP<const CompositeVector> dens0 = S_inter_->GetFieldData("surface_molar_density_liquid");
    Teuchos::RCP<const CompositeVector> cv0 = S_inter_->GetFieldData("surface_cell_volume");
    Teuchos::RCP<const CompositeVector> pd0 = S_inter_->GetFieldData("ponded_depth");
    Teuchos::RCP<const CompositeVector> uf0 = S_inter_->GetFieldData("unfrozen_fraction");

    unsigned int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
    for (std::vector<AmanziMesh::Entity_ID>::const_iterator c0=dc_.begin(); c0!=dc_.end(); ++c0) {
      if (*c0 < ncells) {

        double avail_water = (*dens0)("cell",*c0) * (*uf0)("cell",*c0) * (*pd0)("cell",*c0) * (*cv0)("cell",*c0);
        *out_ << "    available water (mol) = " << avail_water << std::endl;
        *out_ << "    outward flux (mol) = " << (*res)("cell",*c0) * h << std::endl;
        if ((*res)("cell",*c0) * h  > avail_water) {
          *out_ << " ADVECTING ICE!?!?" << std::endl;
          *out_ << "   conductivity (2014,2015) = " << (*cond0)("cell",*c0) << ", " << (*cond0)("cell",*c0) << std::endl;
        }
      }
    }
  }
#endif

  // accumulation term
  AddAccumulation_(res.ptr());
#if DEBUG_FLAG
  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_HIGH, true)) {
    unsigned int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
    for (std::vector<AmanziMesh::Entity_ID>::const_iterator c0=dc_.begin(); c0!=dc_.end(); ++c0) {
      if (*c0 < ncells) {
        *out_ << "  res(" << *c0 << ") (acc): " << (*res)("cell",*c0) << std::endl;
      }
    }
  }
#endif

  // add rhs load value
  AddSourceTerms_(res.ptr());
#if DEBUG_FLAG
  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_HIGH, true)) {
    unsigned int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
    for (std::vector<AmanziMesh::Entity_ID>::const_iterator c0=dc_.begin(); c0!=dc_.end(); ++c0) {
      if (*c0 < ncells) {
        *out_ << "  res(" << *c0 << ") (source): " << (*res)("cell",*c0) << std::endl;
      }
    }
  }
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
void OverlandHeadFlow::precon(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu) {
  // VerboseObject stuff.
  Teuchos::OSTab tab = getOSTab();

#if DEBUG_FLAG
  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_HIGH, true)) {
    *out_ << "Precon application:" << std::endl;
    unsigned int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
    for (std::vector<AmanziMesh::Entity_ID>::const_iterator c0=dc_.begin(); c0!=dc_.end(); ++c0) {
      if (*c0 < ncells) {
        AmanziMesh::Entity_ID_List fnums0;
        std::vector<int> dirs;
        mesh_->cell_get_faces_and_dirs(*c0, &fnums0, &dirs);

        *out_ << "  p(" << *c0 <<"): " << (*u->data())("cell",0,*c0);
        for (unsigned int n=0; n!=fnums0.size(); ++n) *out_ << ",  " << (*u->data())("face",fnums0[n]);
        *out_ << std::endl;
      }
    }
  }
#endif

  // apply the preconditioner
  preconditioner_->ApplyInverse(*u, Pu.ptr());

  // Dump correction
#if DEBUG_FLAG
  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_HIGH, true)) {
    *out_ << "Precon application:" << std::endl;
    unsigned int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
    for (std::vector<AmanziMesh::Entity_ID>::const_iterator c0=dc_.begin(); c0!=dc_.end(); ++c0) {
      if (*c0 < ncells) {
        AmanziMesh::Entity_ID_List fnums0;
        std::vector<int> dirs;
        mesh_->cell_get_faces_and_dirs(*c0, &fnums0, &dirs);

        *out_ << "  PC*p(" << *c0 <<") (pre-var): " << (*Pu->data())("cell",0,*c0);
        for (unsigned int n=0; n!=fnums0.size(); ++n) *out_ << ",  " << (*Pu->data())("face",fnums0[n]);
        *out_ << std::endl;
      }
    }
  }
#endif

  // tack on the variable change
  const Epetra_MultiVector& dh_dp =
      *S_next_->GetFieldData("dponded_depth_d"+key_)->ViewComponent("cell",false);
  Epetra_MultiVector& Pu_c = *Pu->data()->ViewComponent("cell",false);
  unsigned int ncells = Pu_c.MyLength();
  for (unsigned int c=0; c!=ncells; ++c) {
    Pu_c[0][c] /= dh_dp[0][c];
  }

  // Dump correction
#if DEBUG_FLAG
  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_HIGH, true)) {
    *out_ << "Precon application:" << std::endl;
    for (std::vector<AmanziMesh::Entity_ID>::const_iterator c0=dc_.begin(); c0!=dc_.end(); ++c0) {
      if (*c0 < ncells) {
        AmanziMesh::Entity_ID_List fnums0;
        std::vector<int> dirs;
        mesh_->cell_get_faces_and_dirs(*c0, &fnums0, &dirs);

        *out_ << "  PC*p(" << *c0 <<"): " << (*Pu->data())("cell",0,*c0);
        for (unsigned int n=0; n!=fnums0.size(); ++n) *out_ << ",  " << (*Pu->data())("face",fnums0[n]);
        *out_ << std::endl;
      }
    }
  }
#endif
};


// -----------------------------------------------------------------------------
// Update the preconditioner at time t and u = up
// -----------------------------------------------------------------------------
void OverlandHeadFlow::update_precon(double t, Teuchos::RCP<const TreeVector> up, double h) {
  // VerboseObject stuff.
  Teuchos::OSTab tab = getOSTab();
  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_EXTREME, true))
    *out_ << "Precon update at t = " << t << std::endl;

  Teuchos::RCP<const CompositeVector> eff_p = S_next_->GetFieldData("surface_effective_pressure");
  Teuchos::RCP<const CompositeVector> surf_p = S_next_->GetFieldData("surface_pressure");


  // update state with the solution up.
  ASSERT(std::abs(S_next_->time() - t) <= 1.e-4*t);
  PKDefaultBase::solution_to_state(up, S_next_);

  // update boundary conditions
  UpdateBoundaryConditionsMarkers_(S_next_.ptr());

  // update the rel perm according to the scheme of choice
  UpdatePermeabilityData_(S_next_.ptr());

  Teuchos::RCP<const CompositeVector> cond =
    S_next_->GetFieldData("upwind_overland_conductivity");

  // calculating the operator is done in 3 steps:
  // 1. Create all local matrices.
  mfd_preconditioner_->CreateMFDstiffnessMatrices(cond.ptr());
  mfd_preconditioner_->CreateMFDrhsVectors();

  // Patch up BCs in the case of zero conductivity
  FixBCsForPrecon_(S_next_.ptr());

  // 2.a: scale the cell by dh_dp
  S_next_->GetFieldEvaluator("ponded_depth")
      ->HasFieldDerivativeChanged(S_next_.ptr(), name_, key_);
  const Epetra_MultiVector& dh_dp =
      *S_next_->GetFieldData("dponded_depth_d"+key_)->ViewComponent("cell",false);
  const double& p_atm = *S_next_->GetScalarData("atmospheric_pressure");

  // 2.b Update local matrices diagonal with the accumulation terms.
  // -- update the accumulation derivatives
  S_next_->GetFieldEvaluator("surface_water_content")
      ->HasFieldDerivativeChanged(S_next_.ptr(), name_, key_);
  const Epetra_MultiVector& dwc_dp =
      *S_next_->GetFieldData("dsurface_water_content_d"+key_)
      ->ViewComponent("cell",false);

  // -- get other pieces
  const Epetra_MultiVector& head =
      *S_next_->GetFieldData(key_)->ViewComponent("cell",false);

  std::vector<double>& Acc_cells = mfd_preconditioner_->Acc_cells();
  unsigned int ncells = head.MyLength();

  for (unsigned int c=0; c!=ncells; ++c) {
    // - add accumulation terms, treating h as h_bar, which is allowed to go
    // - negative.  This results in a non-zero preconditioner when p < p_atm,
    // - resulting in a correction that should take p >= p_atm.
    Acc_cells[c] += dwc_dp[0][c] / dh_dp[0][c] / h;
  }

  // Assemble and precompute the Schur complement for inversion.
  // Note boundary conditions are in height variables.
  mfd_preconditioner_->ApplyBoundaryConditions(bc_markers_, bc_values_);

  if (coupled_to_subsurface_via_head_ ||
      coupled_to_subsurface_via_flux_ || assemble_preconditioner_) {
    if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_EXTREME, true))
      *out_ << "  assembling..." << std::endl;
    mfd_preconditioner_->AssembleGlobalMatrices();
  }

  if (tpfa_ && full_jacobian_) {
    if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_EXTREME, true))
      *out_ << "    including full Jacobian terms" << std::endl;

    // JACOBIAN?
    // These are already updated for UpdatePerm
    Teuchos::RCP<const CompositeVector> depth =
        S_next_->GetFieldData("ponded_depth");
    Teuchos::RCP<const CompositeVector> pres_elev =
        S_next_->GetFieldData("pres_elev");

    // conducitivity and dcond_dh
    S_next_->GetFieldEvaluator("overland_conductivity")
        ->HasFieldDerivativeChanged(S_next_.ptr(), name_, key_);
    Teuchos::RCP<const CompositeVector> cond =
        S_next_->GetFieldData("overland_conductivity");
    Teuchos::RCP<const CompositeVector> dcond_dp =
        S_next_->GetFieldData("doverland_conductivity_dsurface_pressure");
    CompositeVector dcond_dh(*dcond_dp);
    dcond_dh.ViewComponent("cell",false)->ReciprocalMultiply(1., dh_dp,
            *dcond_dp->ViewComponent("cell",false), 0.);


    // Krel_cell gets n_liq
    Teuchos::RCP<const CompositeVector> uw_cond =
        S_next_->GetFieldData("upwind_overland_conductivity");
    CompositeVector duw_cond_cell_dh(*uw_cond);
    duw_cond_cell_dh.ViewComponent("face",false)->PutScalar(0.);

    S_next_->GetFieldEvaluator("surface_molar_density_liquid")
        ->HasFieldDerivativeChanged(S_next_.ptr(), name_, key_);
    const Epetra_MultiVector& dn_liq_dp =
        *S_next_->GetFieldData("dsurface_molar_density_liquid_dsurface_pressure")
        ->ViewComponent("cell",false);
    duw_cond_cell_dh.ViewComponent("cell",false)
        ->ReciprocalMultiply(1., dh_dp, dn_liq_dp, 0.);

    // Add in the Jacobian
    tpfa_preconditioner_->AnalyticJacobian(*depth, *pres_elev, *cond, dcond_dh,
            *uw_cond, duw_cond_cell_dh);
  }

  // rescale to use as a pressure matrix if used in a coupler
  if (coupled_to_subsurface_via_head_ || coupled_to_subsurface_via_flux_) {
    ASSERT(tpfa_);
    Teuchos::RCP<Operators::MatrixMFD_TPFA> precon_tpfa =
        Teuchos::rcp_dynamic_cast<Operators::MatrixMFD_TPFA>(mfd_preconditioner_);
    ASSERT(precon_tpfa != Teuchos::null);
    Teuchos::RCP<Epetra_FECrsMatrix> Spp = precon_tpfa->TPFA();

    // Scale Spp by -dh/dp
    // if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_EXTREME, true))
    //   *out_ << "  scaling by dh/dp" << std::endl;

    // NOTE: dh/dp to take it to p variable, the negative sign is due to the
    //       equation being K/dz ( p - lambda ) = q = dwc/dt - div q_surf - Q,
    //     ---> K/dz * (p - lambda) - (dwc/dt - div q_surf - Q) = 0,
    //                              ^^ this sign is critical!
    // EpetraExt::RowMatrixToMatlabFile("TPFAbefore.txt", *Spp);
    Epetra_Vector dh_dp0(*dh_dp(0));
    for (unsigned int sc=0; sc!=ncells; ++sc) {
      dh_dp0[sc] = head[0][sc] > p_atm ? dh_dp[0][sc] : 0.;
      //      *out_ << " scaling by = " << dh_dp0[sc] << std::endl;
    }
    int ierr = Spp->RightScale(dh_dp0);
    ASSERT(!ierr);
    // EpetraExt::RowMatrixToMatlabFile("TPFAafter.txt", *Spp);
  }


  if (assemble_preconditioner_) {
    mfd_preconditioner_->ComputeSchurComplement(bc_markers_, bc_values_);
    mfd_preconditioner_->UpdatePreconditioner();
  }

  /*
  // dump the schur complement
  Teuchos::RCP<Epetra_FECrsMatrix> sc = mfd_preconditioner_->Schur();
  std::stringstream filename_s;
  filename_s << "schur_" << S_next_->cycle() << ".txt";
  EpetraExt::RowMatrixToMatlabFile(filename_s.str().c_str(), *sc);
  *out_ << "updated precon " << S_next_->cycle() << std::endl;
  */
};


void OverlandHeadFlow::set_preconditioner(const Teuchos::RCP<Operators::Matrix> precon) {
  preconditioner_ = precon;
  mfd_preconditioner_ = Teuchos::rcp_dynamic_cast<Operators::MatrixMFD>(precon);
  ASSERT(mfd_preconditioner_ != Teuchos::null);
  mfd_preconditioner_->set_symmetric(symmetric_);
  mfd_preconditioner_->SymbolicAssembleGlobalMatrices();
  mfd_preconditioner_->CreateMFDmassMatrices(Teuchos::null);
  mfd_preconditioner_->InitPreconditioner();
}


double OverlandHeadFlow::enorm(Teuchos::RCP<const TreeVector> u,
                       Teuchos::RCP<const TreeVector> du) {
  // Calculate water content at the solution.
  S_next_->GetFieldEvaluator("surface_water_content")
      ->HasFieldChanged(S_next_.ptr(), name_);
  const Epetra_MultiVector& wc = *S_next_->GetFieldData("surface_water_content")
      ->ViewComponent("cell",false);
  const Epetra_MultiVector& flux = *S_next_->GetFieldData("surface_flux")
      ->ViewComponent("face",false);
  double flux_max(0.);
  flux.NormInf(&flux_max);

  Teuchos::RCP<const CompositeVector> res = du->data();
  const Epetra_MultiVector& res_c = *res->ViewComponent("cell",false);
  const Epetra_MultiVector& res_f = *res->ViewComponent("face",false);
  const Epetra_MultiVector& height_f = *u->data()->ViewComponent("face",false);
  double h = S_next_->time() - S_inter_->time();

  // Cell error is based upon error in mass conservation relative to
  // the current water content
  double wc_base = 33.; // 1 cm on .25 m x .25 m blocks @ 55k n_liq ~ 33 mol
  double enorm_cell(0.);
  unsigned int ncells = res_c.MyLength();
  for (unsigned int c=0; c!=ncells; ++c) {
    double tmp = std::abs(h*res_c[0][c])
        / (atol_*wc_base + rtol_*std::abs(wc[0][c]));
    enorm_cell = std::max<double>(enorm_cell, tmp);
  }

  // Face error given by mismatch of flux, so relative to flux.
  double enorm_face(0.);
  unsigned int nfaces = res_f.MyLength();
  for (unsigned int f=0; f!=nfaces; ++f) {
    double tmp = 1.e-4 * std::abs(res_f[0][f]) / (atol_ + rtol_*flux_max);
    enorm_face = std::max<double>(enorm_face, tmp);
  }

  // Write out Inf norms too.
  Teuchos::OSTab tab = getOSTab();
  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_MEDIUM, true)) {
    double infnorm_c(0.), infnorm_f(0.);
    res_c.NormInf(&infnorm_c);
    res_f.NormInf(&infnorm_f);

#ifdef HAVE_MPI
    double buf_c(enorm_cell), buf_f(enorm_face);
    MPI_Allreduce(&buf_c, &enorm_cell, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&buf_f, &enorm_face, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#endif

    *out_ << "ENorm (cells) = " << enorm_cell << " (" << infnorm_c << ")" << std::endl;
    *out_ << "ENorm (faces) = " << enorm_face << " (" << infnorm_f << ")" << std::endl;
  }

  double enorm_val(std::max<double>(enorm_face, enorm_cell));
#ifdef HAVE_MPI
  double buf = enorm_val;
  MPI_Allreduce(&buf, &enorm_val, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#endif
  return enorm_val;
};

}  // namespace Flow
}  // namespace Amanzi
