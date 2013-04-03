/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -----------------------------------------------------------------------------
This is the overland flow component of ATS.
License: BSD
Authors: Ethan Coon (ecoon@lanl.gov)
----------------------------------------------------------------------------- */

#include "Epetra_FECrsMatrix.h"
#include "EpetraExt_RowMatrixOut.h"
#include "boost/math/special_functions/fpclassify.hpp"
#include "matrix_mfd_tpfa.hh"

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

  // bookkeeping
  double h = t_new - t_old;
  ASSERT(std::abs(S_inter_->time() - t_old) < 1.e-4*h);
  ASSERT(std::abs(S_next_->time() - t_new) < 1.e-4*h);

  Teuchos::RCP<CompositeVector> u = u_new->data();

  // zero out residual
  Teuchos::RCP<CompositeVector> res = g->data();
  res->PutScalar(0.0);

#if DEBUG_FLAG
  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_HIGH, true) &&
      c0_ < res->size("cell",false) && c1_ < res->size("cell",false)) {
    AmanziMesh::Entity_ID_List faces, faces0;
    std::vector<int> dirs;
    mesh_->cell_get_faces_and_dirs(c0_, &faces0, &dirs);
    mesh_->cell_get_faces_and_dirs(c1_, &faces, &dirs);

    S_next_->GetFieldEvaluator("pres_elev")->HasFieldChanged(S_next_.ptr(), name_);
    Teuchos::RCP<const CompositeVector> depth= S_next_->GetFieldData("ponded_depth");
    Teuchos::RCP<const CompositeVector> preselev= S_next_->GetFieldData("pres_elev");


    *out_ << "----------------------------------------------------------------" << std::endl;

    *out_ << "Residual calculation: T0 = " << t_old
          << " T1 = " << t_new << " H = " << h << std::endl;
    *out_ << std::setprecision(15);
    *out_ << "  p0: " << (*u)("cell",0,c0_) << std::endl;
    *out_ << "  h0: " << (*depth)("cell",0,c0_) << " "
          << (*depth)("face",0,faces0[0]) << " " << (*depth)("face",0,faces0[1]) << " "
          << (*depth)("face",0,faces0[2]) << " " << (*depth)("face",0,faces0[3]) << std::endl;
    *out_ << "  hz0: " << (*preselev)("cell",0,c0_) << " "
          << (*preselev)("face",0,faces0[0]) << " " << (*preselev)("face",0,faces0[1]) << " "
          << (*preselev)("face",0,faces0[2]) << " " << (*preselev)("face",0,faces0[3]) << std::endl;
    *out_ << " --" << std::endl;
    *out_ << "  p1: " << (*u)("cell",c1_) << std::endl;
    *out_ << "  h1: " << (*depth)("cell",c1_) << " "
          << (*depth)("face",faces[0]) << " " << (*depth)("face",faces[1]) << " "
          << (*depth)("face",faces[2]) << " " << (*depth)("face",faces[3]) << std::endl;
    *out_ << "  hz1: " << (*preselev)("cell",c1_) << " "
          << (*preselev)("face",faces[0]) << " " << (*preselev)("face",faces[1]) << " "
          << (*preselev)("face",faces[2]) << " " << (*preselev)("face",faces[3]) << std::endl;
  }
#endif

  // pointer-copy temperature into state and update any auxilary data
  solution_to_state(u_new, S_next_);

  // update boundary conditions
  bc_pressure_->Compute(t_new);
  bc_head_->Compute(t_new);
  bc_flux_->Compute(t_new);
  UpdateBoundaryConditions_(S_next_.ptr());

  // update the rel perm according to the scheme of choice.
  UpdatePermeabilityData_(S_next_.ptr());

  // update the stiffness matrix
  Teuchos::RCP<const CompositeVector> cond =
    S_next_->GetFieldData("upwind_overland_conductivity", name_);
  matrix_->CreateMFDstiffnessMatrices(cond.ptr());

  // Patch up BCs in the case of zero conductivity
  FixBCsForOperator_(S_next_.ptr());

  // diffusion term, treated implicitly
  ApplyDiffusion_(S_next_.ptr(), res.ptr());

#if DEBUG_FLAG
  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_HIGH, true) &&
      c0_ < res->size("cell",false) && c1_ < res->size("cell",false)) {
    AmanziMesh::Entity_ID_List faces, faces0;
    std::vector<int> dirs;
    mesh_->cell_get_faces_and_dirs(c0_, &faces0, &dirs);
    mesh_->cell_get_faces_and_dirs(c1_, &faces, &dirs);

    *out_ << "  cond0 (diff): " << (*cond)("cell",c0_) << " "
          << (*cond)("face",faces0[0]) << " " << (*cond)("face",faces0[1]) << " "
          << (*cond)("face",faces0[2]) << " " << (*cond)("face",faces0[3]) << std::endl;
    *out_ << "  cond1 (diff): " << (*cond)("cell",c1_) << " "
          << (*cond)("face",faces[0]) << " " << (*cond)("face",faces[1]) << " "
          << (*cond)("face",faces[2]) << " " << (*cond)("face",faces[3]) << std::endl;
    *out_ << "  res0 (diff): " << (*res)("cell",c0_) << " "
             << (*res)("face",faces0[0]) << " " << (*res)("face",faces0[1]) << " "
             << (*res)("face",faces0[2]) << " " << (*res)("face",faces0[3]) << std::endl;
    *out_ << "  res1 (diff): " << (*res)("cell",c1_) << " "
             << (*res)("face",faces[0]) << " " << (*res)("face",faces[1]) << " "
             << (*res)("face",faces[2]) << " " << (*res)("face",faces[3]) << std::endl;
  }

#endif

#if DEBUG_ICE_FLAG
  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_HIGH, true) &&
      c0_ < res->size("cell",false) && c1_ < res->size("cell",false)) {
    *out_ << " WATER BUDGET" << std::endl;
    Teuchos::RCP<const CompositeVector> cond0 = S_next_->GetFieldData("overland_conductivity");
    Teuchos::RCP<const CompositeVector> dens0 = S_inter_->GetFieldData("surface_molar_density_liquid");
    Teuchos::RCP<const CompositeVector> cv0 = S_inter_->GetFieldData("surface_cell_volume");
    Teuchos::RCP<const CompositeVector> pd0 = S_inter_->GetFieldData("ponded_depth");
    Teuchos::RCP<const CompositeVector> uf0 = S_inter_->GetFieldData("unfrozen_fraction");
    double avail_water = (*dens0)("cell",c1_) * (*uf0)("cell",c1_) * (*pd0)("cell",c1_) * (*cv0)("cell",c1_);
    *out_ << "    available water (mol) = " << avail_water << std::endl;
    *out_ << "    outward flux (mol) = " << (*res)("cell",c1_) * h << std::endl;
    if ((*res)("cell",c1_) * h  > avail_water) {
      *out_ << " ADVECTING ICE!?!?" << std::endl;
      *out_ << "   conductivity (2014,2015) = " << (*cond0)("cell",c0_) << ", " << (*cond0)("cell",c0_) << std::endl;
    }
  }
#endif

  // accumulation term
  AddAccumulation_(res.ptr());
#if DEBUG_FLAG
  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_HIGH, true) &&
      c0_ < res->size("cell",false) && c1_ < res->size("cell",false)) {
    *out_ << "  res0 (acc): " << (*res)("cell",c0_) << std::endl;
    *out_ << "  res1 (acc): " << (*res)("cell",c1_) << std::endl;

  }
#endif

  // add rhs load value
  AddSourceTerms_(res.ptr());
#if DEBUG_FLAG
  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_HIGH, true) &&
      c0_ < res->size("cell",false) && c1_ < res->size("cell",false)) {
    *out_ << "  res0 (source): " << (*res)("cell",c0_) << std::endl;
    *out_ << "  res1 (source): " << (*res)("cell",c1_) << std::endl;
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
  //  AmanziMesh::Entity_ID_List cells;
  //  mesh_->face_get_cells(0, AmanziMesh::USED, &cells);

  AmanziMesh::Entity_ID_List faces;
  std::vector<int> dirs;
  mesh_->cell_get_faces_and_dirs(c0_, &faces, &dirs);

  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_HIGH, true) &&
      c0_ < u->data()->size("cell",false) && c1_ < u->data()->size("cell",false)) {
    *out_ << "Precon application:" << std::endl;
    *out_ << "  p0: " << (*u->data())("cell",0,c0_) << " "
          << std::endl;
    //          << (*u->data())("face",0,0) << " " << (*u->data())("cell",0,cells[1]) << std::endl;
  }
#endif

  // apply the preconditioner
  preconditioner_->ApplyInverse(*u, Pu.ptr());

  // Dump correction
#if DEBUG_FLAG
  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_HIGH, true) &&
      c0_ < u->data()->size("cell",false) && c1_ < u->data()->size("cell",false)) {
    *out_ << "  PC*p0, pre-var change: " << (*Pu->data())("cell",0,c0_) << " "
          << std::endl;
        //          << (*Pu->data())("face",0,0) << " " << (*Pu->data())("cell",0,cells[1]) << std::endl;
  }
#endif

  // tack on the variable change
  const Epetra_MultiVector& dh_dp =
      *S_next_->GetFieldData("dponded_depth_d"+key_)->ViewComponent("cell",false);
  Epetra_MultiVector& Pu_c = *Pu->data()->ViewComponent("cell",false);
  int ncells = Pu_c.MyLength();
  for (int c=0; c!=ncells; ++c) {
    Pu_c[0][c] /= dh_dp[0][c];
  }

  // Dump correction
#if DEBUG_FLAG
  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_HIGH, true) &&
      c0_ < u->data()->size("cell",false) && c1_ < u->data()->size("cell",false)) {
    *out_ << "  PC*p0: " << (*Pu->data())("cell",0,c0_) << " "
          << std::endl;
        //          << (*Pu->data())("face",0,0) << " " << (*Pu->data())("cell",0,cells[1]) << std::endl;
  }
#endif
};


// -----------------------------------------------------------------------------
// Update the preconditioner at time t and u = up
// -----------------------------------------------------------------------------
void OverlandHeadFlow::update_precon(double t, Teuchos::RCP<const TreeVector> up, double h) {
  // VerboseObject stuff.
  Teuchos::OSTab tab = getOSTab();

  Teuchos::RCP<const CompositeVector> eff_p = S_next_->GetFieldData("surface_effective_pressure");
  Teuchos::RCP<const CompositeVector> surf_p = S_next_->GetFieldData("surface_pressure");

#if DEBUG_FLAG
  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_EXTREME, true) &&
      c0_ < eff_p->size("cell",false) && c1_ < eff_p->size("cell",false)) {
    *out_ << "Precon update at t = " << t << std::endl;
    *out_ << "  p0: " << (*up->data())("cell",0,c0_) << " "
          << std::endl;
    //              << (*up->data())("face",0,0) << std::endl;
  }
#endif

  // update state with the solution up.
  ASSERT(std::abs(S_next_->time() - t) <= 1.e-4*t);
  PKDefaultBase::solution_to_state(up, S_next_);

  // update boundary conditions
  UpdateBoundaryConditionsMarkers_(S_next_.ptr());

  // update the rel perm according to the scheme of choice
  UpdatePermeabilityData_(S_next_.ptr());

  Teuchos::RCP<const CompositeVector> cond =
    S_next_->GetFieldData("upwind_overland_conductivity");
  mfd_preconditioner_->CreateMFDstiffnessMatrices(cond.ptr());

  // Patch up BCs in the case of zero conductivity
  FixBCsForPrecon_(S_next_.ptr());

#if DEBUG_FLAG
  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_EXTREME, true) &&
      c0_ < cond->size("cell",false) && c1_ < cond->size("cell",false)) {
    *out_ << "  conductivity0: " << (*cond)("cell",0,c0_) << " "
          << (*cond)("face",0,0) << std::endl;
    *out_ << "  bcs: " << std::endl;
    *out_ << "    0: " << bc_markers_[0] << ", " << bc_values_[0] <<std::endl;
  }
#endif

  // calculating the operator is done in 3 steps:
  // 1. Create all local matrices.
  mfd_preconditioner_->CreateMFDrhsVectors();

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
  int ncells = head.MyLength();

  for (int c=0; c!=ncells; ++c) {
    // - add accumulation terms, treating h as h_bar, which is allowed to go
    // - negative.  This results in a non-zero preconditioner when p < p_atm,
    // - resulting in a correction that should take p >= p_atm.
    Acc_cells[c] += dwc_dp[0][c] / dh_dp[0][c] / h;
  }

  // Assemble and precompute the Schur complement for inversion.
  // Note boundary conditions are in height variables.
  mfd_preconditioner_->ApplyBoundaryConditions(bc_markers_, bc_values_);

  if (coupled_to_subsurface_via_head_ || coupled_to_subsurface_via_flux_) {
    mfd_preconditioner_->AssembleGlobalMatrices();

    if (full_jacobian_) {
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

    // TPFA
    Teuchos::RCP<Operators::MatrixMFD_TPFA> precon_tpfa =
        Teuchos::rcp_dynamic_cast<Operators::MatrixMFD_TPFA>(mfd_preconditioner_);
    ASSERT(precon_tpfa != Teuchos::null);
    Teuchos::RCP<Epetra_FECrsMatrix> Spp = precon_tpfa->TPFA();

    // Scale Spp by -dh/dp

    // NOTE: dh/dp to take it to p variable, the negative sign is due to the
    //       equation being K/dz ( p - lambda ) = q = dwc/dt - div q_surf - Q,
    //     ---> K/dz * (p - lambda) - (dwc/dt - div q_surf - Q) = 0,
    //                              ^^ this sign is critical!
    //    EpetraExt::RowMatrixToMatlabFile("TPFAbefore.txt", *Spp);
    Epetra_Vector dh_dp0(*dh_dp(0));
    for (int c=0; c!=ncells; ++c) {
      dh_dp0[c] = head[0][c] >= p_atm ? dh_dp[0][c] : 0.;
      //      *out_ << " scaling by = " << dh_dp0[c] << std::endl;
    }
    int ierr = Spp->RightScale(dh_dp0);
    ASSERT(!ierr);
    //    EpetraExt::RowMatrixToMatlabFile("TPFAafter.txt", *Spp);

  } else if (assemble_preconditioner_) {
    mfd_preconditioner_->AssembleGlobalMatrices();
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


  // print the rel perm
  Teuchos::RCP<const CompositeVector> num_rel_perm =
      S_next_->GetFieldData("upwind_overland_conductivity");
  Teuchos::RCP<const CompositeVector> rel_perm =
      S_next_->GetFieldData("overland_conductivity");
  *out_ << "REL PERM: " << std::endl;
  rel_perm->Print(*out_);
  *out_ << std::endl;
  *out_ << "UPWINDED REL PERM: " << std::endl;
  num_rel_perm->Print(*out_);
  */
};


void OverlandHeadFlow::set_preconditioner(const Teuchos::RCP<Operators::Matrix> precon) {
  preconditioner_ = precon;
  mfd_preconditioner_ = Teuchos::rcp_dynamic_cast<Operators::MatrixMFD>(precon);
  ASSERT(mfd_preconditioner_ != Teuchos::null);
  mfd_preconditioner_->SetSymmetryProperty(symmetric_);
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

  Teuchos::RCP<const CompositeVector> res = du->data();
  const Epetra_MultiVector& res_c = *res->ViewComponent("cell",false);
  const Epetra_MultiVector& res_f = *res->ViewComponent("face",false);
  const Epetra_MultiVector& height_f = *u->data()->ViewComponent("face",false);
  double h = S_next_->time() - S_inter_->time();

  // Cell error is based upon error in mass conservation relative to
  // the current water content
  double wc_base = 33.; // 1 cm on .25 m x .25 m blocks @ 55k n_liq ~ 33 mol
  double enorm_cell(0.);
  int ncells = res_c.MyLength();
  for (int c=0; c!=ncells; ++c) {
    double tmp = std::abs(h*res_c[0][c])
        / (atol_*wc_base + rtol_*std::abs(wc[0][c]));
    enorm_cell = std::max<double>(enorm_cell, tmp);
  }


#if DEBUG_FLAG
  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_EXTREME, true) &&
      c0_ < res->size("cell",false) && c1_ < res->size("cell",false)) {
    *out_ << "ERRORS:" << std::endl;
    *out_ << "  cell0: h*res = " << std::abs(h*res_c[0][c0_])
          << ", abs= " << atol_*wc_base
          << ", rel= " << rtol_*std::abs(wc[0][c0_]) << std::endl;
    *out_ << "  cell1: h*res = " << std::abs(h*res_c[0][c1_])
          << ", abs= " << atol_*wc_base
          << ", rel= " << rtol_*std::abs(wc[0][c1_]) << std::endl;
  }
#endif

  // Face error give by heights?  Loose tolerance
  double enorm_face(0.);
  int nfaces = res_f.MyLength();
  for (int f=0; f!=nfaces; ++f) {
    double tmp = std::abs(res_f[0][f]) / (atol_ + rtol_*(height_f[0][f]));
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

    *out_ << "ENorm (cells) = " << enorm_cell << " (" << infnorm_c << ")  " << std::endl;
    *out_ << "ENorm (faces) = " << enorm_face << " (" << infnorm_f << ")  " << std::endl;
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
