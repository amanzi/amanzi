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

#define DEBUG_FLAG 0
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

#if DEBUG_FLAG
  AmanziMesh::Entity_ID_List faces, faces0;
  std::vector<int> dirs;
  mesh_->cell_get_faces_and_dirs(c0_, &faces0, &dirs);
  mesh_->cell_get_faces_and_dirs(c1_, &faces, &dirs);

  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_HIGH, true)) {
    S_next_->GetFieldEvaluator("pres_elev")->HasFieldChanged(S_next_.ptr(), name_);
    Teuchos::RCP<const CompositeVector> depth= S_next_->GetFieldData("ponded_depth");
    Teuchos::RCP<const CompositeVector> preselev= S_next_->GetFieldData("pres_elev");

    *out_ << "----------------------------------------------------------------" << std::endl;

    *out_ << "OverlandHead Residual calculation: T0 = " << t_old
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

  // zero out residual
  Teuchos::RCP<CompositeVector> res = g->data();
  res->PutScalar(0.0);

  // diffusion term, treated implicitly
  ApplyDiffusion_(S_next_.ptr(), res.ptr());

#if DEBUG_FLAG
  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_HIGH, true)) {
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
  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_HIGH, true)) {
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
  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_HIGH, true)) {
    *out_ << "  res0 (acc): " << (*res)("cell",c0_) << std::endl;
    *out_ << "  res1 (acc): " << (*res)("cell",c1_) << std::endl;

  }
#endif

  // add rhs load value
  AddSourceTerms_(res.ptr());
#if DEBUG_FLAG
  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_HIGH, true)) {
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

  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_HIGH, true)) {
    *out_ << "Precon application:" << std::endl;
    *out_ << "  p0: " << (*u->data())("cell",0,c0_) << " "
          << std::endl;
    //          << (*u->data())("face",0,0) << " " << (*u->data())("cell",0,cells[1]) << std::endl;
  }
#endif

  // apply the preconditioner
  preconditioner_->ApplyInverse(*u->data(), Pu->data().ptr());

  // Dump correction
#if DEBUG_FLAG
  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_HIGH, true)) {
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
  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_HIGH, true)) {
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

#if DEBUG_FLAG
  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_HIGH, true)) {
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
  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_HIGH, true)) {
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

  if (coupled_to_subsurface_via_flux_ || coupled_to_subsurface_via_head_) {
    // Coupled to subsurface, needs derivatives of coupling terms.
    const Epetra_MultiVector& dQ_dp =
        *S_next_->GetFieldData("doverland_source_from_subsurface_dsurface_pressure")
        ->ViewComponent("cell",false);

    for (int c=0; c!=ncells; ++c) {
      *out_ << "Acc diff terms, + " << Acc_cells[c] << std::endl;
      if (head[0][c] >= p_atm) {
        // - add accumulation terms, treating h as h, which goes zero.
        Acc_cells[c] += dwc_dp[0][c] / dh_dp[0][c] / h;
        *out_ << "Acc accum terms, +" << dwc_dp[0][c] / dh_dp[0][c] / h <<std::endl;
      } else {
        *out_ << "Acc accum terms, +" << 0. <<std::endl;
      }

      // add source term
      Acc_cells[c] -= dQ_dp[0][c] / dh_dp[0][c];
      *out_ << "Acc dQ terms, -" << dQ_dp[0][c] / dh_dp[0][c] << std::endl;
      *out_ << "dh_dp = " << dh_dp[0][c] << std::endl;
    }
  } else {
    // no Coupling term
    for (int c=0; c!=ncells; ++c) {
      // - add accumulation terms, treating h as h_bar, which is allowed to go
      // - negative.  This results in a non-zero preconditioner when p < p_atm,
      // - resulting in a correction that should take p >= p_atm.
      Acc_cells[c] += dwc_dp[0][c] / dh_dp[0][c] / h;
      //      *out_ << " adding acc term = " << dwc_dp[0][c] / dh_dp[0][c] / h << std::endl;
    }
  }

  // Assemble and precompute the Schur complement for inversion.
  // Note boundary conditions are in height variables.
  mfd_preconditioner_->ApplyBoundaryConditions(bc_markers_, bc_values_);

  if (coupled_to_subsurface_via_full_) {
    mfd_preconditioner_->AssembleGlobalMatrices();

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

}  // namespace Flow
}  // namespace Amanzi
