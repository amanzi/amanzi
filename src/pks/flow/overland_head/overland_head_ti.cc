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

namespace Amanzi {
namespace Flow {

#define DEBUG_FLAG 1
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
  AmanziMesh::Entity_ID_List faces;
  std::vector<int> dirs;
  mesh_->cell_get_faces_and_dirs(0, &faces, &dirs);

  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_HIGH, true)) {
    S_next_->GetFieldEvaluator("pres_elev")->HasFieldChanged(S_next_.ptr(), name_);
    Teuchos::RCP<const CompositeVector> depth= S_next_->GetFieldData("ponded_depth");
    Teuchos::RCP<const CompositeVector> preselev= S_next_->GetFieldData("pres_elev");

    *out_ << "OverlandHeadFlow Residual calculation:" << std::endl;
    *out_ << std::setprecision(15);
    *out_ << "  p0: " << (*u)("cell",0,0) << " "
          << std::endl;
        //          << (*u)("face",0,faces[0]) << " " << (*u)("face",0,faces[1]) << " "
        //          << (*u)("face",0,faces[2]) << " " << (*u)("face",0,faces[3]) << std::endl;
    *out_ << "  h0: " << (*depth)("cell",0,0) << " "
          << std::endl;
    //          << (*depth)("face",0,faces[0]) << " " << (*depth)("face",0,faces[1]) << " "
    //          << (*depth)("face",0,faces[2]) << " " << (*depth)("face",0,faces[3]) << std::endl;
    *out_ << "  hz0: " << (*preselev)("cell",0,0) << " "
          << std::endl;
    //          << (*preselev)("face",0,faces[0]) << " " << (*preselev)("face",0,faces[1]) << " "
    //          << (*preselev)("face",0,faces[2]) << " " << (*preselev)("face",0,faces[3]) << std::endl;
  }
#endif

  // pointer-copy temperature into state and update any auxilary data
  solution_to_state(u_new, S_next_);

  // update boundary conditions
  bc_pressure_->Compute(t_new);
  bc_flux_->Compute(t_new);
  UpdateBoundaryConditions_(S_next_.ptr());

  // zero out residual
  Teuchos::RCP<CompositeVector> res = g->data();
  res->PutScalar(0.0);

  // diffusion term, treated implicitly
  ApplyDiffusion_(S_next_.ptr(), res.ptr());

#if DEBUG_FLAG
  Teuchos::RCP<const CompositeVector> cond =
      S_next_->GetFieldData("upwind_overland_conductivity", name_);

  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_HIGH, true)) {
    *out_ << "  conductivty0 (diff): " << (*cond)("cell",0,0) << " "
          << (*cond)("face",0,0) << " " << (*cond)("face",0,1) << " "
          << (*cond)("face",0,2) << " " << (*cond)("face",0,3) << std::endl;
    *out_ << "  res0 (diff): " << (*res)("cell",0,0) << " "
          << std::endl;
    //          << (*res)("face",0,faces[0]) << " " << (*res)("face",0,faces[1]) << " "
    //          << (*res)("face",0,faces[2]) << " " << (*res)("face",0,faces[3]) << std::endl;
  }
#endif

  // accumulation term
  AddAccumulation_(res.ptr());
#if DEBUG_FLAG
  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_HIGH, true)) {
    *out_ << "  res0 (acc): " << (*res)("cell",0,0) << " "
          << std::endl;
    //          << (*res)("face",0,faces[0]) << " " << (*res)("face",0,faces[1]) << " "
    //          << (*res)("face",0,faces[2]) << " " << (*res)("face",0,faces[3]) << std::endl;
  }
#endif

  // add rhs load value
  AddSourceTerms_(res.ptr());
#if DEBUG_FLAG
  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_HIGH, true)) {
    *out_ << "  res0 (source): " << (*res)("cell",0,0) << " "
          << std::endl;
    //          << (*res)("face",0,faces[0]) << " " << (*res)("face",0,faces[1]) << " "
    //          << (*res)("face",0,faces[2]) << " " << (*res)("face",0,faces[3]) << std::endl;
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

  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_HIGH, true)) {
    *out_ << "Precon application:" << std::endl;
    *out_ << "  p0: " << (*u->data())("cell",0,0) << " "
          << std::endl;
    //          << (*u->data())("face",0,0) << " " << (*u->data())("cell",0,cells[1]) << std::endl;
  }
#endif

  // apply the preconditioner
  preconditioner_->ApplyInverse(*u->data(), Pu->data().ptr());

  // Dump correction
#if DEBUG_FLAG
  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_HIGH, true)) {
    *out_ << "  PC*p0, pre-var change: " << (*Pu->data())("cell",0,0) << " "
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
    *out_ << "  PC*p0: " << (*Pu->data())("cell",0,0) << " "
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
    *out_ << "  p0: " << (*up->data())("cell",0,0) << " "
          << std::endl;
    //              << (*up->data())("face",0,0) << std::endl;
  }
#endif

  // update state with the solution up.
  ASSERT(std::abs(S_next_->time() - t) <= 1.e-4*t);
  PKDefaultBase::solution_to_state(up, S_next_);

  // update the rel perm according to the scheme of choice
  UpdatePermeabilityData_(S_next_.ptr());

  // update boundary conditions
  bc_pressure_->Compute(S_next_->time());
  bc_flux_->Compute(S_next_->time());
  UpdateBoundaryConditionsNoElev_(S_next_.ptr());
  //  UpdateBoundaryConditions_(S_next_.ptr());

  Teuchos::RCP<const CompositeVector> cond =
    S_next_->GetFieldData("upwind_overland_conductivity");

#if DEBUG_FLAG
  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_HIGH, true)) {
    *out_ << "  conductivity0: " << (*cond)("cell",0,0) << " "
          << (*cond)("face",0,0) << std::endl;
    *out_ << "  bcs: " << std::endl;
    *out_ << "    0: " << bc_markers_[0] << ", " << bc_values_[0] <<std::endl;
  }
#endif

  // calculating the operator is done in 3 steps:
  // 1. Create all local matrices.
  preconditioner_->CreateMFDstiffnessMatrices(cond.ptr());
  preconditioner_->CreateMFDrhsVectors();

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
  const Epetra_MultiVector& cv =
      *S_next_->GetFieldData("surface_cell_volume")->ViewComponent("cell",false);

  std::vector<double>& Acc_cells = preconditioner_->Acc_cells();
  int ncells = cv.MyLength();

  if (!is_coupling_term_) {
    for (int c=0; c!=ncells; ++c) {
      // - add accumulation terms, treating h as h_bar, which is allowed to go
      // - negative.  This results in a non-zero preconditioner when p < p_atm,
      // - resulting in a correction that should take p >= p_atm.
      *out_ << "Acc diff terms, " << Acc_cells[c] <<std::endl;
      Acc_cells[c] += dwc_dp[0][c] / dh_dp[0][c] * cv[0][c] / h;
      *out_ << "Acc accum terms, +" << dwc_dp[0][c] / dh_dp[0][c] * cv[0][c] / h <<std::endl;
      *out_ << "dh_dp = " << dh_dp[0][c] << std::endl;
    }
  } else {
    // - add source terms, which come from the flux from surface to
    // - subsurface.  This vector was filled in by the Richards PK.
    const Epetra_MultiVector& dQ_dp =
        *S_next_->GetFieldData("doverland_source_from_subsurface_dsurface_pressure")
        ->ViewComponent("cell",false);

    for (int c=0; c!=ncells; ++c) {
      *out_ << "Acc diff terms, + " << Acc_cells[c] << std::endl;
      if (head[0][c] >= p_atm) {
        // - add accumulation terms, treating h as h, which goes zero.
        Acc_cells[c] += dwc_dp[0][c] / dh_dp[0][c] * cv[0][c] / h;
        *out_ << "Acc accum terms, +" << dwc_dp[0][c] / dh_dp[0][c] * cv[0][c] / h <<std::endl;
      }

      // add source term
      Acc_cells[c] -= dQ_dp[0][c] / dh_dp[0][c];
      *out_ << "Acc dQ terms, -" << dQ_dp[0][c] / dh_dp[0][c] << std::endl;
      *out_ << "dh_dp = " << dh_dp[0][c] << std::endl;
    }
  }

  // Assemble and precompute the Schur complement for inversion.
  // Note boundary conditions are in height variables.
  preconditioner_->ApplyBoundaryConditions(bc_markers_, bc_values_);

  if (assemble_preconditioner_) {
    preconditioner_->AssembleGlobalMatrices();
    //    preconditioner_->ComputeSchurComplement(bc_markers_, bc_values_);
    //    preconditioner_->UpdatePreconditioner();
  }

  /*
  // dump the schur complement
  Teuchos::RCP<Epetra_FECrsMatrix> sc = preconditioner_->Schur();
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



}  // namespace Flow
}  // namespace Amanzi
