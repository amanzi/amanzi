/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Interface for the derived MPC for water coupling between surface and subsurface.

In this method, a Dirichlet BC is used on the subsurface boundary for
the operator, but the preconditioner is for the residual system with no
extra unknowns.  On the surface, the TPFA is used, resulting in a
subsurface-face-only Schur complement that captures all terms.

------------------------------------------------------------------------- */
#include "EpetraExt_RowMatrixOut.h"

#include "matrix_mfd_surf.hh"
#include "matrix_mfd_tpfa.hh"
#include "pk_physical_bdf_base.hh"

#include "mpc_surface_subsurface_residual_coupler.hh"


namespace Amanzi {

#define DEBUG_FLAG 1

RegisteredPKFactory<MPCSurfaceSubsurfaceResidualCoupler>
MPCSurfaceSubsurfaceResidualCoupler::reg_("surface-subsurface residual coupler");

// -- Setup data.
void MPCSurfaceSubsurfaceResidualCoupler::setup(const Teuchos::Ptr<State>& S) {
  MPCSurfaceSubsurfaceCoupler::setup(S);

  // Get the domain's preconditioner and replace it with a MatrixMFD_Surf.
  Teuchos::RCP<Operators::Matrix> precon = domain_pk_->preconditioner();
  Teuchos::RCP<Operators::MatrixMFD> mfd_precon =
    Teuchos::rcp_dynamic_cast<Operators::MatrixMFD>(precon);
  ASSERT(mfd_precon != Teuchos::null);

  preconditioner_ = Teuchos::rcp(new Operators::MatrixMFD_Surf(*mfd_precon, surf_mesh_));

  // Get the surface's preconditioner and ensure it is TPFA.
  Teuchos::RCP<Operators::Matrix> surf_precon = surf_pk_->preconditioner();
  surf_preconditioner_ =
    Teuchos::rcp_dynamic_cast<Operators::MatrixMFD_TPFA>(surf_precon);
  ASSERT(surf_preconditioner_ != Teuchos::null);

  // set the surface A in the MFD_Surf.
  preconditioner_->set_surface_A(surf_preconditioner_);

  // give the PCs back to the PKs
  domain_pk_->set_preconditioner(preconditioner_);
  surf_pk_->set_preconditioner(surf_preconditioner_);

}

// -------------------------------------------------------------
// Special residual function evaluation process
// -------------------------------------------------------------
void MPCSurfaceSubsurfaceResidualCoupler::fun(double t_old, double t_new,
        Teuchos::RCP<TreeVector> u_old, Teuchos::RCP<TreeVector> u_new,
        Teuchos::RCP<TreeVector> g) {
  // ensure the source term from subsurface -> surface is 0
  Teuchos::RCP<CompositeVector> source =
      S_next_->GetFieldData("overland_source_from_subsurface", domain_pk_name_);
  source->PutScalar(0.);

  // evaluate the residual for overland flow
  Teuchos::RCP<TreeVector> overland_u_new = u_new->SubVector(surf_pk_name_);
  Teuchos::RCP<TreeVector> overland_g = g->SubVector(surf_pk_name_);
  surf_pk_->fun(t_old, t_new, Teuchos::null, overland_u_new, overland_g);

  // Stuff the residual value into the source term.  This is the mass
  // imbalance on the surface, and is used in the subsurface PK.
  *source->ViewComponent("cell",false) = *overland_g->data()->ViewComponent("cell",false);
  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_EXTREME, true)) {
    Teuchos::OSTab tab = getOSTab();
    source->Print(*out_);
  }

  // Evaluate the Richards residual, which will update the source with the true value.
  Teuchos::RCP<TreeVector> richards_u_new = u_new->SubVector(domain_pk_name_);
  Teuchos::RCP<TreeVector> richards_g = g->SubVector(domain_pk_name_);
  domain_pk_->fun(t_old, t_new, Teuchos::null, richards_u_new, richards_g);

  // Set the overland cell residual to 0.
  overland_g->data()->ViewComponent("cell",false)->PutScalar(0.);
}

// applies preconditioner to u and returns the result in Pu
void MPCSurfaceSubsurfaceResidualCoupler::precon(Teuchos::RCP<const TreeVector> u,
        Teuchos::RCP<TreeVector> Pu) {

  Teuchos::RCP<const CompositeVector> surf_u = u->SubVector(surf_pk_name_)->data();
  Teuchos::RCP<const CompositeVector> domain_u = u->SubVector(domain_pk_name_)->data();
  Teuchos::RCP<CompositeVector> domain_Pu = Pu->SubVector(domain_pk_name_)->data();
  Teuchos::RCP<CompositeVector> surf_Pu = Pu->SubVector(surf_pk_name_)->data();

  // Apply the combined preconditioner to the subsurface residual
  preconditioner_->ApplyInverse(*domain_u, domain_Pu.ptr());

#if DEBUG_FLAG
  Teuchos::OSTab tab = getOSTab();
  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_HIGH, true)) {
    AmanziMesh::Entity_ID_List fnums,fnums0;
    std::vector<int> dirs;
    domain_mesh_->cell_get_faces_and_dirs(c0_, &fnums0, &dirs);
    domain_mesh_->cell_get_faces_and_dirs(c1_, &fnums, &dirs);

    *out_ << "Preconditioner application" << std::endl;
    *out_ << " SubSurface precon:" << std::endl;
    *out_ << "  p0: " << (*domain_u)("cell",c0_) << " "
        //          << (*domain_u)("face",fnums0[0])
          << (*domain_u)("face",fnums0[0]) << ", "
          << (*domain_u)("face",fnums0[1]) << ", "
          << (*domain_u)("face",fnums0[2]) << ", "
          << (*domain_u)("face",fnums0[3])
          << std::endl;
    *out_ << "  p1: " << (*domain_u)("cell",c1_) << " "
        //          << (*domain_u)("face",fnums[0])
          << (*domain_u)("face",fnums[0]) << ", "
          << (*domain_u)("face",fnums[1]) << ", "
          << (*domain_u)("face",fnums[2]) << ", "
          << (*domain_u)("face",fnums[3])
          << std::endl;
    *out_ << "  PC*p0: " << (*domain_Pu)("cell",c0_) << " "
        //          << (*domain_Pu)("face",fnums0[0])
          << (*domain_Pu)("face",fnums0[0]) << ", "
          << (*domain_Pu)("face",fnums0[1]) << ", "
          << (*domain_Pu)("face",fnums0[2]) << ", "
          << (*domain_Pu)("face",fnums0[3])
          << std::endl;
    *out_ << "  PC*p1: " << (*domain_Pu)("cell",c1_) << " "
        //          << (*domain_Pu)("face",fnums[0])
          << (*domain_Pu)("face",fnums[0]) << ", "
          << (*domain_Pu)("face",fnums[1]) << ", "
          << (*domain_Pu)("face",fnums[2]) << ", "
          << (*domain_Pu)("face",fnums[3])
          << std::endl;
  }
#endif

  int ncells_surf = surf_mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);

  // Damp surface corrections
  const double& patm = *S_next_->GetScalarData("atmospheric_pressure");
  Epetra_MultiVector& domain_Pu_f = *domain_Pu->ViewComponent("face",false);

  const Epetra_MultiVector& domain_p_f = *S_next_->GetFieldData("pressure")
      ->ViewComponent("face",false);

  if (damping_coef_ > 0.) {
    for (int cs=0; cs!=ncells_surf; ++cs) {
      AmanziMesh::Entity_ID f = surf_mesh_->entity_get_parent(AmanziMesh::CELL, cs);

      double p_old = domain_p_f[0][f];
      double p_new = p_old - domain_Pu_f[0][f];
      if ((p_new > patm) && (std::abs(domain_Pu_f[0][f]) > damping_cutoff_)) {
        domain_Pu_f[0][f] *= damping_coef_;
        std::cout << "  DAMPING: p_old = " << p_old << ", p_new = " << p_new << ", p_damped = " << p_old - domain_Pu_f[0][f] << std::endl;
      }
    }
  }

  // Correction applies to both the domain face and the surface cell.
  // Correct for drift
  Epetra_MultiVector& surf_Pu_c = *surf_Pu->ViewComponent("cell",false);
  const Epetra_MultiVector& surf_p_c = *S_next_->GetFieldData("surface_pressure")
      ->ViewComponent("cell",false);

  for (int cs=0; cs!=ncells_surf; ++cs) {
    AmanziMesh::Entity_ID f = surf_mesh_->entity_get_parent(AmanziMesh::CELL, cs);
    surf_Pu_c[0][cs] = domain_Pu_f[0][f] - domain_p_f[0][f] + surf_p_c[0][cs];
  }

  // Update surface faces.
  // Calculate delta h on the surface
  Teuchos::RCP<CompositeVector> surf_Ph = Teuchos::rcp(new CompositeVector(*surf_Pu));
  surf_Ph->CreateData();
  surf_Ph->PutScalar(0.);

  // old ponded depth
  *surf_Ph->ViewComponent("cell",false) = *S_next_->GetFieldData("ponded_depth")->ViewComponent("cell",false);

  // new ponded depth
  S_next_->GetFieldData("surface_pressure",surf_pk_name_)
      ->ViewComponent("cell",false)->Update(-1., surf_Pu_c, 1.);
  surf_pk_->changed_solution();
  S_next_->GetFieldEvaluator("ponded_depth")
      ->HasFieldChanged(S_next_.ptr(), name_);

  // put delta ponded depth into surf_Ph_cell
  surf_Ph->ViewComponent("cell",false)
      ->Update(-1., *S_next_->GetFieldData("ponded_depth")->ViewComponent("cell",false), 1.);

  // update delta faces
  surf_preconditioner_->UpdateConsistentFaceCorrection(*surf_u, surf_Ph.ptr());
  *surf_Pu->ViewComponent("face",false) = *surf_Ph->ViewComponent("face",false);

  // revert solution so we don't break things
  S_next_->GetFieldData("surface_pressure",surf_pk_name_)
      ->ViewComponent("cell",false)->Update(1., surf_Pu_c, 1.);
  surf_pk_->changed_solution();

#if DEBUG_FLAG
  if (surf_c0_ < surf_u->size("cell",false) && surf_c1_ < surf_u->size("cell",false)) {
    //  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_HIGH, true)) {
     AmanziMesh::Entity_ID_List fnums,fnums0;
     std::vector<int> dirs;
     surf_mesh_->cell_get_faces_and_dirs(surf_c0_, &fnums0, &dirs);
     surf_mesh_->cell_get_faces_and_dirs(surf_c1_, &fnums, &dirs);

    *out_ << " Surface precon:" << std::endl;
    *out_ << "  u0: " << (*surf_u)("cell",surf_c0_) << ", "
          << (*surf_u)("face",fnums0[0]) << std::endl;
    *out_ << "  u1: " << (*surf_u)("cell",surf_c1_) << ", "
          << (*surf_u)("face",fnums[0]) << std::endl;
    *out_ << "  PC*u0: " << (*surf_Pu)("cell",surf_c0_) << ", "
          << (*surf_Pu)("face",fnums0[0]) << std::endl;
    *out_ << "  PC*u1: " << (*surf_Pu)("cell",surf_c1_) << ", "
          << (*surf_Pu)("face",fnums[0]) << std::endl;

  }
#endif
}

// updates the preconditioner
void MPCSurfaceSubsurfaceResidualCoupler::update_precon(double t,
        Teuchos::RCP<const TreeVector> up, double h) {
  MPCSurfaceSubsurfaceCoupler::update_precon(t, up, h);
  preconditioner_->AssembleGlobalMatrices();
  preconditioner_->ComputeSchurComplement(domain_pk_->bc_markers(),
          domain_pk_->bc_values());
  preconditioner_->UpdatePreconditioner();
}

// Modify predictor to ensure lambda and surface cell remain consistent
bool MPCSurfaceSubsurfaceResidualCoupler::modify_predictor(double h, Teuchos::RCP<TreeVector> u) {
  bool changed(false);

  Teuchos::RCP<TreeVector> surf_u = u->SubVector(surf_pk_name_);
  Teuchos::RCP<TreeVector> domain_u = u->SubVector(domain_pk_name_);

  if (modify_predictor_heuristic_) {
    const Epetra_MultiVector& surf_u_prev_c =
        *S_->GetFieldData("surface_pressure")->ViewComponent("cell",false);
    const double& patm = *S_next_->GetScalarData("atmospheric_pressure");

    Epetra_MultiVector& domain_u_f =
        *domain_u->data()->ViewComponent("face",false);
    Epetra_MultiVector& surf_u_c =
        *surf_u->data()->ViewComponent("cell",false);

    int ncells = surf_u_c.MyLength();
    for (int c=0; c!=ncells; ++c) {
      int f = surf_mesh_->entity_get_parent(AmanziMesh::CELL, c);

      double dp = surf_u_c[0][c] - surf_u_prev_c[0][c];
      double pnew = surf_u_c[0][c] - patm;
      double pold = surf_u_prev_c[0][c] - patm;

      if (pnew > 0) {
        if (dp > pnew) {
          std::cout << "CHANGING (first over?): p = " << surf_u_c[0][c] << " to " << patm + .001 << std::endl;
          surf_u_c[0][c] = patm + .001;
          domain_u_f[0][f] = surf_u_c[0][c];

        } else if (pold > 0 && dp > pold) {
          std::cout << "CHANGING (second over?): p = " << surf_u_c[0][c] << " to " << patm + 2*pold << std::endl;
          surf_u_c[0][c] = patm + 2*pold;
          domain_u_f[0][f] = surf_u_c[0][c];
        }
      }
    }
    changed = true;
  }


  changed |= surf_pk_->modify_predictor(h, surf_u);
  changed |= domain_pk_->modify_predictor(h, domain_u);

  return changed;
}

} // namespace

