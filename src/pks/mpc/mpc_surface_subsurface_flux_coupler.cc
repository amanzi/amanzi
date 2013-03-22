/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Interface for the derived MPC for water coupling between surface and subsurface.

In this method, a Dirichlet BC is used on the subsurface boundary for
the operator, but the preconditioner is for the flux system with no
extra unknowns.  On the surface, the TPFA is used, resulting in a
subsurface-face-only Schur complement that captures all terms.

------------------------------------------------------------------------- */
#include "EpetraExt_RowMatrixOut.h"

#include "matrix_mfd_surf.hh"
#include "matrix_mfd_tpfa.hh"
#include "pk_physical_bdf_base.hh"

#include "mpc_surface_subsurface_flux_coupler.hh"


namespace Amanzi {

#define DEBUG_FLAG 1

RegisteredPKFactory<MPCSurfaceSubsurfaceFluxCoupler>
MPCSurfaceSubsurfaceFluxCoupler::reg_("surface-subsurface flux coupler");

// -- Setup data.
void MPCSurfaceSubsurfaceFluxCoupler::setup(const Teuchos::Ptr<State>& S) {
  MPCSurfaceSubsurfaceCoupler::setup(S);

  // get the flux key
  flux_key_ = plist_.get<std::string>("flux key");

  // Get the domain's preconditioner and replace it with a MatrixMFD_Surf.
  Teuchos::RCP<Operators::Matrix> precon = domain_pk_->preconditioner();
  Teuchos::RCP<Operators::MatrixMFD> mfd_precon =
    Teuchos::rcp_dynamic_cast<Operators::MatrixMFD>(precon);
  ASSERT(mfd_precon != Teuchos::null);

  mfd_preconditioner_ = Teuchos::rcp(new Operators::MatrixMFD_Surf(*mfd_precon, surf_mesh_));
  preconditioner_ = mfd_preconditioner_;

  // Get the surface's preconditioner and ensure it is TPFA.
  Teuchos::RCP<Operators::Matrix> surf_precon = surf_pk_->preconditioner();
  surf_preconditioner_ =
    Teuchos::rcp_dynamic_cast<Operators::MatrixMFD_TPFA>(surf_precon);
  ASSERT(surf_preconditioner_ != Teuchos::null);

  // set the surface A in the MFD_Surf.
  mfd_preconditioner_->set_surface_A(surf_preconditioner_);

  // give the PCs back to the PKs
  domain_pk_->set_preconditioner(mfd_preconditioner_);
  surf_pk_->set_preconditioner(surf_preconditioner_);
}


// -------------------------------------------------------------
// Special function evaluation process
// -------------------------------------------------------------
void MPCSurfaceSubsurfaceFluxCoupler::fun(double t_old, double t_new,
        Teuchos::RCP<TreeVector> u_old, Teuchos::RCP<TreeVector> u_new,
        Teuchos::RCP<TreeVector> g) {

  // evaluate the residual of the surface equation
  Teuchos::RCP<TreeVector> surf_u_new = u_new->SubVector(surf_pk_name_);
  Teuchos::RCP<TreeVector> surf_g = g->SubVector(surf_pk_name_);
  surf_pk_->fun(t_old, t_new, Teuchos::null, surf_u_new, surf_g);

  // The residual of the surface equation provides the flux.  This is the mass
  // imbalance on the surface, and is used in the subsurface PK.
  Teuchos::RCP<CompositeVector> source =
      S_next_->GetFieldData(flux_key_, domain_pk_name_);
  *source->ViewComponent("cell",false) = *surf_g->data()->ViewComponent("cell",false);

  // Evaluate the subsurface residual, which uses this flux as a Neumann BC.
  Teuchos::RCP<TreeVector> domain_u_new = u_new->SubVector(domain_pk_name_);
  Teuchos::RCP<TreeVector> domain_g = g->SubVector(domain_pk_name_);
  domain_pk_->fun(t_old, t_new, Teuchos::null, domain_u_new, domain_g);

  // Set the surf cell residual to 0.
  surf_g->data()->ViewComponent("cell",false)->PutScalar(0.);
}


// -------------------------------------------------------------
// Apply preconditioner to u and returns the result in Pu
// -------------------------------------------------------------
void MPCSurfaceSubsurfaceFluxCoupler::precon(Teuchos::RCP<const TreeVector> u,
        Teuchos::RCP<TreeVector> Pu) {
  // Apply the combined preconditioner to the subsurface residual
  PreconApply_(u,Pu);

  // Damp, kluge, hack, etc.
  PreconPostprocess_(u,Pu);

  // Update surface values.
  PreconUpdateSurfaceCells_(u,Pu);
  PreconUpdateSurfaceFaces_(u,Pu);

#if DEBUG_FLAG
  Teuchos::OSTab tab = getOSTab();
  Teuchos::RCP<const CompositeVector> surf_u = u->SubVector(surf_pk_name_)->data();
  Teuchos::RCP<const CompositeVector> domain_u = u->SubVector(domain_pk_name_)->data();
  Teuchos::RCP<CompositeVector> domain_Pu = Pu->SubVector(domain_pk_name_)->data();
  Teuchos::RCP<CompositeVector> surf_Pu = Pu->SubVector(surf_pk_name_)->data();

  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_HIGH, true)) {
    AmanziMesh::Entity_ID_List fnums1,fnums0;
    std::vector<int> dirs;
    domain_mesh_->cell_get_faces_and_dirs(c0_, &fnums0, &dirs);
    domain_mesh_->cell_get_faces_and_dirs(c1_, &fnums1, &dirs);

    *out_ << "Preconditioner application" << std::endl;
    *out_ << " SubSurface precon:" << std::endl;
    *out_ << "  p0: " << (*domain_u)("cell",c0_);
    for (int n=0; n!=fnums0.size(); ++n) *out_ << ", " << (*domain_u)("face",fnums0[n]);
    *out_ << std::endl;
    *out_ << "  p1: " << (*domain_u)("cell",c1_);
    for (int n=0; n!=fnums1.size(); ++n) *out_ << ", " << (*domain_u)("face",fnums1[n]);
    *out_ << std::endl;
    *out_ << "  PC*p0: " << (*domain_Pu)("cell",c0_);
    for (int n=0; n!=fnums0.size(); ++n) *out_ << ", " << (*domain_Pu)("face",fnums0[n]);
    *out_ << std::endl;
    *out_ << "  PC*p1: " << (*domain_Pu)("cell",c1_);
    for (int n=0; n!=fnums1.size(); ++n) *out_ << ", " << (*domain_Pu)("face",fnums1[n]);
    *out_ << std::endl;
  }

  if (surf_c0_ < surf_u->size("cell",false) && surf_c1_ < surf_u->size("cell",false)) {
    //  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_HIGH, true)) {
     AmanziMesh::Entity_ID_List fnums1,fnums0;
     std::vector<int> dirs;
     surf_mesh_->cell_get_faces_and_dirs(surf_c0_, &fnums0, &dirs);
     surf_mesh_->cell_get_faces_and_dirs(surf_c1_, &fnums1, &dirs);

    *out_ << " Surface precon:" << std::endl;
    *out_ << "  u0: " << (*surf_u)("cell",surf_c0_) << ", "
          << (*surf_u)("face",fnums0[0]) << std::endl;
    *out_ << "  u1: " << (*surf_u)("cell",surf_c1_) << ", "
          << (*surf_u)("face",fnums1[0]) << std::endl;
    *out_ << "  PC*u0: " << (*surf_Pu)("cell",surf_c0_) << ", "
          << (*surf_Pu)("face",fnums0[0]) << std::endl;
    *out_ << "  PC*u1: " << (*surf_Pu)("cell",surf_c1_) << ", "
          << (*surf_Pu)("face",fnums1[0]) << std::endl;

  }
#endif
}


void MPCSurfaceSubsurfaceFluxCoupler::PreconApply_(
    Teuchos::RCP<const TreeVector> u,
    Teuchos::RCP<TreeVector> Pu) {
  // preconditioner_->ApplyInverse(*u->SubVector(domain_pk_name_),
  //         Pu->SubVector(domain_pk_name_).ptr());
  domain_pk_->precon(u->SubVector(domain_pk_name_), Pu->SubVector(domain_pk_name_));
};


// -------------------------------------------------------------
// Post-processing in the preconditioner takes the corrections from the
// surface cell's parent faces and uses them for the surface cell, and then
// calculates a correction for the surface faces.
// -------------------------------------------------------------
void MPCSurfaceSubsurfaceFluxCoupler::PreconUpdateSurfaceCells_(
    Teuchos::RCP<const TreeVector> u,
    Teuchos::RCP<TreeVector> Pu) {

  Teuchos::RCP<CompositeVector> domain_Pu = Pu->SubVector(domain_pk_name_)->data();
  Teuchos::RCP<CompositeVector> surf_Pu = Pu->SubVector(surf_pk_name_)->data();

  int ncells_surf = surf_mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);

  // Correction applies to both the domain face and the surface cell.  Also
  // correct for drift.
  const Epetra_MultiVector& domain_Pu_f = *domain_Pu->ViewComponent("face",false);
  const Epetra_MultiVector& domain_p_f = *S_next_->GetFieldData("pressure")
      ->ViewComponent("face",false);
  Epetra_MultiVector& surf_Pu_c = *surf_Pu->ViewComponent("cell",false);
  const Epetra_MultiVector& surf_p_c = *S_next_->GetFieldData("surface_pressure")
      ->ViewComponent("cell",false);

  for (int cs=0; cs!=ncells_surf; ++cs) {
    AmanziMesh::Entity_ID f = surf_mesh_->entity_get_parent(AmanziMesh::CELL, cs);
    surf_Pu_c[0][cs] = domain_Pu_f[0][f] - domain_p_f[0][f] + surf_p_c[0][cs];
  }
}


// -------------------------------------------------------------
// Post-processing in the preconditioner takes the corrections from the
// surface cell's parent faces and uses them for the surface cell, and then
// calculates a correction for the surface faces.
// -------------------------------------------------------------
void MPCSurfaceSubsurfaceFluxCoupler::PreconUpdateSurfaceFaces_(
    Teuchos::RCP<const TreeVector> u,
    Teuchos::RCP<TreeVector> Pu) {
  // update delta faces
  Teuchos::RCP<const CompositeVector> surf_u = u->SubVector(surf_pk_name_)->data();
  Teuchos::RCP<CompositeVector> surf_Pu = Pu->SubVector(surf_pk_name_)->data();
  surf_preconditioner_->UpdateConsistentFaceCorrection(*surf_u, surf_Pu.ptr());
}


// updates the preconditioner
void MPCSurfaceSubsurfaceFluxCoupler::update_precon(double t,
        Teuchos::RCP<const TreeVector> up, double h) {
  MPCSurfaceSubsurfaceCoupler::update_precon(t, up, h);
  mfd_preconditioner_->AssembleGlobalMatrices();
  mfd_preconditioner_->ComputeSchurComplement(domain_pk_->bc_markers(),
          domain_pk_->bc_values());
  mfd_preconditioner_->UpdatePreconditioner();
}


// -----------------------------------------------------------------------------
// Modify the nonlinear predictor to avoid issues.
// This is very delicate code to get the order correct...
// -----------------------------------------------------------------------------
bool MPCSurfaceSubsurfaceFluxCoupler::modify_predictor(double h,
        Teuchos::RCP<TreeVector> up) {
  bool changed(false);

  if (modify_predictor_flux_bc_) {
    // To get this to work, we want to first modify the surface predictor,
    // then call the surface residual function to get the surface flux BCs,
    // then call the subsurface modify_predictor.  The subsurface modify
    // routine will alter the surface face lambdas to make them reasonably
    // consistent with the flux BCs.

    // -- call surface's modify_predictor()
    Teuchos::RCP<TreeVector> surf_u = up->SubVector(surf_pk_name_);
    changed |= surf_pk_->modify_predictor(h, surf_u);

    // -- call the surface residual function to get the surface flux BCs
    Teuchos::RCP<TreeVector> surf_g = Teuchos::rcp(new TreeVector(*surf_u));
    surf_pk_->fun(S_inter_->time(), S_next_->time(), Teuchos::null,
                         surf_u, surf_g);

    // -- set the flux BCs
    Teuchos::RCP<CompositeVector> source =
        S_next_->GetFieldData(flux_key_, domain_pk_name_);
    *source->ViewComponent("cell",false) = *surf_g->data()->ViewComponent("cell",false);

    // -- call the subsurface modify_predictor(), which uses those BCs
    Teuchos::RCP<TreeVector> domain_u = up->SubVector(domain_pk_name_);
    changed |= domain_pk_->modify_predictor(h, domain_u);

    // -- ensure the surface cell's predictor matches the subsurface face
    //    values (this is backwards from MPCSurfaceSubsurfaceCoupler's
    //    version.
    if (changed) {
      const Epetra_MultiVector& domain_u_f = *domain_u->data()
          ->ViewComponent("face",false);
      Epetra_MultiVector& surf_u_c = *surf_u->data()->ViewComponent("cell",false);

      int ncells = surf_u_c.MyLength();
      for (int c=0; c!=ncells; ++c) {
        int f = surf_mesh_->entity_get_parent(AmanziMesh::CELL, c);
        surf_u_c[0][c] = domain_u_f[0][f];
      }
    }

  } else {
    changed |= MPCSurfaceSubsurfaceCoupler::modify_predictor(h, up);
  }

  return changed;
}


} // namespace

