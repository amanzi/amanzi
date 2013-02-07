/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Interface for the derived MPC for water coupling between surface and subsurface.

In this method, a Dirichlet BC is used on the subsurface boundary for
the operator, but the preconditioner is for the full system with no
extra unknowns.  On the surface, the TPFA is used, resulting in a
subsurface-face-only Schur complement that captures all terms.

------------------------------------------------------------------------- */

#include "matrix_mfd_surf.hh"
#include "matrix_mfd_tpfa.hh"
#include "pk_physical_bdf_base.hh"

#include "mpc_surface_subsurface_full_coupler.hh"


namespace Amanzi {

#define DEBUG_FLAG 1

// -- Setup data.
void MPCSurfaceSubsurfaceFullCoupler::setup(const Teuchos::Ptr<State>& S) {
  MPCSurfaceSubsurfaceCoupler::setup(S);

  // Get the domain's preconditioner and replace it with a MatrixMFD_Surf.
  Teuchos::RCP<Operators::Matrix> precon = domain_pk_->preconditioner();
  Teuchos::RCP<Operators::MatrixMFD> mfd_precon =
    Teuchos::rcp_dynamic_cast<Operators::MatrixMFD>(precon);
  ASSERT(mfd_precon != Teuchos::null);

  preconditioner_ = Teuchos::rcp(new Operators::MatrixMFD_Surf(*mfd_precon, surf_mesh_));

  // Get the surface's preconditioner and replace it with a TPFA matrix.
  Teuchos::RCP<Operators::Matrix> surf_precon = surf_pk_->preconditioner();
  Teuchos::RCP<Operators::MatrixMFD> mfd_surf_precon =
    Teuchos::rcp_dynamic_cast<Operators::MatrixMFD>(surf_precon);
  ASSERT(mfd_surf_precon != Teuchos::null);

  Teuchos::RCP<Operators::MatrixMFD_TPFA> surf_preconditioner_ =
    Teuchos::rcp(new Operators::MatrixMFD_TPFA(*mfd_surf_precon));

  // set the surface A in the MFD_Surf.
  preconditioner_->set_surface_A(surf_preconditioner_);

  // give the PCs back to the PKs
  domain_pk_->set_preconditioner(preconditioner_);
  surf_pk_->set_preconditioner(surf_preconditioner_);

}

// applies preconditioner to u and returns the result in Pu
void MPCSurfaceSubsurfaceFullCoupler::precon(Teuchos::RCP<const TreeVector> u,
        Teuchos::RCP<TreeVector> Pu) {
  // Grab the surface cell residuals and stick them in the corresponding
  // subsurface face locations.
  Teuchos::RCP<const CompositeVector> domain_u = u->SubVector(domain_pk_name_)->data();
  Teuchos::RCP<CompositeVector> domain_u_new = Teuchos::rcp(new CompositeVector(*domain_u));
  domain_u_new->CreateData();
  *domain_u_new = *domain_u;

  const Epetra_MultiVector& surf_u_c = *u->SubVector(surf_pk_name_)
      ->data()->ViewComponent("cell",false);
  Epetra_MultiVector& domain_u_new_f = *domain_u_new->ViewComponent("face",false);

  int ncells_surf = surf_u_c.MyLength();
  for (int cs=0; cs!=ncells_surf; ++cs) {
    AmanziMesh::Entity_ID f = surf_mesh_->entity_get_parent(AmanziMesh::CELL, cs);
    domain_u_new_f[0][f] = surf_u_c[0][cs];
  }

  // Apply the combined preconditioner
  Teuchos::RCP<CompositeVector> domain_Pu = Pu->SubVector(domain_pk_name_)->data();
  preconditioner_->ApplyInverse(*domain_u_new, domain_Pu.ptr());

#if DEBUG_FLAG
  Teuchos::OSTab tab = getOSTab();
  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_HIGH, true)) {
    *out_ << "Preconditioner application" << std::endl;
    *out_ << " SubSurface precon:" << std::endl;
    *out_ << "  p0: " << (*domain_u_new)("cell",0) << " " << (*domain_u_new)("face",3)
          << std::endl;
    *out_ << "  p1: " << (*domain_u_new)("cell",99) << " " << (*domain_u_new)("face",500)
          << std::endl;
    *out_ << "  PC*p0: " << (*domain_Pu)("cell",0) << " " << (*domain_Pu)("face",3)
          << std::endl;
    *out_ << "  PC*p1: " << (*domain_Pu)("cell",99) << " " << (*domain_Pu)("face",500)
          << std::endl;
  }
#endif

  // Correction applies to both the domain face and the surface cell.
  const Epetra_MultiVector& domain_Pu_f = *domain_Pu->ViewComponent("face",false);
  Epetra_MultiVector& surf_Pu_c = *Pu->SubVector(surf_pk_name_)->data()
      ->ViewComponent("cell",false);

  for (int cs=0; cs!=ncells_surf; ++cs) {
    AmanziMesh::Entity_ID f = surf_mesh_->entity_get_parent(AmanziMesh::CELL, cs);
    surf_Pu_c[0][cs] = domain_Pu_f[0][f];
  }

  // Update surface with consistent faces
  Teuchos::RCP<CompositeVector> surf_Pu = Pu->SubVector(surf_pk_name_)->data();
  surf_Pu->Update(1., *S_next_->GetFieldData("surface_pressure"), -1.);
  surf_preconditioner_->UpdateConsistentFaceConstraints(surf_Pu.ptr());
  surf_Pu->Update(1., *S_next_->GetFieldData("surface_pressure"), -1.);

#if DEBUG_FLAG
  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_HIGH, true)) {
    Teuchos::RCP<const CompositeVector> surf_u = Pu->SubVector(surf_pk_name_)->data();
    *out_ << " Surface precon:" << std::endl;
    *out_ << "  p0: " << (*surf_u)("cell",0) << std::endl;
    *out_ << "  PC*p0: " << (*surf_Pu)("cell",0) << std::endl;
  }
#endif
}

// updates the preconditioner
void MPCSurfaceSubsurfaceFullCoupler::update_precon(double t,
        Teuchos::RCP<const TreeVector> up, double h) {
  MPCSurfaceSubsurfaceCoupler::update_precon(t, up, h);
  preconditioner_->AssembleGlobalMatrices();
  preconditioner_->ComputeSchurComplement(domain_pk_->bc_markers(),
          domain_pk_->bc_values());
  preconditioner_->UpdatePreconditioner();
}


} // namespace

