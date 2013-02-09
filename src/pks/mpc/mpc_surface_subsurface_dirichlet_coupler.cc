/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Interface for the derived MPC for water coupling between surface and subsurface.
------------------------------------------------------------------------- */

#include "pk_physical_bdf_base.hh"
#include "mpc_surface_subsurface_dirichlet_coupler.hh"

namespace Amanzi {

RegisteredPKFactory<MPCSurfaceSubsurfaceDirichletCoupler>
MPCSurfaceSubsurfaceDirichletCoupler::reg_("surface-subsurface Dirichlet coupler");

// applies preconditioner to u and returns the result in Pu
void MPCSurfaceSubsurfaceDirichletCoupler::precon(Teuchos::RCP<const TreeVector> u,
        Teuchos::RCP<TreeVector> Pu) {
  // Evaluate the overland preconditioner (reverse order), pk 1
  Teuchos::RCP<const TreeVector> surf_u = u->SubVector(surf_pk_name_);
  Teuchos::RCP<TreeVector> surf_Pu = Pu->SubVector(surf_pk_name_);
  surf_pk_->precon(surf_u, surf_Pu);

  // The surface lambda values are Dirichlet, indicating that the residual u
  // is the mismatch to Dirichlet data.  This Dirichlet data is given by the
  // surface mesh's cell value, so the mismatch is the surface mesh's
  // correction Pu.
  // Add u_lambda to Pu_surface.
  Teuchos::RCP<const TreeVector> domain_u = u->SubVector(domain_pk_name_);
  Teuchos::RCP<TreeVector> new_domain_u = Teuchos::rcp(new TreeVector(*domain_u));
  *new_domain_u = *domain_u;

  Epetra_MultiVector& new_domain_u_f = *new_domain_u->data()->ViewComponent("face",false);
  const Epetra_MultiVector& surf_Pu_c = *surf_Pu->data()->ViewComponent("cell",false);

  int ncells = surf_Pu_c.MyLength();
  for (int c=0; c!=ncells; ++c) {
    AmanziMesh::Entity_ID f = surf_mesh_->entity_get_parent(AmanziMesh::CELL, c);
    new_domain_u_f[0][f] += surf_Pu_c[0][c];
  }

  // Call the Richards precon with the updated residual.
  Teuchos::RCP<TreeVector> domain_Pu = Pu->SubVector(domain_pk_name_);
  domain_pk_->precon(new_domain_u, domain_Pu);
}

} // namespace
