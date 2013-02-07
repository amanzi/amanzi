/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Interface for the derived MPC for water coupling between surface and subsurface.
------------------------------------------------------------------------- */

#include "overland_head.hh"
#include "richards.hh"
#include "mpc_surface_subsurface_dirichlet_coupler.hh"

namespace Amanzi {

RegisteredPKFactory<MPCSurfaceSubsurfaceDirichletCoupler>
MPCSurfaceSubsurfaceDirichletCoupler::reg_("surface-subsurface Dirichlet coupler");



void MPCSurfaceSubsurfaceDirichletCoupler::setup(const Teuchos::Ptr<State>& S) {
  // grab the meshes
  std::string sub_mesh_key = plist.get<std::string>("subsurface mesh key", "domain");
  sub_mesh_ = S->GetMesh(sub_mesh_key);
  std::string surf_mesh_key = plist.get<std::string>("surface mesh key", "surface");
  surf_mesh_ = S->GetMesh(surf_mesh_key);

  // grab the PKs

}


// applies preconditioner to u and returns the result in Pu
void MPCSurfaceSubsurfaceDirichletCoupler::precon(Teuchos::RCP<const TreeVector> u,
        Teuchos::RCP<TreeVector> Pu) {
  // Evaluate the overland preconditioner (reverse order), pk 1
  Teuchos::RCP<PKBDFBase> ol_pk = sub_pks_[1];
  Teuchos::RCP<const TreeVector> ol_u = u->SubVector(ol_pk->name());
  Teuchos::RCP<TreeVector> ol_Pu = Pu->SubVector(ol_pk->name());
  ol_pk->precon(ol_u, ol_Pu);

  // The surface lambda values are Dirichlet, indicating that the residual u
  // is the mismatch to Dirichlet data.  This Dirichlet data is given by the
  // surface mesh's cell value, so the mismatch is the surface mesh's
  // correction Pu.
  // Add u_lambda to Pu_surface.
  Teuchos::RCP<PKBDFBase> rich_pk = sub_pks_[0];
  Teuchos::RCP<const TreeVector> rich_u = u->SubVector(rich_pk->name());
  Teuchos::RCP<TreeVector> new_rich_u = Teuchos::rcp(new TreeVector(*rich_u));
  *new_rich_u = *rich_u;

  Epetra_MultiVector& new_rich_u_f = *new_rich_u->data()->ViewComponent("face",false);
  const Epetra_MultiVector& ol_Pu_c = *ol_Pu->data()->ViewComponent("cell",false);

  Teuchos::RCP<const AmanziMesh::Mesh> surface = S_next_->GetMesh("surface");
  int ncells = ol_Pu_c.MyLength();
  for (int c=0; c!=ncells; ++c) {
    AmanziMesh::Entity_ID f = surface->entity_get_parent(AmanziMesh::CELL, c);
    new_rich_u_f[0][f] += ol_Pu_c[0][c];
  }

  // Call the Richards precon with the updated residual.
  Teuchos::RCP<TreeVector> rich_Pu = Pu->SubVector(rich_pk->name());
  rich_pk->precon(new_rich_u, rich_Pu);
}

} // namespace
