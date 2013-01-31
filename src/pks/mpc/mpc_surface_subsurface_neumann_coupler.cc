/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Interface for the derived MPC for neumann coupling between surface and subsurface.
------------------------------------------------------------------------- */

#include "mpc_surface_subsurface_neumann_coupler.hh"

namespace Amanzi {

RegisteredPKFactory<MPCSurfaceSubsurfaceNeumannCoupler>
MPCSurfaceSubsurfaceNeumannCoupler::reg_("surface-subsurface Neumann coupler");


// -------------------------------------------------------------
// Special residual function evaluation process
// -------------------------------------------------------------
void MPCSurfaceSubsurfaceNeumannCoupler::fun(double t_old, double t_new,
        Teuchos::RCP<TreeVector> u_old, Teuchos::RCP<TreeVector> u_new,
        Teuchos::RCP<TreeVector> g) {
  // ensure the source term from subsurface -> surface is 0
  Teuchos::RCP<CompositeVector> source =
      S_next_->GetFieldData("overland_source_from_subsurface", surf_pk_name_);
  source->PutScalar(0.);

  // evaluate the residual for overland flow
  Teuchos::RCP<TreeVector> surf_u_new = u_new->SubVector(surf_pk_name_);
  Teuchos::RCP<TreeVector> surf_u_old = u_old->SubVector(surf_pk_name_);
  Teuchos::RCP<TreeVector> surf_g = g->SubVector(surf_pk_name_);
  surf_pk_->fun(t_old, t_new, surf_u_old, surf_u_new, surf_g);

  // Stuff the residual value into the source term.  This is the mass
  // imbalance on the surface, in mol/s, and is used in the subsurface PK.
  *source->ViewComponent("cell",false) = *surf_g->data()->ViewComponent("cell",false);

  // Since this is now the new source term, the residual is identically zero.
  surf_g->data()->ViewComponent("cell",false)->PutScalar(0.);

  // Evaluate the Richards residual, using this source as the Neumann BC.
  Teuchos::RCP<TreeVector> domain_u_new = u_new->SubVector(domain_pk_name_);
  Teuchos::RCP<TreeVector> domain_u_old = u_old->SubVector(domain_pk_name_);
  Teuchos::RCP<TreeVector> domain_g = g->SubVector(domain_pk_name_);
  domain_pk_->fun(t_old, t_new, Teuchos::null, domain_u_new, domain_g);

}


// applies preconditioner to u and returns the result in Pu
void MPCSurfaceSubsurfaceNeumannCoupler::precon(Teuchos::RCP<const TreeVector> u,
        Teuchos::RCP<TreeVector> Pu) {
  // Precondition evaluation is in reverse order.

  // Evaluate the Richards preconditioner
  Teuchos::RCP<const TreeVector> domain_u = u->SubVector(domain_pk_name_);
  Teuchos::RCP<TreeVector> domain_Pu = Pu->SubVector(domain_pk_name_);
  domain_pk_->precon(domain_u, domain_Pu);

  // Back substitute res_mass_surf <-- res_mass_surf - K/dz*dp
  // -- dQ_dp = K/dz
  const Epetra_MultiVector& dQ_dp =
      *S_next_->GetFieldData("doverland_source_from_subsurface_dsurface_pressure")
      ->ViewComponent("cell",false);

  // -- dp
  Epetra_MultiVector dp_limited(dQ_dp);
  const Epetra_MultiVector& dp = *domain_Pu->data()->ViewComponent("cell",false);

  AmanziMesh::Entity_ID_List cells;
  int ncells = dQ_dp.MyLength();
  for (int c=0; c!=ncells; ++c) {
    AmanziMesh::Entity_ID f = surf_mesh_->entity_get_parent(AmanziMesh::CELL, c);
    domain_mesh_->face_get_cells(f, AmanziMesh::OWNED, &cells);
    ASSERT(cells.size() == 1);
    dp_limited[0][c] = dp[0][cells[0]];
  }

  // -- res_mass_surf = res_mass_surf - 1.0 * K/dz * dp
  Teuchos::RCP<TreeVector> u_new =
      Teuchos::rcp(new TreeVector(*u->SubVector(surf_pk_name_)));
  *u_new = *u->SubVector(surf_pk_name_);

  u_new->data()->ViewComponent("cell",false)->Multiply(-1., dQ_dp, dp_limited, 1.0);

  // Evaluate the Overland preconditioner
  Teuchos::RCP<TreeVector> surf_Pu = Pu->SubVector(surf_pk_name_);
  surf_pk_->precon(u_new, surf_Pu);

}

} // namespace
