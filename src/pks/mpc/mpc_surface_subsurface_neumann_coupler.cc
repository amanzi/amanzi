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
// Initialize PKs
// -------------------------------------------------------------
void MPCSurfaceSubsurfaceNeumannCoupler::setup(const Teuchos::Ptr<State>& S) {
  StrongMPC::setup(S);

  // get the PKs
  if (sub_pks_[0]->name() == "flow") {
    pk_flow_ = sub_pks_[0];
    ASSERT(sub_pks_[1]->name() == "overland flow");
    pk_ol_ = sub_pks_[1];
  } else {
    ASSERT(sub_pks_[0]->name() == "overland flow");
    ASSERT(sub_pks_[1]->name() == "flow");
    pk_flow_ = sub_pks_[1];
    pk_ol_ = sub_pks_[0];
  }

}


// -------------------------------------------------------------
// Special residual function evaluation process
// -------------------------------------------------------------
void MPCSurfaceSubsurfaceNeumannCoupler::fun(double t_old, double t_new,
        Teuchos::RCP<TreeVector> u_old, Teuchos::RCP<TreeVector> u_new,
        Teuchos::RCP<TreeVector> g) {
  // ensure the source term from subsurface -> surface is 0
  Teuchos::RCP<CompositeVector> source =
      S_next_->GetFieldData("overland_source_from_subsurface", pk_ol_->name());
  source->PutScalar(0.);

  // evaluate the residual for overland flow
  Teuchos::RCP<TreeVector> overland_u_new = u_new->SubVector(pk_ol_->name());
  Teuchos::RCP<TreeVector> overland_u_old = u_old->SubVector(pk_ol_->name());
  Teuchos::RCP<TreeVector> overland_g = g->SubVector(pk_ol_->name());
  pk_ol_->fun(t_old, t_new, overland_u_old, overland_u_new, overland_g);

  // Stuff the residual value into the source term.  This is the mass
  // imbalance on the surface, in mol/s, and is used in the subsurface PK.
  *source->ViewComponent("cell",false) = *overland_g->data()->ViewComponent("cell",false);

  // Since this is now the new source term, the residual is identically zero.
  overland_g->data()->ViewComponent("cell",false)->PutScalar(0.);

  // Evaluate the Richards residual, using this source as the Neumann BC.
  Teuchos::RCP<TreeVector> richards_u_new = u_new->SubVector(pk_flow_->name());
  Teuchos::RCP<TreeVector> richards_u_old = u_old->SubVector(pk_flow_->name());
  Teuchos::RCP<TreeVector> richards_g = g->SubVector(pk_flow_->name());
  pk_flow_->fun(t_old, t_new, Teuchos::null, richards_u_new, richards_g);

}


// applies preconditioner to u and returns the result in Pu
void MPCSurfaceSubsurfaceNeumannCoupler::precon(Teuchos::RCP<const TreeVector> u,
        Teuchos::RCP<TreeVector> Pu) {
  // Precondition evaluation is in reverse order.

  // Evaluate the Richards preconditioner
  Teuchos::RCP<const TreeVector> richards_u = u->SubVector(pk_flow_->name());
  Teuchos::RCP<TreeVector> richards_Pu = Pu->SubVector(pk_flow_->name());
  pk_flow_->precon(richards_u, richards_Pu);

  // Back substitute res_mass_surf <-- res_mass_surf - K/dz*dp
  // -- dQ_dp = K/dz
  const Epetra_MultiVector& dQ_dp =
      *S_next_->GetFieldData("doverland_source_from_subsurface_dsurface_pressure")
      ->ViewComponent("cell",false);

  // -- dp
  Epetra_MultiVector dp_limited(dQ_dp);
  const Epetra_MultiVector& dp = *richards_Pu->data()->ViewComponent("cell",false);
  Teuchos::RCP<const AmanziMesh::Mesh> surface = S_next_->GetMesh("surface");
  Teuchos::RCP<const AmanziMesh::Mesh> subsurf = S_next_->GetMesh();

  AmanziMesh::Entity_ID_List cells;
  int ncells = dQ_dp.MyLength();
  for (int c=0; c!=ncells; ++c) {
    AmanziMesh::Entity_ID f = surface->entity_get_parent(AmanziMesh::CELL, c);
    subsurf->face_get_cells(f, AmanziMesh::OWNED, &cells);
    ASSERT(cells.size() == 1);
    dp_limited[0][c] = dp[0][cells[0]];
  }

  // -- res_mass_surf = res_mass_surf - 1.0 * K/dz * dp
  Teuchos::RCP<TreeVector> u_new =
      Teuchos::rcp(new TreeVector(*u->SubVector(pk_ol_->name())));
  *u_new = *u->SubVector(pk_ol_->name());

  u_new->data()->ViewComponent("cell",false)->Multiply(-1., dQ_dp, dp_limited, 1.0);

  // Evaluate the Overland preconditioner
  Teuchos::RCP<TreeVector> ol_Pu = Pu->SubVector(pk_ol_->name());
  pk_ol_->precon(u_new, ol_Pu);

}

} // namespace
