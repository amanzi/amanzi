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


MPCSurfaceSubsurfaceDirichletCoupler::MPCSurfaceSubsurfaceDirichletCoupler(
        Teuchos::ParameterList& plist, const Teuchos::RCP<TreeVector>& soln) :
    PKDefaultBase(plist, soln),
    StrongMPC(plist, soln) {
  surface_mesh_key_ = plist_.get<std::string>("surface mesh name", "surface");
  domain_field_ = plist_.get<std::string>("domain field name", "pressure");
  surface_field_ = plist_.get<std::string>("surface field name", "surface_pressure");
  domain_pk_name_ = plist_.get<std::string>("domain pk name", "flow");
  surface_pk_name_ = plist_.get<std::string>("surface pk name", "overland flow");
}


void MPCSurfaceSubsurfaceDirichletCoupler::fun(double t_old, double t_new,
        Teuchos::RCP<TreeVector> u_old, Teuchos::RCP<TreeVector> u_new,
        Teuchos::RCP<TreeVector> g) {
  StrongMPC::fun(t_old, t_new, u_old, u_new, g);

  // // copy over the residual
  // Teuchos::RCP<const AmanziMesh::Mesh> mesh = S_next_->GetMesh(surface_mesh_key_);
  // const Epetra_MultiVector& surf = *g->SubVector(surface_pk_name_)->data()
  //     ->ViewComponent("cell",false);
  // Epetra_MultiVector& domain = *g->SubVector(domain_pk_name_)->data()
  //     ->ViewComponent("face",false);

  // // copies surface cell values into subsurface face unknowns
  // int ncells = surf.MyLength();
  // for (int c=0; c!=ncells; ++c) {
  //   AmanziMesh::Entity_ID f = mesh->entity_get_parent(AmanziMesh::CELL, c);
  //   domain[0][f] = surf[0][c];
  // }

}


void MPCSurfaceSubsurfaceDirichletCoupler::changed_solution() {
  // Teuchos::RCP<const AmanziMesh::Mesh> mesh = S_next_->GetMesh(surface_mesh_key_);

  // const Epetra_MultiVector& surf = *S_next_->GetFieldData(surface_field_)
  //     ->ViewComponent("cell",false);
  // Epetra_MultiVector& domain = *S_next_->GetFieldData(domain_field_, domain_pk_name_)
  //     ->ViewComponent("face",false);

  // // copies surface cell values into subsurface face unknowns
  // int ncells = surf.MyLength();
  // for (int c=0; c!=ncells; ++c) {
  //   AmanziMesh::Entity_ID f = mesh->entity_get_parent(AmanziMesh::CELL, c);
  //   domain[0][f] = surf[0][c];
  // }

  StrongMPC::changed_solution();
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



// -----------------------------------------------------------------------------
// Update the preconditioner at time t and u = up
// -----------------------------------------------------------------------------
void MPCSurfaceSubsurfaceDirichletCoupler::update_precon(double t, Teuchos::RCP<const TreeVector> up, double h) {
  StrongMPC::update_precon(t,up,h);

  /*
  // test the preconditioner
  double eps = 1.e-3;
  Teuchos::RCP<TreeVector> u = Teuchos::rcp(new TreeVector(*up));
  Teuchos::RCP<TreeVector> fu = Teuchos::rcp(new TreeVector(*up));
  Teuchos::RCP<TreeVector> Pu = Teuchos::rcp(new TreeVector(*up));

  *u = *up;
  S_next_->SetData("pressure","flow",u->SubVector("flow")->data());
  S_next_->SetData("surface_pressure","overland flow",u->SubVector("overland flow")->data());
  changed_solution();

  fun(S_->time(), S_next_->time(), Teuchos::null, u, fu);

  double u_surf0 = (*u->SubVector("overland flow")->data())("cell",0);
  double fu_surf0 = (*fu->SubVector("overland flow")->data())("cell",0);

  double u_sub0 = (*u->SubVector("flow")->data())("cell",99);
  double fu_sub0 = (*fu->SubVector("flow")->data())("cell",99);

  double u_surf1 = u_surf0 + eps;
  (*u->SubVector("overland flow")->data())("cell",0) = u_surf1;
  changed_solution();
  fun(S_->time(), S_next_->time(),
          Teuchos::null, u, fu);
  double fu_surf1 = (*fu->SubVector("overland flow")->data())("cell",0);

  // apply the precon
  Teuchos::RCP<Flow::OverlandHeadFlow> pk_ol = Teuchos::rcp_dynamic_cast<Flow::OverlandHeadFlow>(sub_pks_[1]);
  ASSERT(pk_ol != Teuchos::null);
  fu->PutScalar(0.);
  Teuchos::RCP<CompositeVector> fu_surf = fu->SubVector("overland flow")->data();
  Teuchos::RCP<CompositeVector> u_surf = u->SubVector("overland flow")->data();
  (*fu_surf)("cell",0) = 1.;
  const Epetra_MultiVector& dh_dp =
    *S_next_->GetFieldData("dponded_depth_dsurface_pressure")->ViewComponent("cell",false);
  (*fu_surf)("cell",0) = dh_dp[0][0];
  u_surf->PutScalar(0.);
  pk_ol->preconditioner_->Apply(*fu_surf, u_surf.ptr());

  std::cout << "PRECON CHECK:" << std::endl;
  std::cout << "  fu_surf 0,1 = " << fu_surf0 << ", " << fu_surf1 << std::endl;
  std::cout << "  u_surf = " << u_surf0 << " deriv = " << (fu_surf1 - fu_surf0)/eps << ", " << (*u_surf)("cell",0) << std::endl;

  // now the subsurf
  double u_sub1 = u_sub0 + eps;
  (*u->SubVector("overland flow")->data())("cell",0) = u_surf0;
  (*u->SubVector("flow")->data())("cell",99) = u_sub1;
  changed_solution();
  fun(S_->time(), S_next_->time(),
          Teuchos::null, u, fu);
  double fu_sub1 = (*fu->SubVector("flow")->data())("cell",99);

  // apply the precon
  Teuchos::RCP<Flow::Richards> pk_ri = Teuchos::rcp_dynamic_cast<Flow::Richards>(sub_pks_[0]);
  ASSERT(pk_ri != Teuchos::null);
  fu->PutScalar(0.);
  Teuchos::RCP<CompositeVector> fu_sub = fu->SubVector("flow")->data();
  Teuchos::RCP<CompositeVector> u_sub = u->SubVector("flow")->data();
  (*fu_sub)("cell",99) = 1.;
  pk_ri->preconditioner_->Apply(*fu_sub, u_sub.ptr());
  std::cout << "  u_sub = " << u_sub0 << " deriv = " << (fu_sub1 - fu_sub0)/eps << ", " << (*u_sub)("cell",99) << std::endl;

  ASSERT(0);
  */

}

} // namespace
