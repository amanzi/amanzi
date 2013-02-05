/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Interface for the derived MPC for water coupling between surface and subsurface.
------------------------------------------------------------------------- */

#include "mpc_surface_subsurface_coupler.hh"

namespace Amanzi {

RegisteredPKFactory<MPCSurfaceSubsurfaceCoupler>
MPCSurfaceSubsurfaceCoupler::reg_("surface-subsurface coupler");


MPCSurfaceSubsurfaceCoupler::MPCSurfaceSubsurfaceCoupler(
        Teuchos::ParameterList& plist, const Teuchos::RCP<TreeVector>& soln) :
    PKDefaultBase(plist, soln),
    StrongMPC(plist, soln) {
  surface_mesh_key_ = plist_.get<std::string>("surface mesh name", "surface");
  domain_field_ = plist_.get<std::string>("domain field name", "pressure");
  surface_field_ = plist_.get<std::string>("surface field name", "surface_pressure");
  domain_pk_name_ = plist_.get<std::string>("domain pk name", "flow");
  surface_pk_name_ = plist_.get<std::string>("surface pk name", "overland flow");
}


void MPCSurfaceSubsurfaceCoupler::fun(double t_old, double t_new,
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


void MPCSurfaceSubsurfaceCoupler::changed_solution() {
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
void MPCSurfaceSubsurfaceCoupler::precon(Teuchos::RCP<const TreeVector> u,
        Teuchos::RCP<TreeVector> Pu) {
  StrongMPC::precon(u,Pu);

  // overwrite the subsurface's surface face updates with that from the surface.
  const Epetra_MultiVector& Pu_surf_c = *Pu->SubVector("overland flow")
      ->data()->ViewComponent("cell",false);
  Epetra_MultiVector& Pu_sub_f = *Pu->SubVector("flow")
      ->data()->ViewComponent("face",false);

  Teuchos::RCP<const AmanziMesh::Mesh> surface = S_next_->GetMesh("surface");
  int ncells = Pu_surf_c.MyLength();
  for (int c=0; c!=ncells; ++c) {
    // -- get the surface cell's equivalent subsurface face and neighboring cell
    AmanziMesh::Entity_ID f =
        surface->entity_get_parent(AmanziMesh::CELL, c);
    Pu_sub_f[0][f] = Pu_surf_c[0][c];
  }
}


} // namespace
