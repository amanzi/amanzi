/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Interface for the derived MPC for water coupling between surface and subsurface.

To be used with either MPCSurfaceSubsurfaceDirichletCoupler or
MPCSurfaceSubsurfaceFluxCoupler.

------------------------------------------------------------------------- */

#include "mpc_water_flux_coupler.hh"

namespace Amanzi {

void
MPCWaterFluxCoupler::fun(double t_old, double t_new, Teuchos::RCP<TreeVector> u_old,
                         Teuchos::RCP<TreeVector> u_new, Teuchos::RCP<TreeVector> g) {
  // propagate updated info into state
  solution_to_state(u_new, S_next_);

  // evaluate the residual of the surface equation
  Teuchos::RCP<TreeVector> surf_u_new = u_new->SubVector(surf_pk_index_);
  Teuchos::RCP<TreeVector> surf_g = g->SubVector(surf_pk_index_);
  surf_pk_->fun(t_old, t_new, Teuchos::null, surf_u_new, surf_g);

  // The residual of the surface equation provides the flux.  This is the mass
  // imbalance on the surface, and is used in the subsurface PK.
  Epetra_MultiVector& source = *S_next_->GetFieldData(flux_key_, domain_pk_name_)
      ->ViewComponent("cell",false);
  source = *surf_g->Data()->ViewComponent("cell",false);

  // The exception to this is if the surface unfrozen fraction is 0, in which
  // case the rel perm will be zero, and no flux is allowed.
  Epetra_MultiVector& surf_g_c = *surf_g->Data()->ViewComponent("cell",false);
  const Epetra_MultiVector& uf_frac = *S_next_->GetFieldData("unfrozen_fraction")
      ->ViewComponent("cell",false);
  for (unsigned int sc=0; sc!=uf_frac.MyLength(); ++sc) {
    if (uf_frac[0][sc] == 0.) {
      source[0][sc] = 0.;
    }
  }

  // Evaluate the subsurface residual, which uses this flux as a Neumann BC.
  Teuchos::RCP<TreeVector> domain_u_new = u_new->SubVector(domain_pk_index_);
  Teuchos::RCP<TreeVector> domain_g = g->SubVector(domain_pk_index_);
  domain_pk_->fun(t_old, t_new, Teuchos::null, domain_u_new, domain_g);

  // Clobber the subsurface face's residual, as it gets hit with zero rel perm
  Epetra_MultiVector& domain_g_f = *domain_g->Data()->ViewComponent("face",false);
  for (unsigned int sc=0; sc!=uf_frac.MyLength(); ++sc) {
    if (uf_frac[0][sc] == 0.) {
      AmanziMesh::Entity_ID f = surf_mesh_->entity_get_parent(AmanziMesh::CELL, sc);
      domain_g_f[0][f] = 0.;
    } else { // residual taken by subsurface flux
      surf_g_c[0][sc] = 0.;
    }
  }
}

// bool MPCWaterFluxCoupler::PreconPostprocess_(Teuchos::RCP<const TreeVector> res,
//         Teuchos::RCP<TreeVector> Pu) {
//   bool modified = MPCWaterCoupler<MPCSurfaceSubsurfaceFluxCoupler>::PreconPostprocess_(res,Pu);

//   // check for ice freezing
//   int n_modified_l = 0;
//   const Epetra_MultiVector& domain_p_f = *S_next_->GetFieldData("pressure")
//       ->ViewComponent("face",false);
//   const Epetra_MultiVector& surf_g_c = *res->SubVector(surf_pk_index_)->Data()->ViewComponent("cell",false);
//   Epetra_MultiVector& domain_du_f = *Pu->SubVector(domain_pk_index_)->Data()->ViewComponent("face",false);
//   for (unsigned int sc=0; sc!=surf_g_c.MyLength(); ++sc) {
//     AmanziMesh::Entity_ID f = surf_mesh_->entity_get_parent(AmanziMesh::CELL,sc);
//     if (surf_g_c[0][sc] < 0. && (domain_p_f[0][f] - domain_du_f[0][f] < 101325.)) {
//       domain_du_f[0][f] = domain_p_f[0][f] - (101325. + 100.);
//       n_modified_l++;
//     }
//   }

//   int n_modified;
//   this->domain_mesh_->get_comm()->SumAll(&n_modified_l, &n_modified, 1);
//   if (n_modified > 0) modified = true;
//   return modified;
// }


} // namespace

