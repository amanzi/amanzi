/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Interface for the derived MPC for water coupling between surface and subsurface.
------------------------------------------------------------------------- */

#include "mpc_surface_subsurface_water_coupler.hh"

namespace Amanzi {

RegisteredPKFactory<MPCSurfaceSubsurfaceWaterCoupler>
MPCSurfaceSubsurfaceWaterCoupler::reg_("surface-subsurface water coupler");


// -------------------------------------------------------------
// Initialize PKs
// -------------------------------------------------------------
void MPCSurfaceSubsurfaceWaterCoupler::setup(const Teuchos::Ptr<State>& S) {
  StrongMPC::setup(S);

  // check that the order is as expected
  ASSERT(sub_pks_[1]->name() == std::string("flow"));
  ASSERT(sub_pks_[0]->name() == std::string("overland flow"));
}


// -------------------------------------------------------------
// Special residual function evaluation process
// -------------------------------------------------------------
void MPCSurfaceSubsurfaceWaterCoupler::fun(double t_old, double t_new,
        Teuchos::RCP<TreeVector> u_old, Teuchos::RCP<TreeVector> u_new,
        Teuchos::RCP<TreeVector> g) {
  // get the respective PKs
  Teuchos::RCP<PKBDFBase> richards_pk = sub_pks_[1];
  Teuchos::RCP<PKBDFBase> overland_pk = sub_pks_[0];

  // ensure the source term from subsurface -> surface is 0
  Teuchos::RCP<CompositeVector> source =
      S_next_->GetFieldData("overland_source_from_subsurface", richards_pk->name());
  source->PutScalar(0.);

  // evaluate the residual for overland flow
  Teuchos::RCP<TreeVector> overland_u_new = u_new->SubVector(overland_pk->name());
  Teuchos::RCP<TreeVector> overland_g = g->SubVector(overland_pk->name());
  overland_pk->fun(t_old, t_new, Teuchos::null, overland_u_new, overland_g);

  // Stuff the residual value into the source term.  This is the mass
  // imbalance on the surface, and is used in the subsurface PK.
  *source->ViewComponent("cell",false) = *overland_g->data()->ViewComponent("cell",false);

  // Evaluate the Richards residual, which will update the source with the true value.
  Teuchos::RCP<TreeVector> richards_u_new = u_new->SubVector(richards_pk->name());
  Teuchos::RCP<TreeVector> richards_g = g->SubVector(richards_pk->name());
  richards_pk->fun(t_old, t_new, Teuchos::null, richards_u_new, richards_g);

  // Update the overland residual with this new source term, res = res - Q
  overland_g->data()->ViewComponent("cell",false)->Update(-1.,
          *source->ViewComponent("cell",false), 1.0);
  source->Print(std::cout);

}

void MPCSurfaceSubsurfaceWaterCoupler::update_precon(double t,
        Teuchos::RCP<const TreeVector> up, double h) {
  Teuchos::RCP<TreeVector> up_nc = Teuchos::rcp_const_cast<TreeVector>(up);

  Teuchos::RCP<TreeVector> g = Teuchos::rcp(new TreeVector(*up));
  fun(S_inter_->time(), t, Teuchos::null, up_nc, g);
  StrongMPC::update_precon(t,up,h);
}


} // namespace
