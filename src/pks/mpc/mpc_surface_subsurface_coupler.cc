/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

This class should never be instantiated -- it just provides a base class for
multiple coupler types.
------------------------------------------------------------------------- */

#include "mpc_surface_subsurface_coupler.hh"

namespace Amanzi {

MPCSurfaceSubsurfaceCoupler::MPCSurfaceSubsurfaceCoupler(Teuchos::ParameterList& plist,
        const Teuchos::RCP<TreeVector>& soln) :
    PKDefaultBase(plist, soln),
    StrongMPC(plist,soln) {

  domain_mesh_key_ = plist.get<std::string>("subsurface mesh key");
  surf_mesh_key_ = plist.get<std::string>("surface mesh key");

  surf_pk_name_ = plist.get<std::string>("surface PK name");
  domain_pk_name_ = plist.get<std::string>("subsurface PK name");
};

void MPCSurfaceSubsurfaceCoupler::setup(const Teuchos::Ptr<State>& S) {
  surf_mesh_ = S->GetMesh(surf_mesh_key_);
  domain_mesh_ = S->GetMesh(domain_mesh_key_);

  StrongMPC::setup(S);

  if (sub_pks_[0]->name() == domain_pk_name_) {
    domain_pk_ = sub_pks_[0];
    ASSERT(sub_pks_[1]->name() == surf_pk_name_);
    surf_pk_ = sub_pks_[1];
  } else if (sub_pks_[1]->name() == domain_pk_name_) {
    domain_pk_ = sub_pks_[1];
    ASSERT(sub_pks_[0]->name() == surf_pk_name_);
    surf_pk_ = sub_pks_[0];
  } else {
    ASSERT(0);
  }

}




} // namespace
