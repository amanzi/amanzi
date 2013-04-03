/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

This class should never be instantiated -- it just provides a base class for
multiple coupler types.
------------------------------------------------------------------------- */

#include "pk_physical_bdf_base.hh"

#include "matrix_mfd_surf.hh"
#include "matrix_mfd_tpfa.hh"
#include "mpc_surface_subsurface_coupler.hh"

namespace Amanzi {

MPCSurfaceSubsurfaceCoupler::MPCSurfaceSubsurfaceCoupler(Teuchos::ParameterList& plist,
        const Teuchos::RCP<TreeVector>& soln) :
    PKDefaultBase(plist, soln),
    StrongMPC(plist,soln) {

  domain_mesh_key_ = plist.get<std::string>("subsurface mesh key","domain");
  surf_mesh_key_ = plist.get<std::string>("surface mesh key","surface");

  surf_pk_name_ = plist.get<std::string>("surface PK name");
  domain_pk_name_ = plist.get<std::string>("subsurface PK name");

  surf_c0_ = plist_.get<int>("surface debug cell 0", 0);
  surf_c1_ = plist_.get<int>("surface debug cell 1", 1);

};

void MPCSurfaceSubsurfaceCoupler::setup(const Teuchos::Ptr<State>& S) {
  StrongMPC::setup(S);

  // get the mesh info
  surf_mesh_ = S->GetMesh(surf_mesh_key_);
  domain_mesh_ = S->GetMesh(domain_mesh_key_);

  // get the pk info
  if (sub_pks_[0]->name() == domain_pk_name_) {
    domain_pk_ = Teuchos::rcp_dynamic_cast<PKPhysicalBDFBase>(sub_pks_[0]);
    ASSERT(domain_pk_ != Teuchos::null);
    ASSERT(sub_pks_[1]->name() == surf_pk_name_);
    surf_pk_ = Teuchos::rcp_dynamic_cast<PKPhysicalBDFBase>(sub_pks_[1]);
    ASSERT(surf_pk_ != Teuchos::null);
  } else if (sub_pks_[1]->name() == domain_pk_name_) {
    domain_pk_ = Teuchos::rcp_dynamic_cast<PKPhysicalBDFBase>(sub_pks_[1]);
    ASSERT(domain_pk_ != Teuchos::null);
    ASSERT(sub_pks_[0]->name() == surf_pk_name_);
    surf_pk_ = Teuchos::rcp_dynamic_cast<PKPhysicalBDFBase>(sub_pks_[0]);
    ASSERT(surf_pk_ != Teuchos::null);
  } else {
    ASSERT(0);
  }

  // Get the domain's preconditioner and replace it with a MatrixMFD_Surf.
  Teuchos::RCP<Operators::Matrix> precon = domain_pk_->preconditioner();
  Teuchos::RCP<Operators::MatrixMFD> mfd_precon =
    Teuchos::rcp_dynamic_cast<Operators::MatrixMFD>(precon);
  ASSERT(mfd_precon != Teuchos::null);

  mfd_preconditioner_ =
      Teuchos::rcp(new Operators::MatrixMFD_Surf(*mfd_precon, surf_mesh_));
  preconditioner_ = mfd_preconditioner_;

  // Get the surface's preconditioner and ensure it is TPFA.
  Teuchos::RCP<Operators::Matrix> surf_precon = surf_pk_->preconditioner();
  surf_preconditioner_ =
    Teuchos::rcp_dynamic_cast<Operators::MatrixMFD_TPFA>(surf_precon);
  ASSERT(surf_preconditioner_ != Teuchos::null);

  // set the surface A in the MFD_Surf.
  mfd_preconditioner_->SetSurfaceOperator(surf_preconditioner_);

  // give the PCs back to the PKs
  domain_pk_->set_preconditioner(preconditioner_);
  surf_pk_->set_preconditioner(surf_preconditioner_);

}


bool MPCSurfaceSubsurfaceCoupler::modify_predictor(double h,
        Teuchos::RCP<TreeVector> up) {
  // Do any modifications on the sub PKs
  bool changed = StrongMPC::modify_predictor(h, up);

  // ensure the two match
  if (changed) {
    Epetra_MultiVector& domain_u_f = *up->SubVector(domain_pk_name_)->data()
        ->ViewComponent("face",false);
    const Epetra_MultiVector& surf_u_c = *up->SubVector(surf_pk_name_)->data()
        ->ViewComponent("cell",false);

    int ncells = surf_u_c.MyLength();
    for (int c=0; c!=ncells; ++c) {
      int f = surf_mesh_->entity_get_parent(AmanziMesh::CELL, c);
      domain_u_f[0][f] = surf_u_c[0][c];
    }
  }

  return changed;
}




} // namespace
