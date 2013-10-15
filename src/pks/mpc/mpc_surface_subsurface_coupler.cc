/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

This class should never be instantiated -- it just provides a base class for
multiple coupler types.
------------------------------------------------------------------------- */

#include "pk_physical_bdf_base.hh"

#include "MatrixMFD_Surf.hh"
#include "MatrixMFD_TPFA.hh"
#include "mpc_surface_subsurface_coupler.hh"

namespace Amanzi {

MPCSurfaceSubsurfaceCoupler::MPCSurfaceSubsurfaceCoupler(Teuchos::ParameterList& plist,
        Teuchos::ParameterList& FElist,
        const Teuchos::RCP<TreeVector>& soln) :
    PKDefaultBase(plist, FElist, soln),
    StrongMPC<PKPhysicalBDFBase>(plist, FElist, soln) {

  domain_mesh_key_ = plist.get<std::string>("subsurface mesh key","domain");
  surf_mesh_key_ = plist.get<std::string>("surface mesh key","surface");

  surf_pk_name_ = plist.get<std::string>("surface PK name");
  domain_pk_name_ = plist.get<std::string>("subsurface PK name");

  // ensure the subsurface preconditioner is a Surf
  plist.sublist("PKs").sublist(domain_pk_name_)
      .sublist("Diffusion PC").set("coupled to surface", true);
  // ensure the surface preconditioner is TPFA
  plist.sublist("PKs").sublist(surf_pk_name_)
      .sublist("Diffusion PC").set("TPFA", true);
};

void MPCSurfaceSubsurfaceCoupler::setup(const Teuchos::Ptr<State>& S) {
  StrongMPC<PKPhysicalBDFBase>::setup(S);

  // get the mesh info
  surf_mesh_ = S->GetMesh(surf_mesh_key_);
  domain_mesh_ = S->GetMesh(domain_mesh_key_);

  // get the pk info
  if (sub_pks_[0]->name() == domain_pk_name_) {
    domain_pk_ = sub_pks_[0];
    domain_pk_index_ = 0;
    ASSERT(sub_pks_[1]->name() == surf_pk_name_);
    surf_pk_ = sub_pks_[1];
    surf_pk_index_ = 1;
  } else if (sub_pks_[1]->name() == domain_pk_name_) {
    domain_pk_ = sub_pks_[1];
    domain_pk_index_ = 1;
    ASSERT(sub_pks_[0]->name() == surf_pk_name_);
    surf_pk_ = sub_pks_[0];
    surf_pk_index_ = 0;
  } else {
    ASSERT(0);
  }

  // get the preconditioners for our own use
  mfd_preconditioner_ = Teuchos::rcp_dynamic_cast<Operators::MatrixMFD_Surf>(domain_pk_->preconditioner());
  if (mfd_preconditioner_ == Teuchos::null) {
    Errors::Message msg("MPCSurfaceSubsurfaceCoupler: subsurface operator must be MatrixMFD_Surf.");
    Exceptions::amanzi_throw(msg);
  }

  surf_preconditioner_ = Teuchos::rcp_dynamic_cast<Operators::MatrixMFD_TPFA>(surf_pk_->preconditioner());
  if (surf_preconditioner_ == Teuchos::null) {
    Errors::Message msg("MPCSurfaceSubsurfaceCoupler: surface operator must be MatrixMFD_TPFA.");
    Exceptions::amanzi_throw(msg);
  }

  mfd_preconditioner_->SetSurfaceOperator(surf_preconditioner_);

  // get the PKs debuggers for our own use
  domain_db_ = domain_pk_->debugger();
  surf_db_ = surf_pk_->debugger();
}

bool MPCSurfaceSubsurfaceCoupler::modify_predictor_copy_surf_to_subsurf_(double h,
        Teuchos::RCP<TreeVector> up) {
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_HIGH))
    *vo_->os() << "  Modifying predictor, copying surf -> sub-surf" << std::endl;

  Epetra_MultiVector& domain_u_f = *up->SubVector(domain_pk_index_)->Data()
      ->ViewComponent("face",false);
  const Epetra_MultiVector& surf_u_c = *up->SubVector(surf_pk_index_)->Data()
      ->ViewComponent("cell",false);

  int ncells = surf_u_c.MyLength();
  for (int c=0; c!=ncells; ++c) {
    int f = surf_mesh_->entity_get_parent(AmanziMesh::CELL, c);
    domain_u_f[0][f] = surf_u_c[0][c];
  }

  return true;
}


bool MPCSurfaceSubsurfaceCoupler::modify_predictor_copy_subsurf_to_surf_(double h,
        Teuchos::RCP<TreeVector> up) {
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_HIGH))
    *vo_->os() << "  Modifying predictor, copying sub-surf -> surf" << std::endl;

  const Epetra_MultiVector& domain_u_f = *up->SubVector(domain_pk_index_)->Data()
      ->ViewComponent("face",false);
  Epetra_MultiVector& surf_u_c = *up->SubVector(surf_pk_index_)->Data()
      ->ViewComponent("cell",false);

  int ncells = surf_u_c.MyLength();
  for (int c=0; c!=ncells; ++c) {
    int f = surf_mesh_->entity_get_parent(AmanziMesh::CELL, c);
    surf_u_c[0][c] = domain_u_f[0][f];
  }

  return true;
}

bool MPCSurfaceSubsurfaceCoupler::modify_predictor(double h,
        Teuchos::RCP<TreeVector> up) {

  // In most cases, we first deal with the subsurface, then the surface.
  bool changed = domain_pk_->modify_predictor(h, up->SubVector(domain_pk_index_));

  // Ensure surface face values match subsurface values
  if (changed) modify_predictor_copy_subsurf_to_surf_(h, up);

  // Call the surface modify_predictor() to update surface faces.
  changed |= surf_pk_->modify_predictor(h, up->SubVector(surf_pk_index_));

  return changed;
}
} // namespace
