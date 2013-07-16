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
#include "MatrixMFD_Surf_ScaledConstraint.hh"
#include "MatrixMFD_TPFA_ScaledConstraint.hh"
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

  // cells to debug
  if (plist_.isParameter("debug cells")) {
    Teuchos::Array<int> surf_dc = plist_.get<Teuchos::Array<int> >("surface debug cells");
    for (Teuchos::Array<int>::const_iterator lcv=surf_dc.begin();
         lcv!=surf_dc.end(); ++lcv) {
      surf_dc_.push_back(*lcv);
    }
  } else {
    surf_dc_.push_back(plist_.get<int>("surface debug cell 0",0));
    surf_dc_.push_back(plist_.get<int>("surface debug cell 1",0));
  }
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

  // debugging cells
  if (plist_.get<bool>("debug all surface cells", false)) {
    dc_.clear();
    surf_dc_.clear();

    int ncells = surf_mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
    dc_.resize(ncells);
    surf_dc_.resize(ncells);
    for (int c=0; c!=ncells; ++c) {
      surf_dc_[c] = c;
      AmanziMesh::Entity_ID f = surf_mesh_->entity_get_parent(AmanziMesh::CELL, c);
      AmanziMesh::Entity_ID_List cells;
      domain_mesh_->face_get_cells(f, AmanziMesh::USED, &cells);
      ASSERT(cells.size() == 1);
      dc_[c] = cells[0];
    }
  }

  // Replace the subdomain preconditioners with the necessary versions.
  Teuchos::ParameterList surface_pc_plist = plist_.sublist("PKs")
      .sublist(surf_pk_name_).sublist("Diffusion PC");
  Teuchos::ParameterList subsurface_pc_plist = plist_.sublist("PKs")
      .sublist(domain_pk_name_).sublist("Diffusion PC");
  Teuchos::ParameterList surface_op_plist = plist_.sublist("PKs")
      .sublist(surf_pk_name_).sublist("Diffusion");
  Teuchos::ParameterList subsurface_op_plist = plist_.sublist("PKs")
      .sublist(domain_pk_name_).sublist("Diffusion");

  if (subsurface_op_plist.get<bool>("scaled constraint equation", false)) {
    mfd_preconditioner_ = Teuchos::rcp(new Operators::MatrixMFD_Surf_ScaledConstraint(
        subsurface_pc_plist, domain_mesh_, surf_mesh_));
  } else {
    mfd_preconditioner_ = Teuchos::rcp(new Operators::MatrixMFD_Surf(
        subsurface_pc_plist, domain_mesh_, surf_mesh_));
  }
  preconditioner_ = mfd_preconditioner_;

  if (surface_op_plist.get<bool>("scaled constraint equation", false)) {
    surf_preconditioner_ = Teuchos::rcp(new Operators::MatrixMFD_TPFA_ScaledConstraint(
        surface_pc_plist, surf_mesh_));
  } else {
    surf_preconditioner_ = Teuchos::rcp(new Operators::MatrixMFD_TPFA(
        surface_pc_plist, surf_mesh_));
  }

  // set the surface A in the MFD_Surf.
  mfd_preconditioner_->SetSurfaceOperator(surf_preconditioner_);

  // give the PCs back to the PKs
  domain_pk_->set_preconditioner(preconditioner_);
  surf_pk_->set_preconditioner(surf_preconditioner_);

}

bool MPCSurfaceSubsurfaceCoupler::modify_predictor_copy_surf_to_subsurf_(double h,
        Teuchos::RCP<TreeVector> up) {
  Teuchos::OSTab tab = getOSTab();
  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_HIGH, true))
    *out_ << "  Modifying predictor, copying surf -> sub-surf" << std::endl;

  Epetra_MultiVector& domain_u_f = *up->SubVector(domain_pk_name_)->data()
      ->ViewComponent("face",false);
  const Epetra_MultiVector& surf_u_c = *up->SubVector(surf_pk_name_)->data()
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
  Teuchos::OSTab tab = getOSTab();
  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_HIGH, true))
    *out_ << "  Modifying predictor, copying sub-surf -> surf" << std::endl;

  const Epetra_MultiVector& domain_u_f = *up->SubVector(domain_pk_name_)->data()
      ->ViewComponent("face",false);
  Epetra_MultiVector& surf_u_c = *up->SubVector(surf_pk_name_)->data()
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
  bool changed = domain_pk_->modify_predictor(h, up->SubVector(domain_pk_name_));

  // Ensure surface face values match subsurface values
  if (changed) modify_predictor_copy_subsurf_to_surf_(h, up);

  // Call the surface modify_predictor() to update surface faces.
  changed |= surf_pk_->modify_predictor(h, up->SubVector(surf_pk_name_));

  return changed;
}
} // namespace
