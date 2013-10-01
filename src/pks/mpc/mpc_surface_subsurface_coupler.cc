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

  // cells to debug, subsurface
  if (plist_.isParameter("debug cells")) {
    dc_.clear();
    Teuchos::Array<int> dc = plist_.get<Teuchos::Array<int> >("debug cells");
    for (Teuchos::Array<int>::const_iterator c=dc.begin();
         c!=dc.end(); ++c) dc_.push_back(*c);

    // Enable a vo for each cell, allows parallel printing of debug cells.
    if (plist_.isParameter("debug cell ranks")) {
      Teuchos::Array<int> dc_ranks = plist_.get<Teuchos::Array<int> >("debug cell ranks");
      if (dc.size() != dc_ranks.size()) {
        Errors::Message message("Debug cell and debug cell ranks must be equal length.");
        Exceptions::amanzi_throw(message);
      }
      for (Teuchos::Array<int>::const_iterator dcr=dc_ranks.begin();
           dcr!=dc_ranks.end(); ++dcr) {
        // make a verbose object for each case
        Teuchos::ParameterList vo_plist;
        vo_plist.sublist("VerboseObject");
        vo_plist.sublist("VerboseObject")
            = plist_.sublist("VerboseObject");
        vo_plist.sublist("VerboseObject").set("write on rank", *dcr);

        dcvo_.push_back(Teuchos::rcp(new VerboseObject(domain_mesh_->get_comm(), name_,vo_plist)));

      }
    } else {
      // Simply use the pk's vo
      dcvo_.resize(dc_.size(), vo_);
    }
  }

  // cells to debug, subsurface
  if (plist_.isParameter("surface debug cells")) {
    surf_dc_.clear();
    Teuchos::Array<int> surf_dc = plist_.get<Teuchos::Array<int> >("surface debug cells");
    for (Teuchos::Array<int>::const_iterator c=surf_dc.begin();
         c!=surf_dc.end(); ++c) surf_dc_.push_back(*c);

    // Enable a vo for each cell, allows parallel printing of debug cells.
    if (plist_.isParameter("surface debug cell ranks")) {
      Teuchos::Array<int> surf_dc_ranks = plist_.get<Teuchos::Array<int> >("surface debug cell ranks");
      if (surf_dc.size() != surf_dc_ranks.size()) {
        Errors::Message message("Debug cell and debug cell ranks must be equal length.");
        Exceptions::amanzi_throw(message);
      }
      for (Teuchos::Array<int>::const_iterator surf_dcr=surf_dc_ranks.begin();
           surf_dcr!=surf_dc_ranks.end(); ++surf_dcr) {
        // make a verbose object for each case
        Teuchos::ParameterList vo_plist;
        vo_plist.sublist("VerboseObject");
        vo_plist.sublist("VerboseObject")
            = plist_.sublist("VerboseObject");
        vo_plist.sublist("VerboseObject").set("write on rank", *surf_dcr);

        surf_dcvo_.push_back(Teuchos::rcp(new VerboseObject(surf_mesh_->get_comm(), name_,vo_plist)));

      }
    } else {
      // Simply use the pk's vo
      surf_dcvo_.resize(surf_dc_.size(), vo_);
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

  // CLEAN UP THIS CRUFT!  surf precon should ALWAYS be ScaledConstraint... we
  // should impose this from above, not backhandedly
  surf_preconditioner_ = Teuchos::rcp(new Operators::MatrixMFD_TPFA_ScaledConstraint(
      surface_pc_plist, surf_mesh_));


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

  Epetra_MultiVector& domain_u_f = *up->SubVector(domain_pk_name_)->Data()
      ->ViewComponent("face",false);
  const Epetra_MultiVector& surf_u_c = *up->SubVector(surf_pk_name_)->Data()
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

  const Epetra_MultiVector& domain_u_f = *up->SubVector(domain_pk_name_)->Data()
      ->ViewComponent("face",false);
  Epetra_MultiVector& surf_u_c = *up->SubVector(surf_pk_name_)->Data()
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
