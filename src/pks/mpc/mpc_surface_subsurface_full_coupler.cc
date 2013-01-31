/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Interface for the derived MPC for water coupling between surface and subsurface.

In this method, a Dirichlet BC is used on the subsurface boundary for
the operator, but the preconditioner is for the full system with no
extra unknowns.  On the surface, the TPFA is used, resulting in a
subsurface-face-only Schur complement that captures all terms.

------------------------------------------------------------------------- */

#include "matrix_mfd_surf.hh"
#include "matrix_mfd_tpfa.hh"
#include "pk_physical_bdf_base.hh"

#include "mpc_surface_subsurface_full_coupler.hh"


namespace Amanzi {

#define DEBUG_FLAG 1

RegisteredPKFactory<MPCSurfaceSubsurfaceFullCoupler>
MPCSurfaceSubsurfaceFullCoupler::reg_("surface-subsurface full coupler");

// -- Setup data.
void MPCSurfaceSubsurfaceFullCoupler::setup(const Teuchos::Ptr<State>& S) {
  MPCSurfaceSubsurfaceCoupler::setup(S);

  // Get the domain's preconditioner and replace it with a MatrixMFD_Surf.
  Teuchos::RCP<Operators::Matrix> precon = domain_pk_->preconditioner();
  Teuchos::RCP<Operators::MatrixMFD> mfd_precon =
    Teuchos::rcp_dynamic_cast<Operators::MatrixMFD>(precon);
  ASSERT(mfd_precon != Teuchos::null);

  preconditioner_ = Teuchos::rcp(new Operators::MatrixMFD_Surf(*mfd_precon, surf_mesh_));

  // Get the surface's preconditioner and ensure it is TPFA.
  Teuchos::RCP<Operators::Matrix> surf_precon = surf_pk_->preconditioner();
  surf_preconditioner_ =
    Teuchos::rcp_dynamic_cast<Operators::MatrixMFD_TPFA>(surf_precon);
  ASSERT(surf_preconditioner_ != Teuchos::null);

  // set the surface A in the MFD_Surf.
  preconditioner_->set_surface_A(surf_preconditioner_);

  // give the PCs back to the PKs
  domain_pk_->set_preconditioner(preconditioner_);
  surf_pk_->set_preconditioner(surf_preconditioner_);

}

// applies preconditioner to u and returns the result in Pu
void MPCSurfaceSubsurfaceFullCoupler::precon(Teuchos::RCP<const TreeVector> u,
        Teuchos::RCP<TreeVector> Pu) {
  // Grab the surface cell residuals and stick them in the corresponding
  // subsurface face locations.
  Teuchos::RCP<const CompositeVector> domain_u = u->SubVector(domain_pk_name_)->data();
  Teuchos::RCP<CompositeVector> domain_u_new = Teuchos::rcp(new CompositeVector(*domain_u));
  domain_u_new->CreateData();
  *domain_u_new = *domain_u;

  Teuchos::RCP<const CompositeVector> surf_u = u->SubVector(surf_pk_name_)->data();
  Teuchos::RCP<CompositeVector> domain_Pu = Pu->SubVector(domain_pk_name_)->data();
  Teuchos::RCP<CompositeVector> surf_Pu = Pu->SubVector(surf_pk_name_)->data();


  const Epetra_MultiVector& surf_u_c = *u->SubVector(surf_pk_name_)
      ->data()->ViewComponent("cell",false);
  Epetra_MultiVector& domain_u_new_f = *domain_u_new->ViewComponent("face",false);

  int ncells_surf = surf_u_c.MyLength();
  for (int cs=0; cs!=ncells_surf; ++cs) {
    AmanziMesh::Entity_ID f = surf_mesh_->entity_get_parent(AmanziMesh::CELL, cs);
    domain_u_new_f[0][f] = surf_u_c[0][cs];
  }

  // Apply the combined preconditioner
  preconditioner_->ApplyInverse(*domain_u_new, domain_Pu.ptr());

#if DEBUG_FLAG
  Teuchos::OSTab tab = getOSTab();
  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_HIGH, true)) {
    *out_ << "Preconditioner application" << std::endl;
    *out_ << " SubSurface precon:" << std::endl;
    *out_ << "  p0: " << (*domain_u_new)("cell",0) << " " << (*domain_u_new)("face",3)
          << std::endl;
    *out_ << "  p1: " << (*domain_u_new)("cell",99) << " " << (*domain_u_new)("face",500)
          << std::endl;
    *out_ << "  PC*p0: " << (*domain_Pu)("cell",0) << " " << (*domain_Pu)("face",3)
          << std::endl;
    *out_ << "  PC*p1: " << (*domain_Pu)("cell",99) << " " << (*domain_Pu)("face",500)
          << std::endl;
  }
#endif

  // Correction applies to both the domain face and the surface cell.
  // Additionally, drift in the difference between lambda and p_surf
  // needs to be corrected.
  Epetra_MultiVector& domain_Pu_f = *domain_Pu->ViewComponent("face",false);
  Epetra_MultiVector& surf_Pu_c = *surf_Pu->ViewComponent("cell",false);
  Epetra_MultiVector& surf_Pu_f = *surf_Pu->ViewComponent("face",false);
  const Epetra_MultiVector& p_surf = *S_next_->GetFieldData("surface_pressure")->ViewComponent("cell",false);
  const Epetra_MultiVector& lambda = *S_next_->GetFieldData("pressure")->ViewComponent("face",false);

  for (int cs=0; cs!=ncells_surf; ++cs) {
    AmanziMesh::Entity_ID f = surf_mesh_->entity_get_parent(AmanziMesh::CELL, cs);
    surf_Pu_c[0][cs] = domain_Pu_f[0][f];
    domain_Pu_f[0][f] += lambda[0][f] - p_surf[0][cs];
  }

  // // Update surface with consistent faces
  Teuchos::RCP<CompositeVector> soln =
    S_next_->GetFieldData("surface_pressure",surf_pk_name_);
  Epetra_MultiVector& soln_c = *soln->ViewComponent("cell",false);
  Epetra_MultiVector& soln_f = *soln->ViewComponent("face",false);

  // -- calculate the new soln on cells
  soln_c.Update(-1., surf_Pu_c, 1.);
  surf_pk_->changed_solution();

  // === BEGIN CRAP CODE ===
  // -- enforce BCs
  S_next_->GetFieldEvaluator("ponded_depth")->HasFieldChanged(S_next_.ptr(), surf_pk_name_);
  Teuchos::RCP<const CompositeVector> pd = S_next_->GetFieldData("ponded_depth");
  pd->ScatterMasterToGhosted("cell");
  const Epetra_MultiVector& pd_c = *pd->ViewComponent("cell",true);

  AmanziMesh::Entity_ID_List cells;
  int nfaces_surf = surf_Pu_f.MyLength();
  for (int fs=0; fs!=nfaces_surf; ++fs) {
    surf_mesh_->face_get_cells(fs, AmanziMesh::USED, &cells);
    if (cells.size() == 1) {
      surf_Pu_f[0][fs] = soln_f[0][fs] - pd_c[0][cells[0]];
      std::cout << "Old, new, diff = " << soln_f[0][fs] << ", " << pd_c[0][cells[0]] << ", " << surf_Pu_f[0][fs] << std::endl;
    }
  }

  soln_c.Update(1., surf_Pu_c, 1.);
  // === END CRAP CODE ===


#if DEBUG_FLAG
  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_HIGH, true)) {
    *out_ << " Surface precon:" << std::endl;
    *out_ << "  p0: " << (*surf_u)("cell",0) << std::endl;
    *out_ << "  PC*p0: " << (*surf_Pu)("cell",0) << std::endl;
  }
#endif
}

// updates the preconditioner
void MPCSurfaceSubsurfaceFullCoupler::update_precon(double t,
        Teuchos::RCP<const TreeVector> up, double h) {
  MPCSurfaceSubsurfaceCoupler::update_precon(t, up, h);
  preconditioner_->AssembleGlobalMatrices();
  preconditioner_->ComputeSchurComplement(domain_pk_->bc_markers(),
          domain_pk_->bc_values());
  preconditioner_->UpdatePreconditioner();
}


} // namespace

