/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Interface for the derived MPC for water coupling between surface and subsurface.

In this method, a Dirichlet BC is used on the subsurface boundary for
the operator, but the preconditioner is for the dirichlet system with no
extra unknowns.  On the surface, the TPFA is used, resulting in a
subsurface-face-only Schur complement that captures all terms.

------------------------------------------------------------------------- */
#include "EpetraExt_RowMatrixOut.h"

#include "MatrixMFD_Surf.hh"
#include "MatrixMFD_TPFA.hh"
#include "pk_physical_bdf_base.hh"

#include "mpc_surface_subsurface_dirichlet_coupler.hh"


namespace Amanzi {

#define DEBUG_FLAG 1

RegisteredPKFactory<MPCSurfaceSubsurfaceDirichletCoupler>
MPCSurfaceSubsurfaceDirichletCoupler::reg_("surface-subsurface Dirichlet coupler");

// -------------------------------------------------------------
// Apply preconditioner to u and returns the result in Pu
// -------------------------------------------------------------
void MPCSurfaceSubsurfaceDirichletCoupler::precon(Teuchos::RCP<const TreeVector> u,
        Teuchos::RCP<TreeVector> Pu) {
  Teuchos::OSTab tab = vo_->getOSTab();
  // Apply the combined preconditioner to the subsurface residual
  PreconApply_(u,Pu);

  // Update surface values.
  PreconUpdateSurfaceCells_(Pu);
  PreconUpdateSurfaceFaces_(u,Pu);

  // write residuals
  if (vo_->os_OK(Teuchos::VERB_HIGH)) *vo_->os() << "Preconditioner Application:" << std::endl;
  std::vector<std::string> vnames;
  vnames.push_back("  r_sub"); vnames.push_back("  PC*r_sub"); 
  std::vector< Teuchos::Ptr<const CompositeVector> > vecs;
  vecs.push_back(u->SubVector(0)->Data().ptr()); 
  vecs.push_back(Pu->SubVector(0)->Data().ptr()); 
  domain_db_->WriteVectors(vnames, vecs, true);

  vnames[0] = "  r_surf"; vnames[1] = "  PC*r_surf";
  vecs[0] = u->SubVector(1)->Data().ptr();
  vecs[1] = Pu->SubVector(1)->Data().ptr();
  surf_db_->WriteVectors(vnames, vecs, true);
}


// -------------------------------------------------------------
// Apply preconditioner to u and returns the result in Pu
// -------------------------------------------------------------
void MPCSurfaceSubsurfaceDirichletCoupler::PreconApply_(
    Teuchos::RCP<const TreeVector> u,
    Teuchos::RCP<TreeVector> Pu) {

  // Grab the surface cell residuals and stick them in the corresponding
  // subsurface face locations.
  Teuchos::RCP<const CompositeVector> domain_u = u->SubVector(domain_pk_index_)->Data();
  Teuchos::RCP<CompositeVector> domain_u_new = Teuchos::rcp(new CompositeVector(*domain_u));
  *domain_u_new = *domain_u;

  Teuchos::RCP<const CompositeVector> surf_u = u->SubVector(surf_pk_index_)->Data();
  Teuchos::RCP<CompositeVector> domain_Pu = Pu->SubVector(domain_pk_index_)->Data();
  Teuchos::RCP<CompositeVector> surf_Pu = Pu->SubVector(surf_pk_index_)->Data();


  const Epetra_MultiVector& surf_u_c = *u->SubVector(surf_pk_index_)
      ->Data()->ViewComponent("cell",false);
  Epetra_MultiVector& domain_u_new_f = *domain_u_new->ViewComponent("face",false);

  int ncells_surf = surf_u_c.MyLength();
  for (int cs=0; cs!=ncells_surf; ++cs) {
    AmanziMesh::Entity_ID f = surf_mesh_->entity_get_parent(AmanziMesh::CELL, cs);
    domain_u_new_f[0][f] = surf_u_c[0][cs];
  }

  // Apply the combined preconditioner
  //  mfd_preconditioner_->ApplyInverse(*domain_u_new, domain_Pu.ptr());
  Teuchos::RCP<TreeVector> domain_u_TV = Teuchos::rcp(new TreeVector());
  domain_u_TV->SetData(domain_u_new);
  Teuchos::RCP<TreeVector> domain_Pu_TV = Teuchos::rcp(new TreeVector());
  domain_Pu_TV->SetData(domain_Pu);
  domain_pk_->precon(domain_u_TV, domain_Pu_TV);
}


bool MPCSurfaceSubsurfaceDirichletCoupler::modify_correction(double h,
        Teuchos::RCP<const TreeVector> res,
        Teuchos::RCP<const TreeVector> u,
        Teuchos::RCP<TreeVector> du) {
  bool modified = PreconPostprocess_(res, du);
  if (modified) {
    PreconUpdateSurfaceCells_(du);
    PreconUpdateSurfaceFaces_(res,du);
  }
  return modified;
}



// ------------------------------------------------------------------------------
// Correction applies to both the domain face and the surface cell.
// Additionally, drift in the difference between lambda and p_surf needs to be
// corrected.
// ------------------------------------------------------------------------------
void MPCSurfaceSubsurfaceDirichletCoupler::PreconUpdateSurfaceCells_(
    Teuchos::RCP<TreeVector> Pu) {


  Teuchos::RCP<CompositeVector> domain_Pu = Pu->SubVector(domain_pk_index_)->Data();
  Teuchos::RCP<CompositeVector> surf_Pu = Pu->SubVector(surf_pk_index_)->Data();
  Epetra_MultiVector& domain_Pu_f = *domain_Pu->ViewComponent("face",false);
  Epetra_MultiVector& surf_Pu_c = *surf_Pu->ViewComponent("cell",false);

  int ncells_surf = surf_mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  const Epetra_MultiVector& domain_p_f = *S_next_->GetFieldData("pressure")
      ->ViewComponent("face",false);
  const Epetra_MultiVector& surf_p_c = *S_next_->GetFieldData("surface_pressure")
      ->ViewComponent("cell",false);

  for (int cs=0; cs!=ncells_surf; ++cs) {
    AmanziMesh::Entity_ID f = surf_mesh_->entity_get_parent(AmanziMesh::CELL, cs);
    surf_Pu_c[0][cs] = domain_Pu_f[0][f];
    domain_Pu_f[0][f] += domain_p_f[0][f] - surf_p_c[0][cs];
  }
}


// ------------------------------------------------------------------------------
// Calculate the surface face update.
// ------------------------------------------------------------------------------
void MPCSurfaceSubsurfaceDirichletCoupler::PreconUpdateSurfaceFaces_(
    Teuchos::RCP<const TreeVector> u,
    Teuchos::RCP<TreeVector> Pu) {
  // update delta faces
  Teuchos::RCP<const CompositeVector> surf_u = u->SubVector(surf_pk_index_)->Data();
  Teuchos::RCP<CompositeVector> surf_Pu = Pu->SubVector(surf_pk_index_)->Data();
  surf_preconditioner_->UpdateConsistentFaceCorrection(*surf_u, surf_Pu.ptr());
}


// updates the preconditioner
void MPCSurfaceSubsurfaceDirichletCoupler::update_precon(double t,
        Teuchos::RCP<const TreeVector> up, double h) {
  MPCSurfaceSubsurfaceCoupler::update_precon(t, up, h);

  if (assemble_preconditioner_) {
    mfd_preconditioner_->AssembleGlobalMatrices();
    mfd_preconditioner_->ComputeSchurComplement(domain_pk_->bc_markers(),
            domain_pk_->bc_values());
    mfd_preconditioner_->UpdatePreconditioner();
  }
}

} // namespace

