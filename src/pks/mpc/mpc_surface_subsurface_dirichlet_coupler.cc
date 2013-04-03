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

#include "matrix_mfd_surf.hh"
#include "matrix_mfd_tpfa.hh"
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
  // Apply the combined preconditioner to the subsurface residual
  PreconApply_(u,Pu);

  // Damp, kluge, hack, etc.
  PreconPostprocess_(u,Pu);

  // Update surface values.
  PreconUpdateSurfaceCells_(u,Pu);
  PreconUpdateSurfaceFaces_(u,Pu);

#if DEBUG_FLAG
  Teuchos::OSTab tab = getOSTab();
  Teuchos::RCP<const CompositeVector> surf_u = u->SubVector(surf_pk_name_)->data();
  Teuchos::RCP<const CompositeVector> domain_u = u->SubVector(domain_pk_name_)->data();
  Teuchos::RCP<CompositeVector> domain_Pu = Pu->SubVector(domain_pk_name_)->data();
  Teuchos::RCP<CompositeVector> surf_Pu = Pu->SubVector(surf_pk_name_)->data();

  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_HIGH, true)) {
    AmanziMesh::Entity_ID_List fnums,fnums0;
    std::vector<int> dirs;
    domain_mesh_->cell_get_faces_and_dirs(c0_, &fnums0, &dirs);
    domain_mesh_->cell_get_faces_and_dirs(c1_, &fnums, &dirs);

    *out_ << "Preconditioner application" << std::endl;
    *out_ << " SubSurface precon:" << std::endl;
    *out_ << "  p0: " << (*domain_u)("cell",c0_) << " "
        //          << (*domain_u)("face",fnums0[0])
          << (*domain_u)("face",fnums0[0]) << ", "
          << (*domain_u)("face",fnums0[1]) << ", "
          << (*domain_u)("face",fnums0[2]) << ", "
          << (*domain_u)("face",fnums0[3])
          << std::endl;
    *out_ << "  p1: " << (*domain_u)("cell",c1_) << " "
        //          << (*domain_u)("face",fnums[0])
          << (*domain_u)("face",fnums[0]) << ", "
          << (*domain_u)("face",fnums[1]) << ", "
          << (*domain_u)("face",fnums[2]) << ", "
          << (*domain_u)("face",fnums[3])
          << std::endl;
    *out_ << "  PC*p0: " << (*domain_Pu)("cell",c0_) << " "
        //          << (*domain_Pu)("face",fnums0[0])
          << (*domain_Pu)("face",fnums0[0]) << ", "
          << (*domain_Pu)("face",fnums0[1]) << ", "
          << (*domain_Pu)("face",fnums0[2]) << ", "
          << (*domain_Pu)("face",fnums0[3])
          << std::endl;
    *out_ << "  PC*p1: " << (*domain_Pu)("cell",c1_) << " "
        //          << (*domain_Pu)("face",fnums[0])
          << (*domain_Pu)("face",fnums[0]) << ", "
          << (*domain_Pu)("face",fnums[1]) << ", "
          << (*domain_Pu)("face",fnums[2]) << ", "
          << (*domain_Pu)("face",fnums[3])
          << std::endl;
  }

  if (surf_c0_ < surf_u->size("cell",false) && surf_c1_ < surf_u->size("cell",false)) {
    //  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_HIGH, true)) {
     AmanziMesh::Entity_ID_List fnums,fnums0;
     std::vector<int> dirs;
     surf_mesh_->cell_get_faces_and_dirs(surf_c0_, &fnums0, &dirs);
     surf_mesh_->cell_get_faces_and_dirs(surf_c1_, &fnums, &dirs);

    *out_ << " Surface precon:" << std::endl;
    *out_ << "  u0: " << (*surf_u)("cell",surf_c0_) << ", "
          << (*surf_u)("face",fnums0[0]) << std::endl;
    *out_ << "  u1: " << (*surf_u)("cell",surf_c1_) << ", "
          << (*surf_u)("face",fnums[0]) << std::endl;
    *out_ << "  PC*u0: " << (*surf_Pu)("cell",surf_c0_) << ", "
          << (*surf_Pu)("face",fnums0[0]) << std::endl;
    *out_ << "  PC*u1: " << (*surf_Pu)("cell",surf_c1_) << ", "
          << (*surf_Pu)("face",fnums[0]) << std::endl;

  }
#endif
}


// -------------------------------------------------------------
// Apply preconditioner to u and returns the result in Pu
// -------------------------------------------------------------
void MPCSurfaceSubsurfaceDirichletCoupler::PreconApply_(
    Teuchos::RCP<const TreeVector> u,
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
  //  mfd_preconditioner_->ApplyInverse(*domain_u_new, domain_Pu.ptr());
  Teuchos::RCP<TreeVector> domain_u_TV = Teuchos::rcp(new TreeVector("domain_u_TV"));
  domain_u_TV->set_data(domain_u_new);
  Teuchos::RCP<TreeVector> domain_Pu_TV = Teuchos::rcp(new TreeVector("domain_Pu_TV"));
  domain_Pu_TV->set_data(domain_Pu);
  domain_pk_->precon(domain_u_TV, domain_Pu_TV);
}



// ------------------------------------------------------------------------------
// Correction applies to both the domain face and the surface cell.
// Additionally, drift in the difference between lambda and p_surf needs to be
// corrected.
// ------------------------------------------------------------------------------
void MPCSurfaceSubsurfaceDirichletCoupler::PreconUpdateSurfaceCells_(
    Teuchos::RCP<const TreeVector> u,
    Teuchos::RCP<TreeVector> Pu) {


  Teuchos::RCP<CompositeVector> domain_Pu = Pu->SubVector(domain_pk_name_)->data();
  Teuchos::RCP<CompositeVector> surf_Pu = Pu->SubVector(surf_pk_name_)->data();
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
  Teuchos::RCP<const CompositeVector> surf_u = u->SubVector(surf_pk_name_)->data();
  Teuchos::RCP<CompositeVector> surf_Pu = Pu->SubVector(surf_pk_name_)->data();
  surf_preconditioner_->UpdateConsistentFaceCorrection(*surf_u, surf_Pu.ptr());
}


// updates the preconditioner
void MPCSurfaceSubsurfaceDirichletCoupler::update_precon(double t,
        Teuchos::RCP<const TreeVector> up, double h) {
  MPCSurfaceSubsurfaceCoupler::update_precon(t, up, h);
  mfd_preconditioner_->AssembleGlobalMatrices();
  mfd_preconditioner_->ComputeSchurComplement(domain_pk_->bc_markers(),
          domain_pk_->bc_values());
  mfd_preconditioner_->UpdatePreconditioner();
}

} // namespace

