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

// -- Setup data.
void MPCSurfaceSubsurfaceFullCoupler::setup(const Teuchos::Ptr<State>& S) {
  MPCSurfaceSubsurfaceCoupler::setup(S);

  // Get the domain's preconditioner and replace it with a MatrixMFD_Surf.
  Teuchos::RCP<Operators::Matrix> precon = domain_pk_->preconditioner();
  Teuchos::RCP<Operators::MatrixMFD> mfd_precon =
    Teuchos::rcp_dynamic_cast<Operators::MatrixMFD>(precon);
  ASSERT(mfd_precon != Teuchos::null);

  preconditioner_ = Teuchos::rcp(new Operators::MatrixMFD_Surf(*mfd_precon, surf_mesh_));

  // Get the surface's preconditioner and replace it with a TPFA matrix.
  Teuchos::RCP<Operators::Matrix> surf_precon = surf_pk_->preconditioner();
  Teuchos::RCP<Operators::MatrixMFD> mfd_surf_precon =
    Teuchos::rcp_dynamic_cast<Operators::MatrixMFD>(surf_precon);
  ASSERT(mfd_surf_precon != Teuchos::null);

  Teuchos::RCP<Operators::MatrixMFD_TPFA> new_surf_precon =
    Teuchos::rcp(new Operators::MatrixMFD_TPFA(*mfd_surf_precon));

  // set the surface A in the MFD_Surf.
  preconditioner_->set_surface_A(new_surf_precon);

  // give the PCs back to the PKs
  domain_pk_->set_preconditioner(preconditioner_);
  surf_pk_->set_preconditioner(new_surf_precon);

}

// applies preconditioner to u and returns the result in Pu
void MPCSurfaceSubsurfaceFullCoupler::precon(Teuchos::RCP<const TreeVector> u,
        Teuchos::RCP<TreeVector> Pu) {
  ASSERT(0);
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

