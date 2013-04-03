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

#ifndef PKS_MPC_SURFACE_SUBSURFACE_DIRICHLET_COUPLER_HH_
#define PKS_MPC_SURFACE_SUBSURFACE_DIRICHLET_COUPLER_HH_

#include "mpc_surface_subsurface_coupler.hh"

namespace Amanzi {

class MPCSurfaceSubsurfaceDirichletCoupler : public MPCSurfaceSubsurfaceCoupler {

 public:
  MPCSurfaceSubsurfaceDirichletCoupler(Teuchos::ParameterList& plist,
          const Teuchos::RCP<TreeVector>& soln) :
      PKDefaultBase(plist, soln),
      MPCSurfaceSubsurfaceCoupler(plist, soln) {}

  // applies preconditioner to u and returns the result in Pu
  virtual void precon(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu);

  // updates the preconditioner
  virtual void update_precon(double t, Teuchos::RCP<const TreeVector> up, double h);

 protected:
  // Apply the preconditioner matrix.
  virtual void PreconApply_(Teuchos::RCP<const TreeVector> u,
                            Teuchos::RCP<TreeVector> Pu);

  // Hackery hook for inheriting MPCs.
  virtual void PreconPostprocess_(Teuchos::RCP<const TreeVector> u,
          Teuchos::RCP<TreeVector> Pu) {};

  // Given updates to subsurface, calculate updates to surface cells.
  virtual void PreconUpdateSurfaceCells_(Teuchos::RCP<const TreeVector> u,
          Teuchos::RCP<TreeVector> Pu);

  // Given updates to surface cells, calculate updates to surface faces.
  virtual void PreconUpdateSurfaceFaces_(Teuchos::RCP<const TreeVector> u,
          Teuchos::RCP<TreeVector> Pu);

 private:
  // factory registration
  static RegisteredPKFactory<MPCSurfaceSubsurfaceDirichletCoupler> reg_;

};

} // namespace

#endif

