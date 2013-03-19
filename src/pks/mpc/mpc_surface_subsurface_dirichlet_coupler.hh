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

namespace Operators {
  class MatrixMFD_Surf;
  class MatrixMFD_TPFA;
}

class MPCSurfaceSubsurfaceDirichletCoupler : public MPCSurfaceSubsurfaceCoupler {

 public:
  MPCSurfaceSubsurfaceDirichletCoupler(Teuchos::ParameterList& plist,
          const Teuchos::RCP<TreeVector>& soln) :
      PKDefaultBase(plist, soln),
      MPCSurfaceSubsurfaceCoupler(plist, soln) {
    surf_c0_ = plist_.get<int>("surface debug cell 0", 0);
    surf_c1_ = plist_.get<int>("surface debug cell 1", 1);
  }

  // -- Setup data.
  virtual void setup(const Teuchos::Ptr<State>& S);

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

 protected:
  Teuchos::RCP<Operators::MatrixMFD_Surf> mfd_preconditioner_;
  Teuchos::RCP<Operators::MatrixMFD_TPFA> surf_preconditioner_;

  int surf_c0_;
  int surf_c1_;

 private:
  // factory registration
  static RegisteredPKFactory<MPCSurfaceSubsurfaceDirichletCoupler> reg_;

};

} // namespace

#endif

