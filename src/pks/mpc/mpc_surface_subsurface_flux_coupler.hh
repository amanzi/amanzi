/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Interface for the derived MPC for water coupling between surface and subsurface.

In this method, a Neumann BC is used on the subsurface boundary for
the operator, but the preconditioner is for the flux system with no
extra unknowns.  On the surface, the TPFA is used, resulting in a
subsurface-face-only Schur complement that captures all terms.

------------------------------------------------------------------------- */

#ifndef PKS_MPC_SURFACE_SUBSURFACE_FLUX_COUPLER_HH_
#define PKS_MPC_SURFACE_SUBSURFACE_FLUX_COUPLER_HH_

#include "mpc_surface_subsurface_coupler.hh"

namespace Amanzi {

class MPCSurfaceSubsurfaceFluxCoupler : public MPCSurfaceSubsurfaceCoupler {

 public:
  MPCSurfaceSubsurfaceFluxCoupler(Teuchos::ParameterList& plist,
          const Teuchos::RCP<TreeVector>& soln) :
      PKDefaultBase(plist, soln),
      MPCSurfaceSubsurfaceCoupler(plist, soln) {
    modify_predictor_flux_bc_ =
      plist_.get<bool>("modify predictor for flux BCs", false);
  }

  // -- Setup data.
  virtual void setup(const Teuchos::Ptr<State>& S);

  // evaluate the flux
  virtual void fun(double t_old, double t_new, Teuchos::RCP<TreeVector> u_old,
           Teuchos::RCP<TreeVector> u_new, Teuchos::RCP<TreeVector> g);

  // applies preconditioner to u and returns the result in Pu
  virtual void precon(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu);

  // updates the preconditioner
  virtual void update_precon(double t, Teuchos::RCP<const TreeVector> up, double h);

  virtual bool modify_predictor(double h, Teuchos::RCP<TreeVector> up);

 protected:
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
  Key flux_key_;
  bool modify_predictor_flux_bc_;

 private:
  // factory registration
  static RegisteredPKFactory<MPCSurfaceSubsurfaceFluxCoupler> reg_;

  friend class MPCPermafrost;
};

} // namespace

#endif

