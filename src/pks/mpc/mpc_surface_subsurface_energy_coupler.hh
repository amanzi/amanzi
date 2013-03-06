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

#ifndef PKS_MPC_SURFACE_SUBSURFACE_ENERGY_COUPLER_HH_
#define PKS_MPC_SURFACE_SUBSURFACE_ENERGY_COUPLER_HH_

#include "mpc_surface_subsurface_coupler.hh"

namespace Amanzi {

namespace Operators { class MatrixMFD_Surf; }

class MPCSurfaceSubsurfaceEnergyCoupler : public MPCSurfaceSubsurfaceCoupler {

 public:
  MPCSurfaceSubsurfaceEnergyCoupler(Teuchos::ParameterList& plist,
          const Teuchos::RCP<TreeVector>& soln) :
      PKDefaultBase(plist, soln),
      MPCSurfaceSubsurfaceCoupler(plist, soln) {}

  // -- Setup data.
  virtual void setup(const Teuchos::Ptr<State>& S);

  // applies preconditioner to u and returns the result in Pu
  virtual void precon(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu);

  // updates the preconditioner
  virtual void update_precon(double t, Teuchos::RCP<const TreeVector> up, double h);

  virtual bool modify_predictor(double h, Teuchos::RCP<TreeVector> u);

 protected:
  Teuchos::RCP<Operators::MatrixMFD_Surf> preconditioner_;
  Teuchos::RCP<Operators::MatrixMFD_TPFA> surf_preconditioner_;

 private:
  // factory registration
  static RegisteredPKFactory<MPCSurfaceSubsurfaceEnergyCoupler> reg_;

};

} // namespace

#endif

