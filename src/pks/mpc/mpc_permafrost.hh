/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

MPC for the Coupled Permafrost model.  This MPC sits at the top of the
subtree:

                    MPCPermafrost
                     /          \
                    /            \
                   /              \
         surf/subsurf            surf/subsurf
           water                   energy
         /      \                  /      \
        /        \                /        \
    flow/        flow/         energy/     energy/
  permafrost  icy_overland    threephase    surface_ice

------------------------------------------------------------------------- */

#ifndef MPC_PERMAFROST_HH_
#define MPC_PERMAFROST_HH_

#include "strong_mpc.hh"

namespace Amanzi {

namespace Operators { class MatrixCoupledMFDSurf; }

class MPCPermafrost : public StrongMPC {

 public:
  MPCPermafrost(Teuchos::ParameterList& plist,
                             const Teuchos::RCP<TreeVector>& soln) :
      PKDefaultBase(plist, soln),
      StrongMPC(plist, soln)
  {}


  virtual void setup(const Teuchos::Ptr<State>& S);

  // applies preconditioner to u and returns the result in Pu
  virtual void precon(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu);

  // updates the preconditioner
  virtual void update_precon(double t, Teuchos::RCP<const TreeVector> up, double h);


 protected:
  Key dA_dy2_key_;
  Key dB_dy1_key_;
  Key A_key_;
  Key y2_key_;
  Key B_key_;
  Key y1_key_;

  bool decoupled_;

  Teuchos::RCP<MPCSurfaceSubsurfaceFluxCoupler> coupled_flow_pk_;
  Teuchos::RCP<MPCSurfaceSubsurfaceFluxCoupler> coupled_energy_pk_;

  Teuchos::RCP<Operators::MatrixCoupledMFDSurf> mfd_surf_preconditioner_;

 private:
  // factory registration
  static RegisteredPKFactory<MPCPermafrost> reg_;

};

} // namespace

#endif
