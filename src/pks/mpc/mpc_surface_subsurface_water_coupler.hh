/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Interface for the derived MPC for water coupling between surface and subsurface.
------------------------------------------------------------------------- */

#ifndef PKS_MPC_SURFACE_SUBSURFACE_WATER_COUPLER_HH_
#define PKS_MPC_SURFACE_SUBSURFACE_WATER_COUPLER_HH_

#include "richards.hh"
#include "overland.hh"
#include "strong_mpc.hh"

namespace Amanzi {

class MPCSurfaceSubsurfaceWaterCoupler : public StrongMPC {

 public:
  MPCSurfaceSubsurfaceWaterCoupler(Teuchos::ParameterList& plist,
          const Teuchos::RCP<TreeVector>& soln) :
      PKDefaultBase(plist, soln),
      StrongMPC(plist, soln) {}

  // initialization
  virtual void setup(const Teuchos::Ptr<State>& S);

  // evaluate the residual
  virtual void fun(double t_old, double t_new, Teuchos::RCP<TreeVector> u_old,
           Teuchos::RCP<TreeVector> u_new, Teuchos::RCP<TreeVector> g);

  virtual void update_precon(double t, Teuchos::RCP<const TreeVector> up, double h);

 private:
  // factory registration
  static RegisteredPKFactory<MPCSurfaceSubsurfaceWaterCoupler> reg_;

};

} // namespace


#endif
