/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Interface for the derived MPC for neumann coupling between surface and subsurface.

In this method, a flux BC is used on the subsurface.  The operator looks like:

            p      lamda    p_surf
         -------------------------
p       | A-K/dz |  K/dz |   0    | ( dp  )  =  res_mass_sub
         -------------------------
lambda  |   K/dz | -K/dz |   0    | ( dl  )  =  res_flux BC
         -------------------------
p_surf  |   K/dz |   0   |As-K/dz | ( dps )  =  res_mass_surf
         -------------------------

With this operator, for a precon we can solve the (p,lambda) block, then back
substitute dp and solve the p_surf block.

To do this, we need the flux BC, which is given by the residual for mass on the
surface.


------------------------------------------------------------------------- */

#ifndef PKS_MPC_SURFACE_SUBSURFACE_NEUMANN_COUPLER_HH_
#define PKS_MPC_SURFACE_SUBSURFACE_NEUMANN_COUPLER_HH_

#include "richards.hh"
#include "overland.hh"
#include "strong_mpc.hh"

namespace Amanzi {

class MPCSurfaceSubsurfaceNeumannCoupler : public StrongMPC {

 public:
  MPCSurfaceSubsurfaceNeumannCoupler(Teuchos::ParameterList& plist,
          const Teuchos::RCP<TreeVector>& soln) :
      PKDefaultBase(plist, soln),
      StrongMPC(plist, soln) {}

  // Virtual destructor
  virtual ~MPCSurfaceSubsurfaceNeumannCoupler() {}

  // initialization
  virtual void setup(const Teuchos::Ptr<State>& S);

  // evaluate the residual
  virtual void fun(double t_old, double t_new, Teuchos::RCP<TreeVector> u_old,
           Teuchos::RCP<TreeVector> u_new, Teuchos::RCP<TreeVector> g);

  // applies preconditioner to u and returns the result in Pu
  virtual void precon(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu);

 protected:
  Teuchos::RCP<PKBDFBase> pk_flow_;
  Teuchos::RCP<PKBDFBase> pk_ol_;

 private:
  // factory registration
  static RegisteredPKFactory<MPCSurfaceSubsurfaceNeumannCoupler> reg_;

};

} // namespace


#endif
