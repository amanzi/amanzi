/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Interface for the derived MPC for water coupling between surface and subsurface.

In this method, a Dirichlet BC is used on the subsurface.  The operator looks like:

            p      lamda    p_surf
         -------------------------
p       | A-K/dz |  K/dz |   0    | ( dp  )  =  res_mass_sub
         -------------------------
lambda  |     0  |    1  |   -1   | ( dl  )  =  res_Dirichlet BC
         -------------------------
p_surf  |   K/dz |   0   |As-K/dz | ( dps )  =  res_mass_surf
         -------------------------

This operator is pretty tightly coupled, so we ignore the (p_surf,p) entry
(row 3, col 0 block).  With this zeroed out, we can solve the (3,3) block,
back substitute into res_Dirichlet BC, adding 1*dps, and solve the (p,lambda)
system.

To do this, we need the Dirichlet BC, which is given by the value of pressure
on the surface.

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

 private:
  // factory registration
  static RegisteredPKFactory<MPCSurfaceSubsurfaceDirichletCoupler> reg_;


};

} // namespace


#endif
