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

#include "strong_mpc.hh"

namespace Amanzi {

class MPCSurfaceSubsurfaceDirichletCoupler : public StrongMPC {

 public:
  MPCSurfaceSubsurfaceDirichletCoupler(Teuchos::ParameterList& plist,
          const Teuchos::RCP<TreeVector>& soln);

  void fun(double t_old, double t_new, Teuchos::RCP<TreeVector> u_old,
           Teuchos::RCP<TreeVector> u_new, Teuchos::RCP<TreeVector> g);

  // Virtual destructor
  virtual ~MPCSurfaceSubsurfaceDirichletCoupler() {}

  virtual void changed_solution();

  // applies preconditioner to u and returns the result in Pu
  virtual void precon(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu);

  // updates the preconditioner
  virtual void update_precon(double t, Teuchos::RCP<const TreeVector> up, double h);


 protected:
  Key surface_mesh_key_;
  Key domain_field_;
  Key surface_field_;
  Key domain_pk_name_;
  Key surface_pk_name_;

 private:
  // factory registration
  static RegisteredPKFactory<MPCSurfaceSubsurfaceDirichletCoupler> reg_;


};

} // namespace


#endif
