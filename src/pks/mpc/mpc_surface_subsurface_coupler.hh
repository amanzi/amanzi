/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Interface for the derived MPC for water coupling between surface and subsurface.
------------------------------------------------------------------------- */

#ifndef PKS_MPC_SURFACE_SUBSURFACE_COUPLER_HH_
#define PKS_MPC_SURFACE_SUBSURFACE_COUPLER_HH_

#include "strong_mpc.hh"

namespace Amanzi {

class MPCSurfaceSubsurfaceCoupler : public StrongMPC {

 public:
  MPCSurfaceSubsurfaceCoupler(Teuchos::ParameterList& plist,
          const Teuchos::RCP<TreeVector>& soln);

  void fun(double t_old, double t_new, Teuchos::RCP<TreeVector> u_old,
           Teuchos::RCP<TreeVector> u_new, Teuchos::RCP<TreeVector> g);

  // Virtual destructor
  virtual ~MPCSurfaceSubsurfaceCoupler() {}

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
  static RegisteredPKFactory<MPCSurfaceSubsurfaceCoupler> reg_;


};

} // namespace


#endif
