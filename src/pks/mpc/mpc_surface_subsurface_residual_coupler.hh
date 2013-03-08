/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Interface for the derived MPC for water coupling between surface and subsurface.

In this method, a Neumann BC is used on the subsurface boundary for
the operator, but the preconditioner is for the residual system with no
extra unknowns.  On the surface, the TPFA is used, resulting in a
subsurface-face-only Schur complement that captures all terms.

------------------------------------------------------------------------- */

#ifndef PKS_MPC_SURFACE_SUBSURFACE_RESIDUAL_COUPLER_HH_
#define PKS_MPC_SURFACE_SUBSURFACE_RESIDUAL_COUPLER_HH_

#include "mpc_surface_subsurface_coupler.hh"

namespace Amanzi {

namespace Operators { class MatrixMFD_Surf; }

class MPCSurfaceSubsurfaceResidualCoupler : public MPCSurfaceSubsurfaceCoupler {

 public:
  MPCSurfaceSubsurfaceResidualCoupler(Teuchos::ParameterList& plist,
          const Teuchos::RCP<TreeVector>& soln) :
      PKDefaultBase(plist, soln),
      MPCSurfaceSubsurfaceCoupler(plist, soln) {
    surf_c0_ = plist_.get<int>("surface debug cell 0", 0);
    surf_c1_ = plist_.get<int>("surface debug cell 1", 1);
    damping_coef_ = plist.get<double>("damping coefficient", -1.);
    if (damping_coef_ > 0.) {
      damping_cutoff_ = plist.get<double>("damping cutoff", 0.1);
    }

    modify_predictor_heuristic_ = plist.get<bool>("modify predictor with heuristic", false);
  }

  // -- Setup data.
  virtual void setup(const Teuchos::Ptr<State>& S);

  // evaluate the residual
  virtual void fun(double t_old, double t_new, Teuchos::RCP<TreeVector> u_old,
           Teuchos::RCP<TreeVector> u_new, Teuchos::RCP<TreeVector> g);

  // applies preconditioner to u and returns the result in Pu
  virtual void precon(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu);

  // updates the preconditioner
  virtual void update_precon(double t, Teuchos::RCP<const TreeVector> up, double h);

  virtual bool modify_predictor(double h, Teuchos::RCP<TreeVector> u);

 protected:
  Teuchos::RCP<Operators::MatrixMFD_Surf> preconditioner_;
  Teuchos::RCP<Operators::MatrixMFD_TPFA> surf_preconditioner_;

  int surf_c0_;
  int surf_c1_;
  double damping_coef_;
  double damping_cutoff_;
  bool modify_predictor_heuristic_;

 private:
  // factory registration
  static RegisteredPKFactory<MPCSurfaceSubsurfaceResidualCoupler> reg_;

};

} // namespace

#endif

