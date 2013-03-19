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
#define MPC_FROZEN_PREC_COUPLED_FLOW_ENERGY_HH_

#include "mpc_coupled_cells.hh"

namespace Amanzi {

class PermafrostModel;

class MPCPermafrost : public MPCCoupledCells {

 public:
  MPCPermafrost(Teuchos::ParameterList& plist,
                             const Teuchos::RCP<TreeVector>& soln) :
      PKDefaultBase(plist, soln),
      MPCCoupledCells(plist, soln) {}


  virtual void setup(const Teuchos::Ptr<State>& S);

  // applies preconditioner to u and returns the result in Pu
  virtual void precon(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu);

 protected:
  Teuchos::RCP<MPCSurfaceSubsurfaceFluxCoupler> coupled_flow_pk_;
  Teuchos::RCP<MPCSurfaceSubsurfaceFluxCoupler> coupled_energy_pk_;

 private:
  // factory registration
  static RegisteredPKFactory<MPCPermafrost> reg_;

};
} // namespace

#endif
