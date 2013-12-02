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

#include "mpc_surface_subsurface_flux_coupler.hh"
#include "strong_mpc.hh"

namespace Amanzi {

namespace Operators { class MatrixMFD_Coupled_Surf; }
class MPCDelegateEWC;

class MPCPermafrost : public StrongMPC<MPCSurfaceSubsurfaceFluxCoupler> {

 public:
  MPCPermafrost(const Teuchos::RCP<Teuchos::ParameterList>& plist,
                Teuchos::ParameterList& FElist,
                const Teuchos::RCP<TreeVector>& soln) :
      PKDefaultBase(plist, FElist, soln),
      StrongMPC<MPCSurfaceSubsurfaceFluxCoupler>(plist, FElist, soln) {}


  virtual void setup(const Teuchos::Ptr<State>& S);
  virtual void initialize(const Teuchos::Ptr<State>& S);

  virtual void set_states(const Teuchos::RCP<const State>& S,
                          const Teuchos::RCP<State>& S_inter,
                          const Teuchos::RCP<State>& S_next);

  virtual void commit_state(double dt, const Teuchos::RCP<State>& S);

  // update the predictor to be physically consistent
  virtual bool modify_predictor(double h, Teuchos::RCP<TreeVector> up);

  // modify correction post NKA
  virtual bool modify_correction(double h, Teuchos::RCP<const TreeVector> res,
          Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> du);

  // applies preconditioner to u and returns the result in Pu
  virtual void precon(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu);

  // updates the preconditioner
  virtual void update_precon(double t, Teuchos::RCP<const TreeVector> up, double h);

 protected:
  // update the predictor to be physically consistent
  bool modify_predictor_for_source_on_ice_(double h, Teuchos::RCP<TreeVector> up);

  
 protected:
  Key dA_dy2_key_;
  Key dB_dy1_key_;
  Key A_key_;
  Key y2_key_;
  Key B_key_;
  Key y1_key_;

  bool decoupled_;

  enum PredictorType {
    PREDICTOR_NONE = 0,
    PREDICTOR_EWC,
    PREDICTOR_SMART_EWC
  };
  // prediction methods
  PredictorType predictor_type_;

  Teuchos::RCP<MPCSurfaceSubsurfaceFluxCoupler> coupled_flow_pk_;
  Teuchos::RCP<MPCSurfaceSubsurfaceFluxCoupler> coupled_energy_pk_;

  Teuchos::RCP<Operators::MatrixMFD_Coupled_Surf> mfd_surf_preconditioner_;

  // EWC delegate
  Teuchos::RCP<MPCDelegateEWC> ewc_;

  // meshes
  Teuchos::RCP<const AmanziMesh::Mesh> surf_mesh_;
  Teuchos::RCP<const AmanziMesh::Mesh> domain_mesh_;
  
  // debugger for dumping vectors
  Teuchos::RCP<Debugger> domain_db_;
  Teuchos::RCP<Debugger> surf_db_;
  
 private:
  // factory registration
  static RegisteredPKFactory<MPCPermafrost> reg_;

};

} // namespace

#endif
