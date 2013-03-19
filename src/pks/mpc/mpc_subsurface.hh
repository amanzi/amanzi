/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Interface for the derived MPC for coupling energy and water in the subsurface,
with freezing.

------------------------------------------------------------------------- */

#ifndef MPC_SUBSURFACE_HH_
#define MPC_SUBSURFACE_HH_

#include "mpc_coupled_cells.hh"

namespace Amanzi {

class MPCDelegateEWC;

class MPCSubsurface : public MPCCoupledCells {

 public:

  MPCSubsurface(Teuchos::ParameterList& plist,
                const Teuchos::RCP<TreeVector>& soln) :
      PKDefaultBase(plist, soln),
      MPCCoupledCells(plist,soln) {}

  // -- Initialize owned (dependent) variables.
  virtual void setup(const Teuchos::Ptr<State>& S);
  virtual void initialize(const Teuchos::Ptr<State>& S);

  virtual void set_states(const Teuchos::RCP<const State>& S,
                          const Teuchos::RCP<State>& S_inter,
                          const Teuchos::RCP<State>& S_next);

  virtual void commit_state(double dt, const Teuchos::RCP<State>& S);

  // update the predictor to be physically consistent
  virtual bool modify_predictor(double h, Teuchos::RCP<TreeVector> up);

  // updates the preconditioner
  virtual void update_precon(double t, Teuchos::RCP<const TreeVector> up, double h);

  // preconditioner application
  virtual void precon(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu);

 protected:
  bool modify_predictor_heuristic_(double h, Teuchos::RCP<TreeVector> up);

 protected:

  enum PredictorType {
    PREDICTOR_NONE = 0,
    PREDICTOR_HEURISTIC = 1,
    PREDICTOR_EWC = 2,
    PREDICTOR_SMART_EWC = 3,
  };

  enum PreconditionerType {
    PRECON_NONE = 0,
    PRECON_BLOCK_DIAGONAL = 1,
    PRECON_PICARD = 2,
    PRECON_EWC = 3,
    PRECON_SMART_EWC = 4,
  };

  // preconditioner methods
  PreconditionerType precon_type_;

  // prediction methods
  PredictorType predictor_type_;

  // EWC delegate
  Teuchos::RCP<MPCDelegateEWC> ewc_;

 private:
  // factory registration
  static RegisteredPKFactory<MPCSubsurface> reg_;

};


} // namespace

#endif


