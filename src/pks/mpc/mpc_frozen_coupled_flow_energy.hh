/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Interface for the derived MPC for flow and energy.  This couples using a
block-diagonal coupler.
------------------------------------------------------------------------- */

#ifndef MPC_FROZEN_PREC_COUPLED_FLOW_ENERGY_HH_
#define MPC_FROZEN_PREC_COUPLED_FLOW_ENERGY_HH_

#include "mpc_coupled_cells.hh"

namespace Amanzi {

class PermafrostModel;
namespace Flow { class Richards; }
namespace Energy { class TwoPhase; }

class MPCFrozenCoupledFlowEnergy : public MPCCoupledCells {

public:
  MPCFrozenCoupledFlowEnergy(Teuchos::ParameterList& plist,
                             const Teuchos::RCP<TreeVector>& soln) :
      PKDefaultBase(plist, soln),
      MPCCoupledCells(plist, soln),
      consistent_by_average_(false) {}

  // Virtual destructor
  virtual ~MPCFrozenCoupledFlowEnergy() {}

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
  virtual bool modify_predictor_heuristic_(double h, Teuchos::RCP<TreeVector> up);
  virtual bool modify_predictor_ewc_(double h, Teuchos::RCP<TreeVector> up);
  virtual bool modify_predictor_smart_ewc_(double h, Teuchos::RCP<TreeVector> up);

  void InitializeModels_(const Teuchos::Ptr<State>& S);
  void InitialConditionFromFrozenColumn_(const Teuchos::Ptr<State>& S);

  virtual void precon_ewc_(Teuchos::RCP<const TreeVector> u,
                           Teuchos::RCP<TreeVector> Pu);
  virtual void precon_smart_ewc_(Teuchos::RCP<const TreeVector> u,
          Teuchos::RCP<TreeVector> Pu);
  virtual void update_precon_ewc_(double t, Teuchos::RCP<const TreeVector> up, double h);

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

  // PKs
  Teuchos::RCP<Flow::Richards> flow_pk_;
  Teuchos::RCP<Energy::TwoPhase> energy_pk_;

  // preconditioner methods
  PreconditionerType precon_type_;

  // prediction methods
  PredictorType predictor_type_;
  bool consistent_by_average_;

  // -- heuristic
  bool modify_thaw_to_prev_;

  // -- ewc
  Teuchos::RCP<PermafrostModel> model_;
  double cusp_size_T_freezing_;
  double cusp_size_T_thawing_;
  Teuchos::RCP<State> S_work_;
  std::vector<int> cells_to_track_;
  std::vector<WhetStone::Tensor> jac_;

private:
  // factory registration
  static RegisteredPKFactory<MPCFrozenCoupledFlowEnergy> reg_;

};
} // namespace

#endif
