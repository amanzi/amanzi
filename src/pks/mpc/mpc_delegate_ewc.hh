/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Interface for EWC, a helper class that does projections and preconditioners in
energy/water-content space instead of temperature/pressure space.
------------------------------------------------------------------------- */

#ifndef MPC_DELEGATE_EWC_HH_
#define MPC_DELEGATE_EWC_HH_

#include "tensor.hpp"
#include "state.hh"
#include "tree_vector.hh"

namespace Amanzi {

class EWCModel;

class MPCDelegateEWC {

 public:

  MPCDelegateEWC(Teuchos::ParameterList& plist);

  virtual void setup(const Teuchos::Ptr<State>& S);
  virtual void initialize(const Teuchos::Ptr<State>& S);
  virtual void set_states(const Teuchos::RCP<const State>& S,
                          const Teuchos::RCP<State>& S_inter,
                          const Teuchos::RCP<State>& S_next);
  virtual void commit_state(double dt, const Teuchos::RCP<State>& S);
  virtual bool modify_predictor(double h, Teuchos::RCP<TreeVector> up);
  virtual void update_precon(double t, Teuchos::RCP<const TreeVector> up, double h);
  virtual void precon(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu);

  void set_model(const Teuchos::RCP<EWCModel>& model) { model_ = model; }

 protected:
  virtual bool modify_predictor_ewc_(double h, Teuchos::RCP<TreeVector> up);
  virtual bool modify_predictor_smart_ewc_(double h, Teuchos::RCP<TreeVector> up);

  virtual void precon_smart_ewc_(Teuchos::RCP<const TreeVector> u,
          Teuchos::RCP<TreeVector> Pu);
  virtual void precon_ewc_(Teuchos::RCP<const TreeVector> u,
          Teuchos::RCP<TreeVector> Pu);

  virtual void update_precon_ewc_(double t, Teuchos::RCP<const TreeVector> up, double h);



 protected:
  Teuchos::ParameterList plist_;

  // model
  Teuchos::RCP<EWCModel> model_;

  enum PredictorType {
    PREDICTOR_NONE = 0,
    PREDICTOR_EWC,
    PREDICTOR_SMART_EWC
  };

  enum PreconditionerType {
    PRECON_NONE = 0,
    PRECON_EWC,
    PRECON_SMART_EWC,
  };

  // control flags
  PreconditionerType precon_type_;
  PredictorType predictor_type_;

  // extra data
  std::vector<WhetStone::Tensor> jac_;
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  Teuchos::RCP<Epetra_MultiVector> wc_prev2_;
  Teuchos::RCP<Epetra_MultiVector> e_prev2_;
  double time_prev2_;

  // parameters for heuristic
  double cusp_size_T_freezing_;
  double cusp_size_T_thawing_;

  // states
  Teuchos::RCP<State> S_next_;
  Teuchos::RCP<State> S_inter_;
  Teuchos::RCP<const State> S_;

  // keys
  Key pres_key_;
  Key temp_key_;
  Key e_key_;
  Key wc_key_;
  Key poro_key_;
  Key cv_key_;

};

} // namespace

#endif
