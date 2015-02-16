/*
  License: see $AMANZI_DIR/COPYRIGHT
  Authors: Ethan Coon

  Temporary wrapper converting the Transport_PK, which inherits from 
  BDFFnBase<CompositeVector>, to use TreeVectors.
*/

#ifndef AMANZI_TRANSPORT_PK_WRAPPER_HH_
#define AMANZI_TRANSPORT_PK_WRAPPER_HH_

#include "Teuchos_RCP.hpp"

#include "PK.hh"
#include "Transport_PK.hh"
#include "PK_Factory.hh"

namespace Amanzi {
namespace Transport {

class Transport_PK_Wrapper : public PK {
 public:
  Transport_PK_Wrapper(Teuchos::ParameterList& pk_tree,
                       const Teuchos::RCP<Teuchos::ParameterList>& global_list,
                       const Teuchos::RCP<State>& S,
                       const Teuchos::RCP<TreeVector>& soln);

  // Setup
  virtual void Setup() {
    pk_->InitializeFields();
  }
  
  // Initialize owned (dependent) variables.
  virtual void Initialize() {
    pk_->Initialize(S_.ptr());
  }

  // Choose a time step compatible with physics.
  virtual double get_dt() {
    return pk_->get_dt();
  }

  virtual void set_dt(double dt) {};

  // Advance from state S0 to state S1 at time S0.time + dt.
  // Due to Transport PK / MPC conflict (FIXME when MPC will be upgraded)
  //  virtual int Advance(double dt, double& dt_actual) = 0;
  virtual bool AdvanceStep(double t_old, double t_new);

  // Commit any secondary (dependent) variables.
  virtual void CommitStep(double t_old, double t_new) {
    pk_->CommitState(t_new-t_old, S_.ptr());
  }

  // Calculate any diagnostics prior to doing vis
  virtual void CalculateDiagnostics() {};

  virtual std::string name() {
    return pk_->name();
  }

  Teuchos::RCP<CompositeVector> total_component_concentration() {
    return pk_->total_component_concentration();
  }

 protected:
  std::vector<std::string> comp_names_;
  Teuchos::RCP<Teuchos::ParameterList> glist_;
  Teuchos::ParameterList ti_list_;
  Teuchos::RCP<Transport_PK> pk_;
  Teuchos::RCP<TreeVector> soln_;
  Teuchos::RCP<State> S_;

 private:
  // factory registration
  static RegisteredPKFactory<Transport_PK_Wrapper> reg_;
    
};

}  // namespace Transport
}  // namespace Amanzi

#endif
