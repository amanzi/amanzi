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
  Transport_PK_Wrapper(const Teuchos::RCP<Teuchos::ParameterList>& plist,
                  Teuchos::ParameterList& glist,
                  const Teuchos::RCP<TreeVector>& soln);

  // Setup
  virtual void Setup(const Teuchos::Ptr<State>& S) {}
  
  // Initialize owned (dependent) variables.
  virtual void Initialize(const Teuchos::Ptr<State>& S) {
    pk_->Initialize(S);
  }

  // Choose a time step compatible with physics.
  virtual double get_dt() {
    return pk_->get_dt();
  }

  // Advance from state S0 to state S1 at time S0.time + dt.
  // Due to Transport PK / MPC conflict (FIXME when MPC will be upgraded)
  //  virtual int Advance(double dt, double& dt_actual) = 0;
  virtual bool Advance(double dt);

  // Commit any secondary (dependent) variables.
  virtual void CommitState(double dt, const Teuchos::Ptr<State>& S) {
    pk_->CommitState(dt, S);
  }

  // Calculate any diagnostics prior to doing vis
  virtual void CalculateDiagnostics(const Teuchos::Ptr<State>& S) {}

  virtual void SetState(const Teuchos::RCP<State>& S);

  virtual std::string name() {
    return pk_->name();
  }

 protected:
  std::vector<std::string> comp_names_;
  Teuchos::ParameterList glist_;
  Teuchos::RCP<Transport_PK> pk_;
  Teuchos::RCP<TreeVector> soln_;

 private:
  // factory registration
  static RegisteredPKFactory<Transport_PK_Wrapper> reg_;
    
};

} // namespace
} // namespace

#endif
