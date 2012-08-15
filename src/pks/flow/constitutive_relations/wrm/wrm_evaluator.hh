/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  The WRM Evaluator simply calls the WRM with the correct arguments.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_FLOW_RELATIONS_WRM_EVALUATOR_
#define AMANZI_FLOW_RELATIONS_WRM_EVALUATOR_

#include "secondary_variables_field_evaluator.hh"

namespace Amanzi {
namespace Flow {
namespace FlowRelations {

class WRM; // forward declaration

class WRMEvaluator : public SecondaryVariablesFieldEvaluator {

 public:
  // constructor format for all derived classes
  explicit
  WRMEvaluator(Teuchos::ParameterList& wrm_plist);
  WRMEvaluator(Teuchos::ParameterList& wrm_plist, const Teuchos::RCP<WRM>& wrm);
  WRMEvaluator(const WRMEvaluator& other);

  // Required methods from SecondaryVariableFieldEvaluator
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const std::vector<Teuchos::Ptr<CompositeVector> >& results) = 0;
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const std::vector<Teuchos::Ptr<CompositeVector> > & results) = 0;

  Teuchos::RCP<WRM> get_WRM() { return wrm_; }

 protected:

  Teuchos::ParameterList wrm_plist_;
  Teuchos::RCP<WRM> wrm_;
};

} //namespace
} //namespace
} //namespace

#endif
