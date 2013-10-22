/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  Evaluates the conductivity of surface flow.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_FLOWRELATIONS_OVERLAND_CONDUCTIVITY_EVALUATOR_
#define AMANZI_FLOWRELATIONS_OVERLAND_CONDUCTIVITY_EVALUATOR_

#include "secondary_variable_field_evaluator.hh"

namespace Amanzi {
namespace Flow {
namespace FlowRelations {

class OverlandConductivityModel;

class OverlandConductivityEvaluator : public SecondaryVariableFieldEvaluator {

 public:
  OverlandConductivityEvaluator(Teuchos::ParameterList& plist);
  OverlandConductivityEvaluator(const OverlandConductivityEvaluator& other);
  Teuchos::RCP<FieldEvaluator> Clone() const;

  // Required methods from SecondaryVariableFieldEvaluator
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& result);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result);

  Teuchos::RCP<OverlandConductivityModel> get_Model() { return model_; }

private:
  Teuchos::RCP<OverlandConductivityModel> model_;

  Key depth_key_;
  Key slope_key_;
  Key coef_key_;

};

} //namespace
} //namespace
} //namespace

#endif

