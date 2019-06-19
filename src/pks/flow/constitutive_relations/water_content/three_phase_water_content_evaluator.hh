/*
  The three phase water content evaluator is an algebraic evaluator of a given model.

  Generated via evaluator_generator with:
Water content for a three-phase, gas+liquid+ice evaluator.
    
  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_FLOW_THREE_PHASE_WATER_CONTENT_EVALUATOR_HH_
#define AMANZI_FLOW_THREE_PHASE_WATER_CONTENT_EVALUATOR_HH_

#include "Factory.hh"
#include "secondary_variable_field_evaluator.hh"

namespace Amanzi {
namespace Flow {
namespace Relations {

class ThreePhaseWaterContentModel;

class ThreePhaseWaterContentEvaluator : public SecondaryVariableFieldEvaluator {

 public:
  explicit
  ThreePhaseWaterContentEvaluator(Teuchos::ParameterList& plist);
  ThreePhaseWaterContentEvaluator(const ThreePhaseWaterContentEvaluator& other);

  virtual Teuchos::RCP<FieldEvaluator> Clone() const;

  // Required methods from SecondaryVariableFieldEvaluator
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& result);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result);

  Teuchos::RCP<ThreePhaseWaterContentModel> get_model() { return model_; }

 protected:
  void InitializeFromPlist_();

  Key phi_key_;
  Key sl_key_;
  Key nl_key_;
  Key si_key_;
  Key ni_key_;
  Key sg_key_;
  Key ng_key_;
  Key omega_key_;
  Key cv_key_;

  Teuchos::RCP<ThreePhaseWaterContentModel> model_;

 private:
  static Utils::RegisteredFactory<FieldEvaluator,ThreePhaseWaterContentEvaluator> reg_;

};

} //namespace
} //namespace
} //namespace

#endif