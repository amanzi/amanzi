/*
  The liquid + gas water content evaluator is an algebraic evaluator of a given model.

  Generated via evaluator_generator with:
Water content for a two-phase, liquid+water vapor evaluator.
    
  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_FLOW_LIQUID_GAS_WATER_CONTENT_EVALUATOR_HH_
#define AMANZI_FLOW_LIQUID_GAS_WATER_CONTENT_EVALUATOR_HH_

#include "Factory.hh"
#include "secondary_variable_field_evaluator.hh"

namespace Amanzi {
namespace Flow {
namespace Relations {

class LiquidGasWaterContentModel;

class LiquidGasWaterContentEvaluator : public SecondaryVariableFieldEvaluator {

 public:
  explicit
  LiquidGasWaterContentEvaluator(Teuchos::ParameterList& plist);
  LiquidGasWaterContentEvaluator(const LiquidGasWaterContentEvaluator& other);

  virtual Teuchos::RCP<FieldEvaluator> Clone() const;

  // Required methods from SecondaryVariableFieldEvaluator
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& result);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result);

  Teuchos::RCP<LiquidGasWaterContentModel> get_model() { return model_; }

 protected:
  void InitializeFromPlist_();

  Key phi_key_;
  Key sl_key_;
  Key nl_key_;
  Key sg_key_;
  Key ng_key_;
  Key omega_key_;
  Key cv_key_;

  Teuchos::RCP<LiquidGasWaterContentModel> model_;

 private:
  static Utils::RegisteredFactory<FieldEvaluator,LiquidGasWaterContentEvaluator> reg_;

};

} //namespace
} //namespace
} //namespace

#endif