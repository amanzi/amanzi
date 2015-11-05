/*
  The richards water content with vapor evaluator is an algebraic evaluator of a given model.

  Generated via evaluator_generator with:
Richards water content with vapor.
    
  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_FLOW_RICHARDS_WATER_CONTENT_WITH_VAPOR_EVALUATOR_HH_
#define AMANZI_FLOW_RICHARDS_WATER_CONTENT_WITH_VAPOR_EVALUATOR_HH_

#include "factory.hh"
#include "secondary_variable_field_evaluator.hh"

namespace Amanzi {
namespace Flow {
namespace Relations {

class RichardsWaterContentWithVaporModel;

class RichardsWaterContentWithVaporEvaluator : public SecondaryVariableFieldEvaluator {

 public:
  explicit
  RichardsWaterContentWithVaporEvaluator(Teuchos::ParameterList& plist);
  RichardsWaterContentWithVaporEvaluator(const RichardsWaterContentWithVaporEvaluator& other);

  virtual Teuchos::RCP<FieldEvaluator> Clone() const;

  // Required methods from SecondaryVariableFieldEvaluator
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& result);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result);

  Teuchos::RCP<RichardsWaterContentWithVaporModel> get_model() { return model_; }

 protected:
  void InitializeFromPlist_();

  Key phi_key_;
  Key sl_key_;
  Key sg_key_;
  Key nl_key_;
  Key ng_key_;
  Key omega_key_;
  Key cv_key_;

  Teuchos::RCP<RichardsWaterContentWithVaporModel> model_;

 private:
  static Utils::RegisteredFactory<FieldEvaluator,RichardsWaterContentWithVaporEvaluator> reg_;

};

} //namespace
} //namespace
} //namespace

#endif
