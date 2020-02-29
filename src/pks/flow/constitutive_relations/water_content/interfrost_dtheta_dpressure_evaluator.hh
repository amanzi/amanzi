/*
  The interfrost dtheta_dpressure evaluator is an algebraic evaluator of a given model.

  Generated via evaluator_generator with:
Interfrost water content portion sl.
    
  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_FLOW_INTERFROST_DTHETA_DPRESSURE_EVALUATOR_HH_
#define AMANZI_FLOW_INTERFROST_DTHETA_DPRESSURE_EVALUATOR_HH_

#include "Factory.hh"
#include "secondary_variable_field_evaluator.hh"

namespace Amanzi {
namespace Flow {
namespace Relations {

class InterfrostDthetaDpressureModel;

class InterfrostDthetaDpressureEvaluator : public SecondaryVariableFieldEvaluator {

 public:
  explicit
  InterfrostDthetaDpressureEvaluator(Teuchos::ParameterList& plist);
  InterfrostDthetaDpressureEvaluator(const InterfrostDthetaDpressureEvaluator& other);

  virtual Teuchos::RCP<FieldEvaluator> Clone() const;

  // Required methods from SecondaryVariableFieldEvaluator
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& result);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result);

  Teuchos::RCP<InterfrostDthetaDpressureModel> get_model() { return model_; }

 protected:
  void InitializeFromPlist_();

  Key nl_key_;
  Key sl_key_;
  Key phi_key_;

  Teuchos::RCP<InterfrostDthetaDpressureModel> model_;

 private:
  static Utils::RegisteredFactory<FieldEvaluator,InterfrostDthetaDpressureEvaluator> reg_;

};

} //namespace
} //namespace
} //namespace

#endif
