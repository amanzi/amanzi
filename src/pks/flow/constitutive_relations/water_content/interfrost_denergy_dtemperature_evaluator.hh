/*
  The interfrost denergy_dtemperature evaluator is an algebraic evaluator of a given model.

  Generated via evaluator_generator with:
Interfrost water content portion sl.
    
  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_FLOW_INTERFROST_DENERGY_DTEMPERATURE_EVALUATOR_HH_
#define AMANZI_FLOW_INTERFROST_DENERGY_DTEMPERATURE_EVALUATOR_HH_

#include "Factory.hh"
#include "secondary_variable_field_evaluator.hh"

namespace Amanzi {
namespace Flow {
namespace Relations {

class InterfrostDenergyDtemperatureModel;

class InterfrostDenergyDtemperatureEvaluator : public SecondaryVariableFieldEvaluator {

 public:
  explicit
  InterfrostDenergyDtemperatureEvaluator(Teuchos::ParameterList& plist);
  InterfrostDenergyDtemperatureEvaluator(const InterfrostDenergyDtemperatureEvaluator& other);

  virtual Teuchos::RCP<FieldEvaluator> Clone() const;

  // Required methods from SecondaryVariableFieldEvaluator
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& result);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result);

  Teuchos::RCP<InterfrostDenergyDtemperatureModel> get_model() { return model_; }

 protected:
  void InitializeFromPlist_();

  Key phi_key_;
  Key sl_key_;
  Key nl_key_;
  Key si_key_;
  Key ni_key_;
  Key rhos_key_;
  Key T_key_;

  Teuchos::RCP<InterfrostDenergyDtemperatureModel> model_;

 private:
  static Utils::RegisteredFactory<FieldEvaluator,InterfrostDenergyDtemperatureEvaluator> reg_;

};

} //namespace
} //namespace
} //namespace

#endif
