/*
  The interfrost sl water content evaluator is an algebraic evaluator of a given model.

  Generated via evaluator_generator with:
Interfrost water content portion sl.
    
  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_FLOW_INTERFROST_SL_WC_EVALUATOR_HH_
#define AMANZI_FLOW_INTERFROST_SL_WC_EVALUATOR_HH_

#include "Factory.hh"
#include "secondary_variable_field_evaluator.hh"

namespace Amanzi {
namespace Flow {
namespace Relations {

class InterfrostSlWcModel;

class InterfrostSlWcEvaluator : public SecondaryVariableFieldEvaluator {

 public:
  explicit
  InterfrostSlWcEvaluator(Teuchos::ParameterList& plist);
  InterfrostSlWcEvaluator(const InterfrostSlWcEvaluator& other);

  virtual Teuchos::RCP<FieldEvaluator> Clone() const;

  // Required methods from SecondaryVariableFieldEvaluator
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& result);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result);

  Teuchos::RCP<InterfrostSlWcModel> get_model() { return model_; }

 protected:
  void InitializeFromPlist_();

  Key phi_key_;
  Key sl_key_;
  Key nl_key_;
  Key ni_key_;
  Key cv_key_;

  Teuchos::RCP<InterfrostSlWcModel> model_;

 private:
  static Utils::RegisteredFactory<FieldEvaluator,InterfrostSlWcEvaluator> reg_;

};

} //namespace
} //namespace
} //namespace

#endif
