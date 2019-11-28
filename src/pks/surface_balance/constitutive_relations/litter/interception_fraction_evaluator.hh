/*
  The interception fraction evaluator is an algebraic evaluator of a given model.

  Generated via evaluator_generator with:

    
  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_SURFACEBALANCE_INTERCEPTION_FRACTION_EVALUATOR_HH_
#define AMANZI_SURFACEBALANCE_INTERCEPTION_FRACTION_EVALUATOR_HH_

#include "Factory.hh"
#include "secondary_variable_field_evaluator.hh"

namespace Amanzi {
namespace SurfaceBalance {
namespace Relations {

class InterceptionFractionModel;

class InterceptionFractionEvaluator : public SecondaryVariableFieldEvaluator {

 public:
  explicit
  InterceptionFractionEvaluator(Teuchos::ParameterList& plist);
  InterceptionFractionEvaluator(const InterceptionFractionEvaluator& other);

  virtual Teuchos::RCP<FieldEvaluator> Clone() const;

  // Required methods from SecondaryVariableFieldEvaluator
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& result);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result);

  Teuchos::RCP<InterceptionFractionModel> get_model() { return model_; }

 protected:
  void InitializeFromPlist_();

  Key ai_key_;

  Teuchos::RCP<InterceptionFractionModel> model_;

 private:
  static Utils::RegisteredFactory<FieldEvaluator,InterceptionFractionEvaluator> reg_;

};

} //namespace
} //namespace
} //namespace

#endif