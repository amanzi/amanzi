/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  Evaluates the unfrozen fraction model.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_FLOWRELATIONS_UNFROZEN_FRACTION_EVALUATOR_
#define AMANZI_FLOWRELATIONS_UNFROZEN_FRACTION_EVALUATOR_

#include "Factory.hh"
#include "secondary_variable_field_evaluator.hh"

namespace Amanzi {
namespace Flow {

class UnfrozenFractionModel;

class UnfrozenFractionEvaluator : public SecondaryVariableFieldEvaluator {

 public:
  UnfrozenFractionEvaluator(Teuchos::ParameterList& plist);
  UnfrozenFractionEvaluator(const UnfrozenFractionEvaluator& other);
  Teuchos::RCP<FieldEvaluator> Clone() const;

  // Required methods from SecondaryVariableFieldEvaluator
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& result);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result);

  Teuchos::RCP<const UnfrozenFractionModel> get_Model() const { return model_; }
  Teuchos::RCP<UnfrozenFractionModel> get_Model() { return model_; }

protected:
  Teuchos::RCP<UnfrozenFractionModel> model_;
  Key temp_key_;

private:
  static Utils::RegisteredFactory<FieldEvaluator,UnfrozenFractionEvaluator> fac_;


};

} //namespace
} //namespace

#endif

