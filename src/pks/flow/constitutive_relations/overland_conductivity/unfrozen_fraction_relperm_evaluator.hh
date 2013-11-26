/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  Evaluates the unfrozen fraction model.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_FLOWRELATIONS_UNFROZEN_FRACTION_KR_EVALUATOR_
#define AMANZI_FLOWRELATIONS_UNFROZEN_FRACTION_KR_EVALUATOR_

#include "factory.hh"
#include "secondary_variable_field_evaluator.hh"

namespace Amanzi {
namespace Flow {
namespace FlowRelations {

class UnfrozenFractionRelPermModel;

class UnfrozenFractionRelPermEvaluator : public SecondaryVariableFieldEvaluator {

 public:
  UnfrozenFractionRelPermEvaluator(Teuchos::ParameterList& plist);
  UnfrozenFractionRelPermEvaluator(const UnfrozenFractionRelPermEvaluator& other);
  Teuchos::RCP<FieldEvaluator> Clone() const;

  // Required methods from SecondaryVariableFieldEvaluator
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& result);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result);

  Teuchos::RCP<UnfrozenFractionRelPermModel> get_Model() { return model_; }

protected:
  Teuchos::RCP<UnfrozenFractionRelPermModel> model_;
  Key uf_key_;
  Key h_key_;

private:
  static Utils::RegisteredFactory<FieldEvaluator,UnfrozenFractionRelPermEvaluator> fac_;


};

} //namespace
} //namespace
} //namespace

#endif

