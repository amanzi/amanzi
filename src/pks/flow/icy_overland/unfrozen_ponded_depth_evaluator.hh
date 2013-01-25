/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  Unfrozen ponded depth

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_FLOWRELATIONS_UNFROZEN_PONDED_DEPTH_EVALUATOR_
#define AMANZI_FLOWRELATIONS_UNFROZEN_PONDED_DEPTH_EVALUATOR_

#include "secondary_variable_field_evaluator.hh"


namespace Amanzi {
namespace Flow {
namespace FlowRelations {

class UnfrozenPondedDepthEvaluator : public SecondaryVariableFieldEvaluator {
 public:
  explicit
  UnfrozenPondedDepthEvaluator(Teuchos::ParameterList& plist);
  UnfrozenPondedDepthEvaluator(const UnfrozenPondedDepthEvaluator& other);
  virtual Teuchos::RCP<FieldEvaluator> Clone() const;

  // Required methods from SecondaryVariableFieldEvaluator
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& result);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result);

 protected:
  Key depth_key_;
  Key temp_key_;
};

} //namespace
} //namespace
} //namespace

#endif

