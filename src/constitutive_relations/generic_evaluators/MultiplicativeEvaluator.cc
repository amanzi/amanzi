/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  MultiplicativeEvaluator is the generic evaluator for multipying two vectors.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "MultiplicativeEvaluator.hh"

namespace Amanzi {
namespace Relations {

MultiplicativeEvaluator::MultiplicativeEvaluator(Teuchos::ParameterList& plist) :
    SecondaryVariableFieldEvaluator(plist) {}

MultiplicativeEvaluator(const MultiplicativeEvaluator& other) :
  virtual Teuchos::RCP<FieldEvaluator> Clone() const;

  // Required methods from SecondaryVariableFieldEvaluator
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& result);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result);


} // namespace
} // namespace
