/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  Evaluates the unfrozen effective depth.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "unfrozen_effective_depth_evaluator.hh"

namespace Amanzi {
namespace Flow {
namespace FlowRelations {

// registry of method
Utils::RegisteredFactory<FieldEvaluator,UnfrozenEffectiveDepthEvaluator> UnfrozenEffectiveDepthEvaluator::fac_("unfrozen effective depth");


UnfrozenEffectiveDepthEvaluator::UnfrozenEffectiveDepthEvaluator(Teuchos::ParameterList& plist) :
    SecondaryVariableFieldEvaluator(plist) {

  depth_key_ = plist_.get<std::string>("depth key", "ponded_depth");
  dependencies_.insert(depth_key_);

  uf_key_ = plist_.get<std::string>("unfrozen fraction key", "unfrozen_fraction");
  dependencies_.insert(uf_key_);

  if (my_key_ == std::string("")) {
    my_key_ = plist_.get<std::string>("unfrozen effective depth key",
            "unfrozen_effective_depth");
  }
  setLinePrefix(my_key_+std::string(" evaluator"));

}


UnfrozenEffectiveDepthEvaluator::UnfrozenEffectiveDepthEvaluator(const UnfrozenEffectiveDepthEvaluator& other) :
    SecondaryVariableFieldEvaluator(other),
    depth_key_(other.depth_key_),
    uf_key_(other.uf_key_) {}


Teuchos::RCP<FieldEvaluator>
UnfrozenEffectiveDepthEvaluator::Clone() const {
  return Teuchos::rcp(new UnfrozenEffectiveDepthEvaluator(*this));
}


// Required methods from SecondaryVariableFieldEvaluator
void UnfrozenEffectiveDepthEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result) {
  Teuchos::RCP<const CompositeVector> depth = S->GetFieldData(depth_key_);
  Teuchos::RCP<const CompositeVector> uf = S->GetFieldData(uf_key_);

  result->Multiply(1., *depth, *uf, 0.);
}


void UnfrozenEffectiveDepthEvaluator::EvaluateFieldPartialDerivative_(
    const Teuchos::Ptr<State>& S,
    Key wrt_key, const Teuchos::Ptr<CompositeVector>& result) {

  Teuchos::RCP<const CompositeVector> depth = S->GetFieldData(depth_key_);
  Teuchos::RCP<const CompositeVector> uf = S->GetFieldData(uf_key_);

  if (wrt_key == depth_key_) {
    *result = *uf;
  } else if (wrt_key == uf_key_) {
    *result = *depth;
  } else {
    ASSERT(0);
  }
}


} //namespace
} //namespace
} //namespace

