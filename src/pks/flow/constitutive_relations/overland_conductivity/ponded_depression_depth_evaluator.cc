/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  The elevation evaluator gets the surface elevation, slope, and updates pres + elev.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "ponded_depression_depth_evaluator.hh"

namespace Amanzi {
namespace Flow {
namespace FlowRelations {

PondedDepressionDepthEvaluator::PondedDepressionDepthEvaluator(Teuchos::ParameterList& plist) :
    SecondaryVariableFieldEvaluator(plist) {

  Key domain = Keys::getDomain(my_key_);
  /*
  if(!my_key_.empty())
  else if (my_key_.empty())
    my_key_ = plist_.get<std::string>("ponded depresssion depth key", "surface_star-ponded_depression_depth");
  */
  
  pd_key_ = Keys::readKey(plist_, domain, "height", "ponded_depth");
  dependencies_.insert(pd_key_);

  depr_depth_key_ = Keys::readKey(plist_, domain, "depression depth", "depression_depth");
  dependencies_.insert(depr_depth_key_);
}


PondedDepressionDepthEvaluator::PondedDepressionDepthEvaluator(const PondedDepressionDepthEvaluator& other) :
    SecondaryVariableFieldEvaluator(other),
    pd_key_(other.pd_key_),
    depr_depth_key_(other.depr_depth_key_) {};

Teuchos::RCP<FieldEvaluator>
PondedDepressionDepthEvaluator::Clone() const {
  return Teuchos::rcp(new PondedDepressionDepthEvaluator(*this));
}


void PondedDepressionDepthEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result) {

  auto& res = *result->ViewComponent("cell",false);
  const auto& d_depth = *S->GetFieldData(pd_key_)->ViewComponent("cell",false);
  const auto& depr_depth_v = *S->GetFieldData(depr_depth_key_)->ViewComponent("cell", false);

  int ncells = res.MyLength();
  res.Update(1., d_depth, -1., depr_depth_v, 0.);
}


// This is hopefully never called?
void PondedDepressionDepthEvaluator::EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
        Key wrt_key, const Teuchos::Ptr<CompositeVector>& result) {
  if (wrt_key == pd_key_) {
    result->PutScalar(1.0);
  } else if (wrt_key == depr_depth_key_) {
    result->PutScalar(-1.);
  } else {
    AMANZI_ASSERT(false);
  }
}

} //namespace
} //namespace
} //namespace
