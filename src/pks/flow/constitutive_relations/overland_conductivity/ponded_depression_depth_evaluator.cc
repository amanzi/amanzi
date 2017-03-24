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

  Key domain;
  
  if(!my_key_.empty())
    domain = getDomain(my_key_);
  else if (my_key_.empty())
    my_key_ = plist_.get<std::string>("ponded depresssion depth key", "surface_star-ponded_depression_depth");

  pd_key_ = plist_.get<std::string>("height key", getKey(domain,"ponded_depth"));
  dependencies_.insert(pd_key_);

  depr_depth_key_ = plist_.get<std::string>("depression depth key", getKey(domain,"depression_depth"));
  dependencies_.insert(depr_depth_key_);
  // depr_depth_ = plist_.get<double>("depression depth");//, 0.043);

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

  Epetra_MultiVector& res = *result->ViewComponent("cell",false);
  const Epetra_MultiVector& d_depth = *S->GetFieldData(pd_key_)->ViewComponent("cell",false);

  //  Teuchos::RCP<const CompositeVector> p_depth = S->GetFieldData(pdd_key_);
  //  const double dep& d_depth = *S->GetScalarData(dd_key_);
  
  const Epetra_MultiVector& depr_depth_v = *S->GetFieldData(depr_depth_key_)->ViewComponent("cell", false);

  assert(depr_depth_v[0][3] > 0.);

  int ncells = res.MyLength();
  for (int c=0; c!=ncells; ++c) {
    res[0][c] = d_depth[0][c] - depr_depth_v[0][c];
      }
  //  result->Update(1.0, *p_depth, -1.0, *d_depth, 0.0);
}


// This is hopefully never called?
void PondedDepressionDepthEvaluator::EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
        Key wrt_key, const Teuchos::Ptr<CompositeVector>& result) {

  result->PutScalar(1.0);
}

} //namespace
} //namespace
} //namespace
