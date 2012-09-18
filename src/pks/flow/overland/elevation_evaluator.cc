/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  The elevation evaluator gets the surface elevation, slope, and updates pres + elev.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "elevation_evaluator.hh"

namespace Amanzi {
namespace Flow {
namespace FlowRelations {

ElevationEvaluator::ElevationEvaluator(Teuchos::ParameterList& plist) :
    SecondaryVariablesFieldEvaluator(plist) {
  my_keys_.push_back(plist_.get<std::string>("elevation key", "elevation"));
  my_keys_.push_back(plist_.get<std::string>("slope magnitude key", "slope_magnitude"));
  my_keys_.push_back(plist_.get<std::string>("potential key", "pres_elev"));
  setLinePrefix(my_keys_[0]+std::string(" evaluator"));

  pres_key_ = plist_.get<std::string>("height key", "ponded_depth");
  dependencies_.insert(pres_key_);
}

void ElevationEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const std::vector<Teuchos::Ptr<CompositeVector> >& results) {
  // If they haven't been done yet, update slope and elevation.
  if (!updated_once_) {
    EvaluateElevationAndSlope_(S, results);
    updated_once_ = true;
  }

  // update pressure + elevation
  Teuchos::RCP<const CompositeVector> pres = S->GetFieldData(pres_key_);
  results[2]->Update(1.0, *results[0], 1.0, *pres, 0.0);
}


// This is hopefully never called?
void ElevationEvaluator::EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
        Key wrt_key, const std::vector<Teuchos::Ptr<CompositeVector> >& results) {
  ASSERT(0);

  ASSERT(wrt_key == pres_key_);

  results[0]->PutScalar(0.0);
  results[1]->PutScalar(0.0);
  results[2]->PutScalar(1.0);
}

} //namespace
} //namespace
} //namespace
