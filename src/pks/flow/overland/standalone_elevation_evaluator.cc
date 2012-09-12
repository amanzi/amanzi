/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  An elevation evaluator getting values from the volumetric mesh.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "composite_vector_function_factory.hh"
#include "standalone_elevation_evaluator.hh"

namespace Amanzi {
namespace Flow {
namespace FlowRelations {

StandaloneElevationEvaluator::StandaloneElevationEvaluator(Teuchos::ParameterList& elev_plist) :
    ElevationEvaluator(),
    elev_plist_(elev_plist) {};

StandaloneElevationEvaluator::StandaloneElevationEvaluator(const StandaloneElevationEvaluator& other) :
    ElevationEvaluator(),
    elev_plist_(other.elev_plist_) {};

Teuchos::RCP<FieldEvaluator> StandaloneElevationEvaluator::Clone() const {
  return Teuchos::rcp(new StandaloneElevationEvaluator(*this));
}

void StandaloneElevationEvaluator::EvaluateElevationAndSlope_(const Teuchos::Ptr<State>& S,
        const std::vector<Teuchos::Ptr<CompositeVector> >& results) {

  Teuchos::Ptr<CompositeVector> elev = results[0];
  Teuchos::Ptr<CompositeVector> slope = results[1];

  // If necessary, create the functions from paramater lists.
  if (elevation_function_ == Teuchos::null) {
    Teuchos::ParameterList elev_plist = elev_plist_.sublist("elevation function");
    elevation_function_ = Functions::CreateCompositeVectorFunction(elev_plist, *elev);
  }

  if (slope_function_ == Teuchos::null) {
    Teuchos::ParameterList slope_plist = elev_plist_.sublist("slope function");
    slope_function_ = Functions::CreateCompositeVectorFunction(slope_plist, *slope);
  }

  // Evaluate the functions.
  elevation_function_->Compute(S->time(), elev);
  slope_function_->Compute(S->time(), slope);
};


} //namespace
} //namespace
} //namespace
