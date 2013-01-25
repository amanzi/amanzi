/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  Unfrozen ponded depth

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "unfrozen_ponded_depth_evaluator.hh"

namespace Amanzi {
namespace Flow {
namespace FlowRelations {

UnfrozenPondedDepthEvaluator::UnfrozenPondedDepthEvaluator(Teuchos::ParameterList& plist) :
    SecondaryVariableFieldEvaluator(plist) {
  my_key_ = plist_.get<std::string>("unfrozen ponded depth key",
          "unfrozen_ponded_depth");
  setLinePrefix(my_key_+std::string(" evaluator"));

  depth_key_ = plist_.get<std::string>("ponded depth key", "ponded_depth");
  dependencies_.insert(depth_key_);

  temp_key_ = plist_.get<std::string>("surface temperature key",
          "surface_temperature");
  dependencies_.insert(temp_key_);


}


UnfrozenPondedDepthEvaluator::UnfrozenPondedDepthEvaluator(const UnfrozenPondedDepthEvaluator& other) :
    SecondaryVariableFieldEvaluator(other),
    depth_key_(other.depth_key_),
    temp_key_(other.temp_key_) {}


Teuchos::RCP<FieldEvaluator>
UnfrozenPondedDepthEvaluator::Clone() const {
  return Teuchos::rcp(new UnfrozenPondedDepthEvaluator(*this));
}



void UnfrozenPondedDepthEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result) {

  Teuchos::RCP<const CompositeVector> depth = S->GetFieldData(depth_key_);
  Teuchos::RCP<const CompositeVector> temp = S->GetFieldData(temp_key_);

  const Epetra_MultiVector& res_v = *result->ViewComponent("cell",false);

  // Evaluate the evaluator to calculate sat.
  for (WRMRegionPairList::iterator region=wrms_->begin();
       region!=wrms_->end(); ++region) {
    std::string name = region->first;
    int ncells = sat->mesh()->get_set_size(name,
            AmanziMesh::CELL, AmanziMesh::OWNED);
    AmanziMesh::Entity_ID_List cells(ncells);
    sat->mesh()->get_set_entities(name,
            AmanziMesh::CELL, AmanziMesh::OWNED, &cells);

    // use the wrm to evaluate saturation on each cell in the region
    for (AmanziMesh::Entity_ID_List::iterator c=cells.begin(); c!=cells.end(); ++c) {
      double pc = region->second->capillaryPressure((*sat)("cell", *c));
      res_v[0][*c] = std::max(region->second->k_relative(pc), min_val_);
    }
  }
}


void UnfrozenPondedDepthEvaluator::EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
        Key wrt_key, const Teuchos::Ptr<CompositeVector>& results) {
  // not implemented, not yet needed --etc
  ASSERT(0);
}



} //namespace
} //namespace
} //namespace
