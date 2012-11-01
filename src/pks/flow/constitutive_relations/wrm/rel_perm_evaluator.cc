/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  Rel perm( pc ( sat ) ).

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "rel_perm_evaluator.hh"
#include "wrm_factory.hh"

namespace Amanzi {
namespace Flow {
namespace FlowRelations {

// RelPermEvaluator::RelPermEvaluator(Teuchos::ParameterList& plist) :
//     SecondaryVariableFieldEvaluator(),
//     plist_(plist) {
//   ASSERT(plist.isSublist("WRM parameters"));
//   Teuchos::ParameterList sublist = plist.sublist("WRM parameters");
//   WRMFactory fac;
//   wrm_ = fac.createWRM(sublist);

//   InitializeFromPlist_();
// }

RelPermEvaluator::RelPermEvaluator(Teuchos::ParameterList& plist,
        const Teuchos::RCP<WRMRegionPairList>& wrms) :
    SecondaryVariableFieldEvaluator(plist),
    wrms_(wrms) {
  InitializeFromPlist_();
}

RelPermEvaluator::RelPermEvaluator(const RelPermEvaluator& other) :
    SecondaryVariableFieldEvaluator(other),
    wrms_(other.wrms_),
    sat_key_(other.sat_key_) {}


Teuchos::RCP<FieldEvaluator>
RelPermEvaluator::Clone() const {
  return Teuchos::rcp(new RelPermEvaluator(*this));
}


void RelPermEvaluator::InitializeFromPlist_() {
  // my keys are for saturation and rel perm.
  my_key_ = plist_.get<string>("rel perm key", "relative_permeability");
  setLinePrefix(my_key_+std::string(" evaluator"));

  // my dependencies are just saturation.
  sat_key_ = plist_.get<string>("saturation key", "saturation_liquid");
  dependencies_.insert(sat_key_);
}


void RelPermEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result) {

  Teuchos::RCP<const CompositeVector> sat = S->GetFieldData(sat_key_);

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
      (*result)("cell", *c) = region->second->k_relative(pc);
    }
  }
}


void RelPermEvaluator::EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
        Key wrt_key, const Teuchos::Ptr<CompositeVector>& results) {
  // not implemented, not yet needed --etc
  ASSERT(0);
}



} //namespace
} //namespace
} //namespace
