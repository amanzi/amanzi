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

RelPermEvaluator::RelPermEvaluator(Teuchos::ParameterList& wrm_plist) :
    SecondaryVariableFieldEvaluator(),
    wrm_plist_(wrm_plist) {
  ASSERT(wrm_plist.isSublist("WRM parameters"));
  Teuchos::ParameterList sublist = wrm_plist.sublist("WRM parameters");
  WRMFactory fac;
  wrm_ = fac.createWRM(sublist);

  InitializeFromPlist_();
}

RelPermEvaluator::RelPermEvaluator(Teuchos::ParameterList& wrm_plist, const Teuchos::RCP<WRM>& wrm) :
    SecondaryVariableFieldEvaluator(),
    wrm_plist_(wrm_plist),
    wrm_(wrm) {
  InitializeFromPlist_();
}

RelPermEvaluator::RelPermEvaluator(const RelPermEvaluator& other) :
    SecondaryVariableFieldEvaluator(other),
    wrm_plist_(other.wrm_plist_),
    wrm_(other.wrm_),
    sat_key_(other.sat_key_) {}


Teuchos::RCP<FieldEvaluator>
RelPermEvaluator::Clone() const {
  return Teuchos::rcp(new RelPermEvaluator(*this));
}


void RelPermEvaluator::InitializeFromPlist_() {
  // my keys are for saturation and rel perm.
  my_key_ = wrm_plist_.get<string>("rel perm key", "relative_permeability");

  // my dependencies are just saturation.
  sat_key_ = wrm_plist_.get<string>("saturation key", "saturation_liquid");
  dependencies_.insert(sat_key_);
}


void RelPermEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result) {

  Teuchos::RCP<const CompositeVector> sat = S->GetFieldData(sat_key_);

  // Loop over names in the target and then owned entities in that name,
  // evaluating the evaluator to calculate sat.
  for (CompositeVector::name_iterator comp=sat->begin();
       comp!=sat->end(); ++comp) {
    for (int id=0; id!=sat->size(*comp); ++id) {
      (*result)(*comp, id) = wrm_->k_relative(wrm_->capillaryPressure((*sat)(*comp, id)));
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
