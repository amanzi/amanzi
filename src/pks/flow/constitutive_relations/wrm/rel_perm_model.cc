/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  The WRM Model simply calls the WRM with the correct arguments.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "rel_perm_model.hh"
#include "wrm_factory.hh"

namespace Amanzi {
namespace Flow {
namespace FlowRelations {

RelPermModel::RelPermModel(Teuchos::ParameterList& wrm_plist) :
    SecondaryVariableFieldModel(),
    wrm_plist_(wrm_plist) {
  ASSERT(wrm_plist.isSublist("WRM parameters"));
  Teuchos::ParameterList sublist = wrm_plist.sublist("WRM parameters");
  WRMFactory fac;
  wrm_ = fac.createWRM(sublist);

  InitializeFromPlist_();
}

RelPermModel::RelPermModel(Teuchos::ParameterList& wrm_plist, const Teuchos::RCP<WRM>& wrm) :
    SecondaryVariableFieldModel(),
    wrm_plist_(wrm_plist),
    wrm_(wrm) {
  InitializeFromPlist_();
}

RelPermModel::RelPermModel(const RelPermModel& other) :
    SecondaryVariableFieldModel(other),
    wrm_plist_(other.wrm_plist_),
    wrm_(other.wrm_) {}


Teuchos::RCP<FieldModel>
RelPermModel::Clone() const {
  return Teuchos::rcp(new RelPermModel(*this));
}


void RelPermModel::InitializeFromPlist_() {
  // my keys are for saturation and rel perm.
  my_key_ = wrm_plist_.get<string>("rel perm key", "relative_permeability");

  // my dependencies are just saturation.
  sat_key_ = wrm_plist_.get<string>("saturation key", "saturation_liquid");
  dependencies_.insert(sat_key_);
}


void RelPermModel::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result) {

  Teuchos::RCP<const CompositeVector> sat = S->GetFieldData(sat_key_);

  // Loop over names in the target and then owned entities in that name,
  // evaluating the model to calculate sat.
  for (CompositeVector::name_iterator comp=sat->begin();
       comp!=sat->end(); ++comp) {
    for (int id=0; id!=sat->size(*comp); ++id) {
      (*result)(*comp, id) = wrm_->k_relative(wrm_->capillaryPressure((*sat)(*comp, id)));
    }
  }
}


void RelPermModel::EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
        Key wrt_key, const Teuchos::Ptr<CompositeVector>& results) {
  // not implemented
  ASSERT(0);
}



} //namespace
} //namespace
} //namespace
