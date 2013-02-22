/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  Evaluates the conductivity of surface flow.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "independent_variable_field_evaluator.hh"
#include "overland_conductivity_evaluator.hh"
#include "manning_conductivity_model.hh"

namespace Amanzi {
namespace Flow {
namespace FlowRelations {

OverlandConductivityEvaluator::OverlandConductivityEvaluator(Teuchos::ParameterList& plist) :
    SecondaryVariableFieldEvaluator(plist) {

  depth_key_ = plist_.get<std::string>("height key", "ponded_depth");
  dependencies_.insert(depth_key_);

  slope_key_ = plist_.get<std::string>("slope key", "slope_magnitude");
  dependencies_.insert(slope_key_);

  coef_key_ = plist_.get<std::string>("coefficient key", "manning_coefficient");
  dependencies_.insert(coef_key_);

  my_key_ = "overland_conductivity";
  setLinePrefix(my_key_+std::string(" evaluator"));

  // create the model, hard-coded until we have a 2nd model
  ASSERT(plist_.isSublist("overland conductivity model"));
  Teuchos::ParameterList sublist = plist_.sublist("overland conductivity model");
  model_ = Teuchos::rcp(new ManningConductivityModel(sublist));
}


OverlandConductivityEvaluator::OverlandConductivityEvaluator(const OverlandConductivityEvaluator& other) :
    SecondaryVariableFieldEvaluator(other),
    depth_key_(other.depth_key_),
    slope_key_(other.slope_key_),
    coef_key_(other.coef_key_),
    model_(other.model_) {}


Teuchos::RCP<FieldEvaluator>
OverlandConductivityEvaluator::Clone() const {
  return Teuchos::rcp(new OverlandConductivityEvaluator(*this));
}


// Required methods from SecondaryVariableFieldEvaluator
void OverlandConductivityEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result) {

  Teuchos::RCP<const CompositeVector> depth = S->GetFieldData(depth_key_);
  Teuchos::RCP<const CompositeVector> slope = S->GetFieldData(slope_key_);
  Teuchos::RCP<const CompositeVector> coef = S->GetFieldData(coef_key_);

  for (CompositeVector::name_iterator comp=result->begin();
       comp!=result->end(); ++comp) {
    const Epetra_MultiVector& depth_v = *depth->ViewComponent(*comp,false);
    const Epetra_MultiVector& slope_v = *slope->ViewComponent(*comp,false);
    const Epetra_MultiVector& coef_v = *coef->ViewComponent(*comp,false);
    Epetra_MultiVector& result_v = *result->ViewComponent(*comp,false);

    int ncomp = result->size(*comp, false);
    for (int i=0; i!=ncomp; ++i) {
      result_v[0][i] = model_->Conductivity(depth_v[0][i], slope_v[0][i], coef_v[0][i]);
    }
  }
}


void OverlandConductivityEvaluator::EvaluateFieldPartialDerivative_(
    const Teuchos::Ptr<State>& S,
    Key wrt_key, const Teuchos::Ptr<CompositeVector>& result) {

  ASSERT(0);
  // not implemented, likely not needed.
}


} //namespace
} //namespace
} //namespace

