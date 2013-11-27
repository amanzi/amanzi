/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  Evaluates the conductivity of surface flow.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "unfrozen_fraction_relperm_model.hh"
#include "unfrozen_fraction_relperm_evaluator.hh"

namespace Amanzi {
namespace Flow {
namespace FlowRelations {

UnfrozenFractionRelPermEvaluator::UnfrozenFractionRelPermEvaluator(Teuchos::ParameterList& plist) :
    SecondaryVariableFieldEvaluator(plist) {

  uf_key_ = plist_.get<std::string>("unfrozen fraction key", "unfrozen_fraction");
  dependencies_.insert(uf_key_);

  h_key_ = plist_.get<std::string>("ponded depth key", "ponded_depth");
  dependencies_.insert(h_key_);
  
  if (my_key_ == std::string("")) {
    my_key_ = "unfrozen_fraction_relperm";
  }

  // create the model, hard-coded until we have a 2nd model
  //  ASSERT(plist_.isSublist("unfrozen fraction rel perm model"));
  Teuchos::ParameterList sublist = plist_.sublist("unfrozen fraction rel perm model");
  model_ = Teuchos::rcp(new UnfrozenFractionRelPermModel(sublist));
}


UnfrozenFractionRelPermEvaluator::UnfrozenFractionRelPermEvaluator(const UnfrozenFractionRelPermEvaluator& other) :
    SecondaryVariableFieldEvaluator(other),
    uf_key_(other.uf_key_),
    h_key_(other.h_key_),
    model_(other.model_) {}


Teuchos::RCP<FieldEvaluator>
UnfrozenFractionRelPermEvaluator::Clone() const {
  return Teuchos::rcp(new UnfrozenFractionRelPermEvaluator(*this));
}


// Required methods from SecondaryVariableFieldEvaluator
void UnfrozenFractionRelPermEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result) {
  Teuchos::RCP<const CompositeVector> uf = S->GetFieldData(uf_key_);
  Teuchos::RCP<const CompositeVector> h = S->GetFieldData(h_key_);

  for (CompositeVector::name_iterator comp=result->begin();
       comp!=result->end(); ++comp) {
    const Epetra_MultiVector& uf_v = *uf->ViewComponent(*comp,false);
    const Epetra_MultiVector& h_v = *h->ViewComponent(*comp,false);
    Epetra_MultiVector& result_v = *result->ViewComponent(*comp,false);

    int ncomp = result->size(*comp, false);
    for (int i=0; i!=ncomp; ++i) {
      result_v[0][i] = model_->UnfrozenFractionRelPerm(uf_v[0][i], h_v[0][i]);
    }
  }
}


void UnfrozenFractionRelPermEvaluator::EvaluateFieldPartialDerivative_(
    const Teuchos::Ptr<State>& S,
    Key wrt_key, const Teuchos::Ptr<CompositeVector>& result) {

  ASSERT(wrt_key == uf_key_);
  Teuchos::RCP<const CompositeVector> uf = S->GetFieldData(uf_key_);
  Teuchos::RCP<const CompositeVector> h = S->GetFieldData(h_key_);

  for (CompositeVector::name_iterator comp=result->begin();
       comp!=result->end(); ++comp) {
    const Epetra_MultiVector& uf_v = *uf->ViewComponent(*comp,false);
    const Epetra_MultiVector& h_v = *h->ViewComponent(*comp,false);
    Epetra_MultiVector& result_v = *result->ViewComponent(*comp,false);

    int ncomp = result->size(*comp, false);
    for (int i=0; i!=ncomp; ++i) {
      result_v[0][i] = model_->DUnfrozenFractionRelPermDUnfrozenFraction(uf_v[0][i], h_v[0][i]);
    }
  }
}


} //namespace
} //namespace
} //namespace

