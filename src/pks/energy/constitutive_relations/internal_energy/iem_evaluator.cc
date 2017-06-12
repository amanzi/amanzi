/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  The WRM Evaluator simply calls the WRM with the correct arguments.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "iem_evaluator.hh"
#include "iem_factory.hh"

namespace Amanzi {
namespace Energy {
namespace EnergyRelations {


IEMEvaluator::IEMEvaluator(Teuchos::ParameterList& plist) :
    SecondaryVariableFieldEvaluator(plist) {

  ASSERT(plist_.isSublist("IEM parameters"));
  Teuchos::ParameterList sublist = plist_.sublist("IEM parameters");
  IEMFactory fac;
  iem_ = fac.createIEM(sublist);

  InitializeFromPlist_();
}


IEMEvaluator::IEMEvaluator(Teuchos::ParameterList& plist, const Teuchos::RCP<IEM>& iem) :
    SecondaryVariableFieldEvaluator(plist),
    iem_(iem) {

  InitializeFromPlist_();
}


IEMEvaluator::IEMEvaluator(const IEMEvaluator& other) :
    SecondaryVariableFieldEvaluator(other),
    iem_(other.iem_),
    temp_key_(other.temp_key_) {}


Teuchos::RCP<FieldEvaluator>
IEMEvaluator::Clone() const {
  return Teuchos::rcp(new IEMEvaluator(*this));
}


void IEMEvaluator::InitializeFromPlist_() {
  if (my_key_.empty()) {
    my_key_ = plist_.get<std::string>("internal energy key");
  }

  // Set up my dependencies.
  std::string domain_name = getDomain(my_key_);

  // -- temperature
  temp_key_ = plist_.get<std::string>("temperature key",
          getKey(domain_name, "temperature"));
  dependencies_.insert(temp_key_);
}


void IEMEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result) {
  Teuchos::RCP<const CompositeVector> temp = S->GetFieldData(temp_key_);

  for (CompositeVector::name_iterator comp=result->begin();
       comp!=result->end(); ++comp) {
    const Epetra_MultiVector& temp_v = *temp->ViewComponent(*comp,false);
    Epetra_MultiVector& result_v = *result->ViewComponent(*comp,false);

    int ncomp = result->size(*comp, false);
    for (int i=0; i!=ncomp; ++i) {
      result_v[0][i] = iem_->InternalEnergy(temp_v[0][i]);
    }
  }
}


void IEMEvaluator::EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
        Key wrt_key, const Teuchos::Ptr<CompositeVector>& result) {
  ASSERT(wrt_key == temp_key_);
  Teuchos::RCP<const CompositeVector> temp = S->GetFieldData(temp_key_);

  for (CompositeVector::name_iterator comp=result->begin();
       comp!=result->end(); ++comp) {
    const Epetra_MultiVector& temp_v = *temp->ViewComponent(*comp,false);
    Epetra_MultiVector& result_v = *result->ViewComponent(*comp,false);

    int ncomp = result->size(*comp, false);
    for (int i=0; i!=ncomp; ++i) {
      result_v[0][i] = iem_->DInternalEnergyDT(temp_v[0][i]);
    }
  }
}


} //namespace
} //namespace
} //namespace
