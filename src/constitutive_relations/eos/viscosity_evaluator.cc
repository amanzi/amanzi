/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  ViscosityEvaluator is the interface between state/data and the model, a VPM.

  License: BSD
  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "viscosity_relation_factory.hh"
#include "viscosity_evaluator.hh"

namespace Amanzi {
namespace Relations {

// registry of method
Utils::RegisteredFactory<FieldEvaluator,ViscosityEvaluator> ViscosityEvaluator::factory_("viscosity");

ViscosityEvaluator::ViscosityEvaluator(Teuchos::ParameterList& plist) :
    SecondaryVariableFieldEvaluator(plist) {

  // my keys
  if (my_key_ == std::string("")) {
    my_key_ = plist_.get<std::string>("viscosity key", "viscosity_liquid");
  }
  setLinePrefix(my_key_+std::string(" evaluator"));

  // Set up my dependencies.
  std::size_t end = my_key_.find_first_of("_");
  std::string domain_name = my_key_.substr(0,end);
  if (domain_name == std::string("viscosity")) {
    domain_name = std::string("");
  } else {
    domain_name = domain_name+std::string("_");
  }

  // -- temperature
  temp_key_ = plist_.get<std::string>("temperature key",
          domain_name+std::string("temperature"));
  dependencies_.insert(temp_key_);

  // Construct my Viscosity model
  ASSERT(plist_.isSublist("viscosity model parameters"));
  ViscosityRelationFactory visc_fac;
  visc_ = visc_fac.createViscosity(plist_.sublist("viscosity model parameters"));
};


ViscosityEvaluator::ViscosityEvaluator(const ViscosityEvaluator& other) :
    SecondaryVariableFieldEvaluator(other),
    visc_(other.visc_),
    temp_key_(other.temp_key_) {}


Teuchos::RCP<FieldEvaluator> ViscosityEvaluator::Clone() const {
  return Teuchos::rcp(new ViscosityEvaluator(*this));
}


void ViscosityEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
                         const Teuchos::Ptr<CompositeVector>& result) {
  // Pull dependencies out of state.
  Teuchos::RCP<const CompositeVector> temp = S->GetFieldData(temp_key_);

  // evaluate p_s / p_atm
  for (CompositeVector::name_iterator comp=result->begin();
       comp!=result->end(); ++comp) {
    const Epetra_MultiVector& temp_v = *(temp->ViewComponent(*comp,false));
    Epetra_MultiVector& result_v = *(result->ViewComponent(*comp,false));

    int count = result->size(*comp);
    for (int id=0; id!=count; ++id) {
      result_v[0][id] = visc_->Viscosity(temp_v[0][id]);
    }
  }
}


void ViscosityEvaluator::EvaluateFieldPartialDerivative_(
    const Teuchos::Ptr<State>& S, Key wrt_key,
    const Teuchos::Ptr<CompositeVector>& result) {
  ASSERT(wrt_key == temp_key_);

  // Pull dependencies out of state.
  Teuchos::RCP<const CompositeVector> temp = S->GetFieldData(temp_key_);

  // evaluate d/dT( p_s / p_atm )
  for (CompositeVector::name_iterator comp=result->begin();
       comp!=result->end(); ++comp) {
    const Epetra_MultiVector& temp_v = *(temp->ViewComponent(*comp,false));
    Epetra_MultiVector& result_v = *(result->ViewComponent(*comp,false));

    int count = result->size(*comp);
    for (int id=0; id!=count; ++id) {
      result_v[0][id] = visc_->DViscosityDT(temp_v[0][id]);
    }
  }
}

} // namespace
} // namespace
