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

ViscosityEvaluator::ViscosityEvaluator(
    Teuchos::ParameterList& visc_plist) :
    visc_plist_(visc_plist) {

  // my keys
  my_key_ = visc_plist_.get<std::string>("molar fraction of water key");

  // Set up my dependencies.
  std::size_t end = my_key_.find_first_of("_");
  std::string domain_name = my_key_.substr(0,end);
  if (domain_name == std::string("viscosity")) {
    domain_name = std::string("");
  } else {
    domain_name = domain_name+std::string("_");
  }

  // -- temperature
  temp_key_ = visc_plist_.get<std::string>("temperature key",
          domain_name+std::string("temperature"));
  dependencies_.insert(temp_key_);

  // Construct my Viscosity model
  ASSERT(visc_plist_.isSublist("viscosity model parameters"));
  ViscosityRelationFactory visc_fac;
  visc_ = visc_fac.createViscosity(visc_plist_.sublist("viscosity model parameters"));
};


ViscosityEvaluator::ViscosityEvaluator(const ViscosityEvaluator& other) :
    SecondaryVariableFieldEvaluator(other),
    visc_plist_(other.visc_plist_),
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
    for (int id=0; id!=result->size(*comp); ++id) {
      (*result)(*comp, id) =
          visc_->Viscosity((*temp)(*comp, id));
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
    for (int id=0; id!=result->size(*comp); ++id) {
      (*result)(*comp, id) =
          visc_->DViscosityDT((*temp)(*comp, id));
    }
  }
}

} // namespace
} // namespace
