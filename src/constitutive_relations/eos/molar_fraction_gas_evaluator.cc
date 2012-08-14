/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  Determining the molar fraction of a gas component within a gas mixture.

  License: BSD
  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "vapor_pressure_relation_factory.hh"
#include "molar_fraction_gas_evaluator.hh"

namespace Amanzi {
namespace Relations {

// registry of method
Utils::RegisteredFactory<FieldEvaluator,MolarFractionGasEvaluator> MolarFractionGasEvaluator::factory_("molar fraction gas");


MolarFractionGasEvaluator::MolarFractionGasEvaluator(Teuchos::ParameterList& mfg_plist) :
    mfg_plist_(mfg_plist) {

  // set up the actual model
  ASSERT(mfg_plist_.isSublist("vapor pressure model parameters"));
  VaporPressureRelationFactory vpm_fac;
  sat_vapor_model_ = vpm_fac.createVaporPressure(
      mfg_plist_.sublist("vapor pressure model parameters"));

  // process the list for my provided field.
  if (mfg_plist_.isParameter("molar fraction key")) {
    my_key_ = mfg_plist_.get<string>("molar fraction key");
  } else {
    std::string name = mfg_plist_.name();
    std::size_t start = name.find_last_of(">");
    my_key_ = name.substr(start+1);
  }

  // set up dependencies
  std::size_t end = my_key_.find_first_of("_");
  std::string domain_name = my_key_.substr(0,end);
  if (domain_name == std::string("mol") ||
      domain_name == std::string("molar")) {
    domain_name = std::string("");
  } else {
    domain_name = domain_name+std::string("_");
  }

  temp_key_ = mfg_plist_.get<std::string>("temperature key",
          domain_name+std::string("temperature"));
  dependencies_.insert(temp_key_);
}


MolarFractionGasEvaluator::MolarFractionGasEvaluator(const MolarFractionGasEvaluator& other) :
    SecondaryVariableFieldEvaluator(other),
    mfg_plist_(other.mfg_plist_),
    sat_vapor_model_(other.sat_vapor_model_),
    temp_key_(other.temp_key_) {}


Teuchos::RCP<FieldEvaluator>
MolarFractionGasEvaluator::Clone() const {
  return Teuchos::rcp(new MolarFractionGasEvaluator(*this));
}


void MolarFractionGasEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result) {
  // Pull dependencies out of state.
  Teuchos::RCP<const CompositeVector> temp = S->GetFieldData(temp_key_);
  Teuchos::RCP<const double> p_atm = S->GetScalarData("atmospheric_pressure");

  // evaluate p_s / p_atm
  for (CompositeVector::name_iterator comp=result->begin();
       comp!=result->end(); ++comp) {
    for (int id=0; id!=result->size(*comp); ++id) {
      double svp = sat_vapor_model_->SaturatedVaporPressure((*temp)(*comp, id));
      (*result)(*comp, id) = svp / (*p_atm);
    }
  }
}


void MolarFractionGasEvaluator::EvaluateFieldPartialDerivative_(
    const Teuchos::Ptr<State>& S, Key wrt_key,
    const Teuchos::Ptr<CompositeVector>& result) {
  ASSERT(wrt_key == temp_key_);

  // Pull dependencies out of state.
  Teuchos::RCP<const CompositeVector> temp = S->GetFieldData(temp_key_);
  Teuchos::RCP<const double> p_atm = S->GetScalarData("atmospheric_pressure");

  // evaluate d/dT( p_s / p_atm )
  for (CompositeVector::name_iterator comp=result->begin();
       comp!=result->end(); ++comp) {
    for (int id=0; id!=result->size(*comp); ++id) {
      double dsvp = sat_vapor_model_->DSaturatedVaporPressureDT((*temp)(*comp, id));
      (*result)(*comp, id) = dsvp / (*p_atm);
    }
  }
}


} // namespace
} // namespace

