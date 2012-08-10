/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  Determining the molar fraction of a gas component within a gas mixture.

  License: BSD
  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "vapor_pressure_model_factory.hh"
#include "molar_fraction_gas.hh"

namespace Amanzi {
namespace Flow {
namespace FlowRelations {

// registry of method
Utils::RegisteredFactory<FieldModel,MolarFractionGas> MolarFractionGas::factory_("molar fraction gas");


MolarFractionGas::MolarFractionGas(Teuchos::ParameterList& mfg_plist) :
    mfg_plist_(mfg_plist) {

  // set up the actual model
  VaporPressureModelFactory fac;
  sat_vapor_model_ = fac.createVaporPressureModel(mfg_plist);

  // process the list for my providied field.
  if (mfg_plist_.isParameter("molar fraction key")) {
    my_key_ = mfg_plist_.get<string>("molar fraction key");
  } else {
    std::string name = mfg_plist_.name();
    std::size_t start = name.find_last_of(">");
    my_key_ = name.substr(start+1);
  }

  // set up dependencies
  temp_key_ = mfg_plist_.get<string>("temperature key", "temperature");
  dependencies_.insert(temp_key_);
}


MolarFractionGas::MolarFractionGas(const MolarFractionGas& other) :
    SecondaryVariableFieldModel(other),
    mfg_plist_(other.mfg_plist_),
    sat_vapor_model_(other.sat_vapor_model_),
    temp_key_(other.temp_key_) {}


Teuchos::RCP<FieldModel>
MolarFractionGas::Clone() const {
  return Teuchos::rcp(new MolarFractionGas(*this));
}


void MolarFractionGas::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result) {

  Teuchos::RCP<const CompositeVector> temp = S->GetFieldData(temp_key_);
  Teuchos::RCP<const double> p_atm = S->GetScalarData("atmospheric_pressure");

  // Loop over names in the target and then owned entities in that name,
  // evaluating the model.
  for (CompositeVector::name_iterator comp=result->begin();
       comp!=result->end(); ++comp) {
    for (int id=0; id!=result->size(*comp); ++id) {
      double svp = sat_vapor_model_->SaturatedVaporPressure((*temp)(*comp, id));
      (*result)(*comp, id) = svp / (*p_atm);
    }
  }
}


void MolarFractionGas::EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
        Key wrt_key, const Teuchos::Ptr<CompositeVector>& result) {
  // note we do NOT evaluate dKrel/d wrt_key.  Not sure if this breaks things or not...
  ASSERT(wrt_key == temp_key_);

  Teuchos::RCP<const CompositeVector> temp = S->GetFieldData(temp_key_);
  Teuchos::RCP<const double> p_atm = S->GetScalarData("atmospheric_pressure");

  // Loop over names in the target and then owned entities in that name,
  // evaluating the model.
  for (CompositeVector::name_iterator comp=result->begin();
       comp!=result->end(); ++comp) {
    for (int id=0; id!=result->size(*comp); ++id) {
      double svp = sat_vapor_model_->DSaturatedVaporPressureDT((*temp)(*comp, id));
      (*result)(*comp, id) = svp / (*p_atm);
    }
  }
}


} // namespace
} // namespace
} // namespace

