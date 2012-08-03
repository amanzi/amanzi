/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  Interface for a thermal conductivity model with three phases.

  License: BSD
  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "dbc.hh"
#include "thermal_conductivity_threephase_factory.hh"
#include "thermal_conductivity_threephase_evaluator.hh"

namespace Amanzi {
namespace Energy {
namespace EnergyRelations {

ThermalConductivityThreePhaseEvaluator::ThermalConductivityThreePhaseEvaluator(
      Teuchos::ParameterList& tc_plist) :
    tc_plist_(tc_plist) {
  my_key_ = tc_plist_.get<std::string>("thermal conductivity key", "thermal_conductivity");

  poro_key_ = tc_plist_.get<std::string>("porosity key", "porosity");
  dependencies_.insert(poro_key_);

  sat_key_ = tc_plist_.get<std::string>("saturation key", "saturation_liquid");
  dependencies_.insert(sat_key_);

  sat2_key_ = tc_plist_.get<std::string>("second saturation key", "saturation_ice");
  dependencies_.insert(sat2_key_);

  ASSERT(tc_plist_.isSublist("thermal conductivity parameters"));
  Teuchos::ParameterList sublist = tc_plist_.sublist("thermal conductivity parameters");
  ThermalConductivityThreePhaseFactory fac;
  tc_ = fac.createThermalConductivityModel(sublist);
}


ThermalConductivityThreePhaseEvaluator::ThermalConductivityThreePhaseEvaluator(
      const ThermalConductivityThreePhaseEvaluator& other) :
    SecondaryVariableFieldEvaluator(other),
    tc_plist_(other.tc_plist_),
    poro_key_(other.poro_key_),
    sat_key_(other.sat_key_),
    sat2_key_(other.sat2_key_),
    tc_(other.tc_) {}

Teuchos::RCP<FieldEvaluator>
ThermalConductivityThreePhaseEvaluator::Clone() const {
  return Teuchos::rcp(new ThermalConductivityThreePhaseEvaluator(*this));
}


void ThermalConductivityThreePhaseEvaluator::EvaluateField_(
      const Teuchos::Ptr<State>& S,
      const Teuchos::Ptr<CompositeVector>& result) {
  // pull out the dependencies
  Teuchos::RCP<const CompositeVector> poro = S->GetFieldData(poro_key_);
  Teuchos::RCP<const CompositeVector> sat = S->GetFieldData(sat_key_);
  Teuchos::RCP<const CompositeVector> sat2 = S->GetFieldData(sat2_key_);

  for (CompositeVector::name_iterator comp=poro->begin();
       comp!=poro->end(); ++comp) {
    for (int i=0; i!=poro->size(*comp); ++i) {
      (*result)(*comp, i) = tc_->ThermalConductivity((*poro)(*comp, i),
              (*sat)(*comp, i), (*sat2)(*comp,i));
    }
  }
}


void ThermalConductivityThreePhaseEvaluator::EvaluateFieldPartialDerivative_(
      const Teuchos::Ptr<State>& S, Key wrt_key,
      const Teuchos::Ptr<CompositeVector>& result) {
  ASSERT(0); // not implemented, not yet needed
}


} //namespace
} //namespace
} //namespace
