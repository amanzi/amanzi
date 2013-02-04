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
      Teuchos::ParameterList& plist) :
    SecondaryVariableFieldEvaluator(plist) {
  if (my_key_ == std::string("")) {
    my_key_ = plist_.get<std::string>("thermal conductivity key", "thermal_conductivity");
  }
  setLinePrefix(my_key_+std::string(" evaluator"));

  poro_key_ = plist_.get<std::string>("porosity key", "porosity");
  dependencies_.insert(poro_key_);

  temp_key_ = plist_.get<std::string>("temperature key", "temperature");
  dependencies_.insert(temp_key_);

  sat_key_ = plist_.get<std::string>("saturation key", "saturation_liquid");
  dependencies_.insert(sat_key_);

  sat2_key_ = plist_.get<std::string>("second saturation key", "saturation_ice");
  dependencies_.insert(sat2_key_);

  ASSERT(plist_.isSublist("thermal conductivity parameters"));
  Teuchos::ParameterList sublist = plist_.sublist("thermal conductivity parameters");
  ThermalConductivityThreePhaseFactory fac;
  tc_ = fac.createThermalConductivityModel(sublist);
}


ThermalConductivityThreePhaseEvaluator::ThermalConductivityThreePhaseEvaluator(
      const ThermalConductivityThreePhaseEvaluator& other) :
    SecondaryVariableFieldEvaluator(other),
    poro_key_(other.poro_key_),
    temp_key_(other.temp_key_),
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
  Teuchos::RCP<const CompositeVector> temp = S->GetFieldData(temp_key_);
  Teuchos::RCP<const CompositeVector> sat = S->GetFieldData(sat_key_);
  Teuchos::RCP<const CompositeVector> sat2 = S->GetFieldData(sat2_key_);

  for (CompositeVector::name_iterator comp=poro->begin();
       comp!=poro->end(); ++comp) {
    // much more efficient to pull out vectors first
    const Epetra_MultiVector& poro_v = *poro->ViewComponent(*comp,false);
    const Epetra_MultiVector& temp_v = *temp->ViewComponent(*comp,false);
    const Epetra_MultiVector& sat_v = *sat->ViewComponent(*comp,false);
    const Epetra_MultiVector& sat2_v = *sat2->ViewComponent(*comp,false);
    Epetra_MultiVector& result_v = *result->ViewComponent(*comp,false);

    int ncomp = result->size(*comp, false);
    for (int i=0; i!=ncomp; ++i) {
      result_v[0][i] = tc_->ThermalConductivity(poro_v[0][i],
              sat_v[0][i], sat2_v[0][i], temp_v[0][i]);
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
