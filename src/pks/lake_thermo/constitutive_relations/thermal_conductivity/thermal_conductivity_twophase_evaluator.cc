/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  Interface for a thermal conductivity model with two phases.

  License: BSD
  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "dbc.hh"
#include "thermal_conductivity_twophase_factory.hh"
#include "thermal_conductivity_twophase_evaluator.hh"

namespace Amanzi {
namespace Energy {

ThermalConductivityTwoPhaseEvaluator::ThermalConductivityTwoPhaseEvaluator(
      Teuchos::ParameterList& plist) :
    SecondaryVariableFieldEvaluator(plist) {
  if (my_key_ == std::string("")) {
    my_key_ = plist_.get<std::string>("thermal conductivity key", "thermal_conductivity");
  }

  poro_key_ = plist_.get<std::string>("porosity key", "porosity");
  dependencies_.insert(poro_key_);

  sat_key_ = plist_.get<std::string>("saturation key", "saturation_liquid");
  dependencies_.insert(sat_key_);

  AMANZI_ASSERT(plist_.isSublist("thermal conductivity parameters"));
  Teuchos::ParameterList sublist = plist_.sublist("thermal conductivity parameters");
  ThermalConductivityTwoPhaseFactory fac;
  tc_ = fac.createThermalConductivityModel(sublist);
}


ThermalConductivityTwoPhaseEvaluator::ThermalConductivityTwoPhaseEvaluator(
      const ThermalConductivityTwoPhaseEvaluator& other) :
    SecondaryVariableFieldEvaluator(other),
    poro_key_(other.poro_key_),
    sat_key_(other.sat_key_),
    tc_(other.tc_) {}

Teuchos::RCP<FieldEvaluator>
ThermalConductivityTwoPhaseEvaluator::Clone() const {
  return Teuchos::rcp(new ThermalConductivityTwoPhaseEvaluator(*this));
}


void ThermalConductivityTwoPhaseEvaluator::EvaluateField_(
      const Teuchos::Ptr<State>& S,
      const Teuchos::Ptr<CompositeVector>& result) {
  // pull out the dependencies
  Teuchos::RCP<const CompositeVector> poro = S->GetFieldData(poro_key_);
  Teuchos::RCP<const CompositeVector> sat = S->GetFieldData(sat_key_);

  for (CompositeVector::name_iterator comp=result->begin();
       comp!=result->end(); ++comp) {
    // much more efficient to pull out vectors first
    const Epetra_MultiVector& poro_v = *poro->ViewComponent(*comp,false);
    const Epetra_MultiVector& sat_v = *sat->ViewComponent(*comp,false);
    Epetra_MultiVector& result_v = *result->ViewComponent(*comp,false);

    int ncomp = result->size(*comp, false);
    for (int i=0; i!=ncomp; ++i) {
      result_v[0][i] = tc_->ThermalConductivity(poro_v[0][i], sat_v[0][i]);
    }
  }
  result->Scale(1.e-6); // convert to MJ
}


void ThermalConductivityTwoPhaseEvaluator::EvaluateFieldPartialDerivative_(
      const Teuchos::Ptr<State>& S, Key wrt_key,
      const Teuchos::Ptr<CompositeVector>& result) {
  AMANZI_ASSERT(0); // not implemented, not yet needed
  result->Scale(1.e-6); // convert to MJ
}

} //namespace
} //namespace
