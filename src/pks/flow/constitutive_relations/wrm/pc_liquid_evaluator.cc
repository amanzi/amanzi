/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  PCLiquidEvaluator is the interface between state/data and the model, a PC relation.

  License: BSD
  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "pc_liq_atm.hh"
#include "pc_liquid_evaluator.hh"

namespace Amanzi {
namespace Flow {
namespace FlowRelations {

// registry of method
Utils::RegisteredFactory<FieldEvaluator,PCLiquidEvaluator> PCLiquidEvaluator::factory_("capillary pressure, atmospheric gas over liquid");

PCLiquidEvaluator::PCLiquidEvaluator(Teuchos::ParameterList& plist) :
    SecondaryVariableFieldEvaluator(plist) {

  // my keys
  my_key_ = plist_.get<std::string>("gas-liquid capillary pressure key",
          "capillary_pressure_gas_liq");
  setLinePrefix(my_key_+std::string(" evaluator"));

  // -- pressure
  pres_key_ = plist_.get<std::string>("pressure key", "pressure");
  dependencies_.insert(pres_key_);

  p_atm_key_ = plist_.get<std::string>("atmospheric pressure key",
          "atmospheric_pressure");

  // Construct my PCLiquid model
  model_ = Teuchos::rcp(new PCLiqAtm(plist_.sublist("capillary pressure model parameters")));
};


PCLiquidEvaluator::PCLiquidEvaluator(const PCLiquidEvaluator& other) :
    SecondaryVariableFieldEvaluator(other),
    model_(other.model_),
    pres_key_(other.pres_key_),
    p_atm_key_(other.p_atm_key_) {}


Teuchos::RCP<FieldEvaluator> PCLiquidEvaluator::Clone() const {
  return Teuchos::rcp(new PCLiquidEvaluator(*this));
}


void PCLiquidEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
                         const Teuchos::Ptr<CompositeVector>& result) {
  // Pull dependencies out of state.
  Teuchos::RCP<const CompositeVector> pres = S->GetFieldData(pres_key_);
  Teuchos::RCP<const double> p_atm = S->GetScalarData(p_atm_key_);

  // evaluate pc
  for (CompositeVector::name_iterator comp=result->begin();
       comp!=result->end(); ++comp) {
    const Epetra_MultiVector& pres_v = *(pres->ViewComponent(*comp,false));
    Epetra_MultiVector& result_v = *(result->ViewComponent(*comp,false));

    int count = result->size(*comp);
    for (int id=0; id!=count; ++id) {
      result_v[0][id] = model_->CapillaryPressure(pres_v[0][id], *p_atm);
    }
  }
}


void PCLiquidEvaluator::EvaluateFieldPartialDerivative_(
    const Teuchos::Ptr<State>& S, Key wrt_key,
    const Teuchos::Ptr<CompositeVector>& result) {
  ASSERT(wrt_key == pres_key_);

  // Pull dependencies out of state.
  Teuchos::RCP<const CompositeVector> pres = S->GetFieldData(pres_key_);
  Teuchos::RCP<const double> p_atm = S->GetScalarData(p_atm_key_);

  // evaluate d/dT( p_s / p_atm )
  for (CompositeVector::name_iterator comp=result->begin();
       comp!=result->end(); ++comp) {
    const Epetra_MultiVector& pres_v = *(pres->ViewComponent(*comp,false));
    Epetra_MultiVector& result_v = *(result->ViewComponent(*comp,false));

    int count = result->size(*comp);
    for (int id=0; id!=count; ++id) {
      result_v[0][id] = model_->DCapillaryPressureDp(pres_v[0][id], *p_atm);
    }
  }
}

} // namespace
} // namespace
} // namespace
