/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  Interface for a heat flux at the surface of soil model.

  License: BSD
  Authors: Svetlana Tokareva (tokareva@lanl.gov)
 */

#include "soil_heat_flux_bc_evaluator.hh"

namespace Amanzi {
namespace SoilThermo {

SoilHeatFluxBCEvaluator::SoilHeatFluxBCEvaluator(
    Teuchos::ParameterList& plist) :
        SecondaryVariableFieldEvaluator(plist) {
  if (my_key_ == std::string("")) {
    my_key_ = plist_.get<std::string>("soil heat flux bc key",
        "surface-heat_flux_bc");
  }

  Key domain = Keys::getDomain(my_key_);

  // Set up my dependencies.
  std::string domain_name = Keys::getDomain(my_key_);

  // -- temperature
  temperature_key_ = Keys::readKey(plist_, domain_name, "temperature", "temperature");
  dependencies_.insert(temperature_key_);

  AMANZI_ASSERT(plist_.isSublist("soil heat flux bc parameters"));
  Teuchos::ParameterList sublist = plist_.sublist("soil heat flux bc parameters");

  // later: read these parameters from xml
  SS = 0;
  alpha = 0;
  E_a = 0;
  E_s = 0;
  H = 0;
  LE = 0;

}


SoilHeatFluxBCEvaluator::SoilHeatFluxBCEvaluator(
    const SoilHeatFluxBCEvaluator& other) :
        SecondaryVariableFieldEvaluator(other),
        SS(other.SS),
        alpha(other.alpha),
        E_a(other.E_a),
        E_s(other.E_s),
        H(other.H),
        LE(other.LE),
        temperature_key_(other.temperature_key_){}


Teuchos::RCP<FieldEvaluator>
SoilHeatFluxBCEvaluator::Clone() const {
  return Teuchos::rcp(new SoilHeatFluxBCEvaluator(*this));
}

void SoilHeatFluxBCEvaluator::EvaluateField_(
    const Teuchos::Ptr<State>& S,
    const Teuchos::Ptr<CompositeVector>& result) {

  ice_cover_ = false; // first always assume that there is no ice

  // get temperature
  Teuchos::RCP<const CompositeVector> temp = S->GetFieldData(temperature_key_);

  for (CompositeVector::name_iterator comp=result->begin();
      comp!=result->end(); ++comp) {

    const Epetra_MultiVector& temp_v = *temp->ViewComponent(*comp,false);
    Epetra_MultiVector& result_v = *result->ViewComponent(*comp,false);

    int ncomp = result->size(*comp, false);

    for (int i=0; i!=ncomp; ++i) {
      result_v[0][i] = SS*(1.-alpha) + E_a - E_s - H - LE;
    } // i

  }

}

void SoilHeatFluxBCEvaluator::EvaluateFieldPartialDerivative_(
    const Teuchos::Ptr<State>& S, Key wrt_key,
    const Teuchos::Ptr<CompositeVector>& result) {
  std::cout<<"HEAT FLUX BC: Derivative not implemented yet!"<<wrt_key<<"\n";
  AMANZI_ASSERT(0); // not implemented, not yet needed
  result->Scale(1.e-6); // convert to MJ
}

} //namespace
} //namespace
