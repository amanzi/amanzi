/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/* -----------------------------------------------------------------------------
ATS

Authors: Svetlana Tokareva (tokareva@lanl.gov)

FieldEvaluator for water density.
----------------------------------------------------------------------------- */


#include "soil_energy_evaluator.hh"

namespace Amanzi {
namespace SoilThermo {

SoilEnergyEvaluator::SoilEnergyEvaluator(Teuchos::ParameterList& plist) :
    SecondaryVariableFieldEvaluator(plist) {
  if (my_key_.empty()) {

    my_key_ = plist_.get<std::string>("energy key", "surface-energy");
  }

  // Set up my dependencies.
  std::string domain_name = Keys::getDomain(my_key_);

  // -- temperature
  temperature_key_ = Keys::readKey(plist_, domain_name, "temperature", "temperature");
  dependencies_.insert(temperature_key_);

  // -- density
  density_key_ = Keys::readKey(plist_, domain_name, "density", "density");
  dependencies_.insert(density_key_);

  // -- heat capacity
  heat_capacity_key_ = Keys::readKey(plist_, domain_name, "heat capacity", "heat_capacity");
  dependencies_.insert(heat_capacity_key_);

  // -- pressure
  pres_key_ = Keys::readKey(plist_, domain_name, "pressure", "pressure");
  dependencies_.insert(pres_key_);

};

SoilEnergyEvaluator::SoilEnergyEvaluator(const SoilEnergyEvaluator& other) :
    SecondaryVariableFieldEvaluator(other),
    temperature_key_(other.temperature_key_),
    density_key_(other.density_key_),
    heat_capacity_key_(other.heat_capacity_key_),
    pres_key_(other.pres_key_){};

Teuchos::RCP<FieldEvaluator>
SoilEnergyEvaluator::Clone() const {
  return Teuchos::rcp(new SoilEnergyEvaluator(*this));
};


void SoilEnergyEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result) {
  Teuchos::RCP<const CompositeVector> temp = S->GetFieldData(temperature_key_);

  double rho0 = 1200.;
  double cp0 = 800./rho0;

  S->GetFieldEvaluator(density_key_)->HasFieldChanged(S.ptr(), "soil thermo");
  // evaluate density
  const Epetra_MultiVector& rho =
  *S->GetFieldData(density_key_)->ViewComponent("cell",false);

  S->GetFieldEvaluator(heat_capacity_key_)->HasFieldChanged(S.ptr(), "soil thermo");

  // evaluate heat capacity
  const Epetra_MultiVector& cp =
  *S->GetFieldData(heat_capacity_key_)->ViewComponent("cell",false);

  for (CompositeVector::name_iterator comp=result->begin();
       comp!=result->end(); ++comp) {
    const Epetra_MultiVector& temp_v = *temp->ViewComponent(*comp,false);
    Epetra_MultiVector& result_v = *result->ViewComponent(*comp,false);

    int ncomp = result->size(*comp, false);
    for (int i=0; i!=ncomp; ++i) {
      double T = temp_v[0][i];
      result_v[0][i] = rho[0][i]*cp[0][i]*T;
//      result_v[0][i] = rho0*cp0*T;
    }
  }
};


void SoilEnergyEvaluator::EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
        Key wrt_key, const Teuchos::Ptr<CompositeVector>& result) {

  result->PutScalar(0.);

  if (wrt_key == temperature_key_) {
    Teuchos::RCP<const CompositeVector> temp = S->GetFieldData(temperature_key_);

    double rho0 = 1200.;
    double cp0 = 800./rho0;

    // evaluate density
    const Epetra_MultiVector& rho =
    *S->GetFieldData(density_key_)->ViewComponent("cell",false);

    // evaluate heat capacity
    const Epetra_MultiVector& cp =
    *S->GetFieldData(heat_capacity_key_)->ViewComponent("cell",false);

    for (CompositeVector::name_iterator comp=result->begin();
         comp!=result->end(); ++comp) {
      const Epetra_MultiVector& temp_v = *temp->ViewComponent(*comp,false);
      Epetra_MultiVector& result_v = *result->ViewComponent(*comp,false);

      int ncomp = result->size(*comp, false);
      for (int i=0; i!=ncomp; ++i) {
        double T = temp_v[0][i];
        result_v[0][i] = rho[0][i]*cp[0][i];
//        result_v[0][i] = rho0*cp0;
      }
    }
  }

  if (wrt_key == pres_key_) {
      result->PutScalar(0.);
  }

};


} //namespace
} //namespace
