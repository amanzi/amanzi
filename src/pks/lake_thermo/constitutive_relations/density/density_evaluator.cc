/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/* -----------------------------------------------------------------------------
ATS

Authors: Svetlana Tokareva (tokareva@lanl.gov)

FieldEvaluator for water density.
----------------------------------------------------------------------------- */


#include "density_evaluator.hh"

namespace Amanzi {
namespace LakeThermo {

DensityEvaluator::DensityEvaluator(Teuchos::ParameterList& plist) :
    SecondaryVariableFieldEvaluator(plist) {
  if (my_key_.empty()) {

    my_key_ = plist_.get<std::string>("density key", "surface-density");
  }

  // Set up my dependencies.
  std::string domain_name = Keys::getDomain(my_key_);

  // -- temperature
  temperature_key_ = Keys::readKey(plist_, domain_name, "temperature", "temperature");
  dependencies_.insert(temperature_key_);

};

DensityEvaluator::DensityEvaluator(const DensityEvaluator& other) :
    SecondaryVariableFieldEvaluator(other),
    temperature_key_(other.temperature_key_) {};

Teuchos::RCP<FieldEvaluator>
DensityEvaluator::Clone() const {
  return Teuchos::rcp(new DensityEvaluator(*this));
};


void DensityEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result) {

  Teuchos::RCP<const CompositeVector> temp = S->GetFieldData(temperature_key_);

//  double rho0 = 1.;
  double rho0 = 1000.;

  for (CompositeVector::name_iterator comp=result->begin();
       comp!=result->end(); ++comp) {
    const Epetra_MultiVector& temp_v = *temp->ViewComponent(*comp,false);
    Epetra_MultiVector& result_v = *result->ViewComponent(*comp,false);

    int ncomp = result->size(*comp, false);
    for (int i=0; i!=ncomp; ++i) {
      double T = temp_v[0][i]-273.15;
      result_v[0][i] = rho0 * (1.+8.0*1.e-5 + 5.88*1.e-5*T - 8.11*1.e-6*T*T + 4.77*1.e-8*T*T*T);
    }
  }
};


void DensityEvaluator::EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
        Key wrt_key, const Teuchos::Ptr<CompositeVector>& result) {
  // not implemented
  if (wrt_key == temperature_key_) {
    Teuchos::RCP<const CompositeVector> temp = S->GetFieldData(temperature_key_);

    double rho0 = 1000.;

    for (CompositeVector::name_iterator comp=result->begin();
         comp!=result->end(); ++comp) {
      const Epetra_MultiVector& temp_v = *temp->ViewComponent(*comp,false);
      Epetra_MultiVector& result_v = *result->ViewComponent(*comp,false);

      int ncomp = result->size(*comp, false);
      for (int i=0; i!=ncomp; ++i) {
        double T = temp_v[0][i]-273.15;
        if (T > 0.) { //water}
          result_v[0][i] = rho0 * (5.88*1.e-5 - 2.*8.11*1.e-6*T + 3.*4.77*1.e-8*T*T);
        }
        else { // ice 
          result_v[0][i] = 917.;
        }
      }
    }
  }
};


} //namespace
} //namespace
