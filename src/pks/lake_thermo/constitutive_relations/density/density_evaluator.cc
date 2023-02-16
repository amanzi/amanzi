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
  double a0 = 800.969e-7;
  double a1 = 588.194e-7;
  double a2 = 811.465e-8;
  double a3 = 476.600e-10;

  double temp0 = 3.85;
  double const_ampl = 1.9549e-5;
  double const_power = 1.68;

  for (CompositeVector::name_iterator comp=result->begin();
       comp!=result->end(); ++comp) {
    const Epetra_MultiVector& temp_v = *temp->ViewComponent(*comp,false);
    Epetra_MultiVector& result_v = *result->ViewComponent(*comp,false);

    int ncomp = result->size(*comp, false);
    std::vector<double> rho(ncomp);
    for (int i=0; i!=ncomp; ++i) {
      double T = temp_v[0][i]-273.15;
      if (T > 0.) { //water
        rho[i] = rho0*(1+a0+a1*T-a2*T*T+a3*T*T*T); // UNESCO
        // result_v[0][i] = rho0*(1.-const_ampl*std::pow(std::fabs(T-temp0),const_power)); // Hostetler
      }
      else { // ice 
        rho[i] = 917.;
      }
      // rho[i] = rho0*(1+a0+a1*T-a2*T*T+a3*T*T*T); 
    }

    for (int i=0; i!=ncomp; ++i) {
      result_v[0][i] = rho[i];
    }

    // result_v[0][0] = rho[0];
    // for (int i=1; i!=ncomp; ++i) {
    //   result_v[0][i] = 0.5*(rho[i]+rho[i-1]);
    // }

  }
};


void DensityEvaluator::EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
        Key wrt_key, const Teuchos::Ptr<CompositeVector>& result) {
  
  result->PutScalar(0.);

  if (wrt_key == temperature_key_) {
    Teuchos::RCP<const CompositeVector> temp = S->GetFieldData(temperature_key_);

    double rho0 = 1000.;
    double a0 = 800.969e-7;
    double a1 = 588.194e-7;
    double a2 = 811.465e-8;
    double a3 = 476.600e-10;

    double temp0 = 3.85;
    double const_ampl = 1.9549e-5;
    double const_power = 1.68;

    for (CompositeVector::name_iterator comp=result->begin();
         comp!=result->end(); ++comp) {
      const Epetra_MultiVector& temp_v = *temp->ViewComponent(*comp,false);
      Epetra_MultiVector& result_v = *result->ViewComponent(*comp,false);

      int ncomp = result->size(*comp, false);
      std::vector<double> drho(ncomp);
      for (int i=0; i!=ncomp; ++i) {
        double T = temp_v[0][i]-273.15;
        if (T > 0.) { //water
          drho[i] = rho0*(a1 - 2.*a2*T + 3.*a3*T*T);
          // double x = T-temp0;
          // double sign = (x > 0) ? 1 : ((x < 0) ? -1 : 0);
          // result_v[0][i] = rho0*(-const_power*const_ampl*std::pow(std::fabs(T-temp0),const_power-1)*sign);
        }
        else { // ice 
          drho[i] = 0.;
        }
        // drho[i] = rho0*(a1 - 2.*a2*T + 3.*a3*T*T);
      }

      for (int i=0; i!=ncomp; ++i) {
        result_v[0][i] = drho[i];
      }

      // result_v[0][0] = drho[0];
      // for (int i=1; i!=ncomp; ++i) {
      //   result_v[0][i] = 0.5*(drho[i]+drho[i-1]);
      // }

    }

  } // if
};


} //namespace
} //namespace
