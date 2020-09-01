/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  Interface for a thermal conductivity of lake model.

  License: BSD
  Authors: Svetlana Tokareva (tokareva@lanl.gov)
*/

#include "thermal_conductivity_evaluator.hh"

namespace Amanzi {
namespace LakeThermo {

ThermalConductivityEvaluator::ThermalConductivityEvaluator(
      Teuchos::ParameterList& plist) :
    SecondaryVariableFieldEvaluator(plist) {
  if (my_key_ == std::string("")) {
    my_key_ = plist_.get<std::string>("thermal conductivity key",
            "surface-thermal_conductivity");
  }

  Key domain = Keys::getDomain(my_key_);

  // Set up my dependencies.
  std::string domain_name = Keys::getDomain(my_key_);

  // -- temperature
  temperature_key_ = Keys::readKey(plist_, domain_name, "temperature", "temperature");
  dependencies_.insert(temperature_key_);

  AMANZI_ASSERT(plist_.isSublist("thermal conductivity parameters"));
  Teuchos::ParameterList sublist = plist_.sublist("thermal conductivity parameters");
//  K_liq_ = sublist.get<double>("thermal conductivity of water [W/(m-K)]", 0.58);
//  K_ice_ = sublist.get<double>("thermal conductivity of ice [W/(m-K)]", 2.18);
//  min_K_ = sublist.get<double>("minimum thermal conductivity", 1.e-14);

  // later: read these parameters from xml
  K_max_    = 150; // [W/(m * K)]
  K_0_      = 50.;
  V_wind_   = 10.; // [m/s]
  V_wind_0_ = 20.;
}


ThermalConductivityEvaluator::ThermalConductivityEvaluator(
      const ThermalConductivityEvaluator& other) :
    SecondaryVariableFieldEvaluator(other),
    K_max_(other.K_max_),
    K_0_(other.K_0_),
    V_wind_(other.V_wind_),
    V_wind_0_(other.V_wind_0_) {}
//    uf_key_(other.uf_key_),
//    height_key_(other.height_key_),
//    K_liq_(other.K_liq_),
//    K_ice_(other.K_ice_),
//    min_K_(other.min_K_) {}


Teuchos::RCP<FieldEvaluator>
ThermalConductivityEvaluator::Clone() const {
  return Teuchos::rcp(new ThermalConductivityEvaluator(*this));
}

void ThermalConductivityEvaluator::EvaluateField_(
      const Teuchos::Ptr<State>& S,
      const Teuchos::Ptr<CompositeVector>& result) {

  std::cout << "In ThermalConductivityEvaluator::EvaluateField_" << std::endl;

  for (CompositeVector::name_iterator comp=result->begin();
         comp!=result->end(); ++comp) {
      // much more efficient to pull out vectors first
//      const Epetra_MultiVector& eta_v = *eta->ViewComponent(*comp,false);
//      const Epetra_MultiVector& height_v = *height->ViewComponent(*comp,false);
      Epetra_MultiVector& result_v = *result->ViewComponent(*comp,false);

      int ncomp = result->size(*comp, false);
      std::cout << "ncomp = " << ncomp << std::endl;
      for (int i=0; i!=ncomp; ++i) {
        if (ice_cover_) {
            result_v[0][i] = 1.5;
        } else {
            result_v[0][i] = 10.*K_0_ + V_wind_/V_wind_0_*(K_max_ - K_0_);
        }
        std::cout << "conductivity: " << result_v[0][i] << std::endl;
      }
    }

}


void ThermalConductivityEvaluator::EvaluateFieldPartialDerivative_(
      const Teuchos::Ptr<State>& S, Key wrt_key,
      const Teuchos::Ptr<CompositeVector>& result) {
  std::cout<<"THERMAL CONDUCITIVITY: Derivative not implemented yet!"<<wrt_key<<"\n";
  AMANZI_ASSERT(0); // not implemented, not yet needed
  result->Scale(1.e-6); // convert to MJ
}

} //namespace
} //namespace
