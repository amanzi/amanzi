/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#include "sorption_isotherm_linear.hh"

#include <iostream>
#include <iomanip>

namespace amanzi {
namespace chemistry {

SorptionIsothermLinear::SorptionIsothermLinear()
    : SorptionIsotherm("linear"),
      KD_(0.) {
}  // end SorptionIsothermLinear() constructor

SorptionIsothermLinear::SorptionIsothermLinear(const double KD)
    : SorptionIsotherm("linear"),
      KD_(KD) {
}  // end SorptionIsothermLinear() constructor

SorptionIsothermLinear::~SorptionIsothermLinear() {
}  // end SorptionIsothermLinear() destructor

void SorptionIsothermLinear::Init(const double KD) {
  set_KD(KD);
}

std::vector<double> SorptionIsothermLinear::GetParameters(void) const {
  std::vector<double> params;
  params.resize(1, KD());
  return params;
}  // end GetParameters()

void SorptionIsothermLinear::SetParameters(const std::vector<double>& params) {
  set_KD(params.at(0));
}  // end SetParameters()


double SorptionIsothermLinear::Evaluate(const Species& primarySpecies ) {
  // Csorb = KD * activity
  // Units:
  // sorbed_concentration [mol/m^3 bulk] = KD [kg water/m^3 bulk] * 
  //   activity [mol/kg water]
  return KD() * primarySpecies.activity();
}  // end Evaluate()

double SorptionIsothermLinear::EvaluateDerivative(const Species& primarySpecies) {
  // Csorb = KD * activity
  // dCsorb/dCaq = KD * activity_coef
  // Units:
  //  KD [kg water/m^3 bulk]
  return KD() * primarySpecies.act_coef();
}  // end EvaluateDerivative()

void SorptionIsothermLinear::Display(void) const {
  std::cout << std::setw(5) << "KD:"
            << std::scientific << std::setprecision(5)
            << std::setw(15) << KD() << std::endl;
}  // end Display()

}  // namespace chemistry
}  // namespace amanzi
