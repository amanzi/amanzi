/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#include "sorption_isotherm_freundlich.hh"

#include <cmath>
#include <vector>
#include <iostream>
#include <iomanip>

#include "sorption_isotherm.hh"

namespace amanzi {
namespace chemistry {

SorptionIsothermFreundlich::SorptionIsothermFreundlich()
    : SorptionIsotherm("freundlich", SorptionIsotherm::FREUNDLICH),
      KD_(0.), 
      n_(0.),
      params_(2, 0.0) {
}  // end SorptionIsothermLangmuir() constructor

SorptionIsothermFreundlich::SorptionIsothermFreundlich(const double KD, 
                                                       const double n)
    : SorptionIsotherm("freundlich", SorptionIsotherm::FREUNDLICH),
      KD_(KD), 
      n_(n),
      params_(2, 0.0) {
}  // end SorptionIsothermLangmuir() constructor

SorptionIsothermFreundlich::~SorptionIsothermFreundlich() {
}  // end SorptionIsothermLangmuir() destructor

double SorptionIsothermFreundlich::Evaluate(const Species& primarySpecies) {
  // Csorb = KD * activity^(n)
  // Units: The units don't make a whole lot of sense.
  return KD() * std::pow(primarySpecies.activity(), n());
}  // end Evaluate()

const std::vector<double>& SorptionIsothermFreundlich::GetParameters(void) {
  params_.at(0) = KD();
  params_.at(1) = n();
  return params_;
}  // end GetParameters()

void SorptionIsothermFreundlich::SetParameters(const std::vector<double>& params) {
  set_KD(params.at(0));
  set_n(params.at(1));
}  // end SetParameters()

double SorptionIsothermFreundlich::EvaluateDerivative(
    const Species& primarySpecies) {
  // Csorb = KD * activity^(n)
  // dCsorb/dCaq = KD * n * activity^(n-1) * activity_coef
  //             = Csorb * n / molality
  // Units:
  //  KD [kg water/m^3 bulk]
  double C_sorb = KD() * std::pow(primarySpecies.activity(), n());
  return C_sorb * n() / primarySpecies.molality();
}  // end EvaluateDerivative()

void SorptionIsothermFreundlich::Display(void) const {
  std::cout << std::setw(5) << "KD:"
            << std::scientific << std::setprecision(5)
            << std::setw(15) << KD() 
            << std::setw(5) << "n:"
            << std::scientific << std::setprecision(5)
            << std::setw(15) << n() 
            << std::endl;
}  // end Display()

}  // namespace chemistry
}  // namespace amanzi
