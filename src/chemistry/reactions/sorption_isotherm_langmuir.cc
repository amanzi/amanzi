/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#include "sorption_isotherm_langmuir.hh"

#include <vector>
#include <iostream>
#include <iomanip>

#include "sorption_isotherm.hh"

namespace amanzi {
namespace chemistry {

SorptionIsothermLangmuir::SorptionIsothermLangmuir()
    : SorptionIsotherm("langmuir", SorptionIsotherm::LANGMUIR),
      K_(0.), 
      b_(0.),
      params_(2, 0.0) {
}  // end SorptionIsothermLangmuir() constructor

SorptionIsothermLangmuir::SorptionIsothermLangmuir(const double K, 
                                                   const double b)
    : SorptionIsotherm("langmuir", SorptionIsotherm::LANGMUIR),
      K_(K), 
      b_(b),
      params_(2, 0.0) {
}  // end SorptionIsothermLangmuir() constructor

SorptionIsothermLangmuir::~SorptionIsothermLangmuir() {
}  // end SorptionIsothermLangmuir() destructor

void SorptionIsothermLangmuir::Init(const double K, const double b) {
  set_K(K);
  set_b(b);
}

const std::vector<double>& SorptionIsothermLangmuir::GetParameters(void) {
  params_.at(0) = K();
  params_.at(1) = b();
  return params_;
}  // end GetParameters()

void SorptionIsothermLangmuir::SetParameters(const std::vector<double>& params) {
  set_K(params.at(0));
  set_b(params.at(1));
}  // end SetParameters()


double SorptionIsothermLangmuir::Evaluate(const Species& primarySpecies) {
  // Csorb = K * activity * b / (1 + K * activity)
  // Units:
  // sorbed_concentration [mol/m^3 bulk] = 
  //   K [kg water/mol] * activity [mol/kg water] * b [mol/m^3 bulk] /
  //     (1. + K [kg water/mol] * activity [mol/kg water])
  //
  // NOTE(bandre): need to be careful with the variable names
  // here. Looking at Langmuir (1997), would lead one to expect:
  // Csorb = K * b * activity / (1 + b * activity)
  double K_activity = K() * primarySpecies.activity(); // temporary variable
  return K_activity * b() / (1. + K_activity);
}  // end Evaluate()

double SorptionIsothermLangmuir::EvaluateDerivative(
    const Species& primarySpecies) {
  // Csorb = K * activity * b / (1 + K * activity)
  // dCsorb/dCaq = (K * activity_coef * b / (1 + K * activity)) - 
  //               (K * activity * b / (1 + K * activity)^2 * K * activity_coef)
  // Units:
  //  KD [kg water/m^3 bulk]
  double K_activity = K() * primarySpecies.activity(); // temporary variable
  double C_sorb = K_activity * b() / (1. + K_activity);
  return C_sorb / primarySpecies.molality() - 
           (C_sorb / (1. + K_activity) * K_activity / 
             primarySpecies.molality());
}  // end EvaluateDerivative()

void SorptionIsothermLangmuir::Display(void) const {
  std::cout << std::setw(5) << "K:"
            << std::scientific << std::setprecision(5)
            << std::setw(15) << K() 
            << std::setw(5) << "b:"
            << std::scientific << std::setprecision(5)
            << std::setw(15) << b() 
            << std::endl;
}  // end Display()

}  // namespace chemistry
}  // namespace amanzi
