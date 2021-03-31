/*
  Chemistry 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Class for Freundlich isotherm
*/

#include <cmath>
#include <iostream>
#include <iomanip>
#include <vector>

#include "sorption_isotherm.hh"
#include "sorption_isotherm_freundlich.hh"

namespace Amanzi {
namespace AmanziChemistry {

SorptionIsothermFreundlich::SorptionIsothermFreundlich()
    : SorptionIsotherm("freundlich", SorptionIsotherm::FREUNDLICH),
      KD_(0.), 
      n_(0.),
      params_(2, 0.0) {
}


SorptionIsothermFreundlich::SorptionIsothermFreundlich(double KD, double n)
    : SorptionIsotherm("freundlich", SorptionIsotherm::FREUNDLICH),
      KD_(KD), 
      n_(n),
      params_(2, 0.0) {
}


double SorptionIsothermFreundlich::Evaluate(const Species& primarySpecies) {
  // Csorb = KD * activity^(n)
  // Units: The units don't make a whole lot of sense.
  return KD_ * std::pow(primarySpecies.activity(), n_);
}


const std::vector<double>& SorptionIsothermFreundlich::GetParameters(void) {
  params_.at(0) = KD_;
  params_.at(1) = n_;
  return params_;
}


void SorptionIsothermFreundlich::SetParameters(const std::vector<double>& params) {
  set_KD(params.at(0));
  set_n(params.at(1));
}


double SorptionIsothermFreundlich::EvaluateDerivative(const Species& primarySpecies) {
  // Csorb = KD * activity^(n)
  // dCsorb/dCaq = KD * n * activity^(n-1) * activity_coef
  //             = Csorb * n / molality
  // Units:
  //  KD [kg water/m^3 bulk]
  double C_sorb = KD_ * std::pow(primarySpecies.activity(), n_);
  return C_sorb * n_ / primarySpecies.molality();
}


void SorptionIsothermFreundlich::Display() const {
  std::cout << std::setw(5) << "KD:"
            << std::scientific << std::setprecision(5)
            << std::setw(15) << KD_
            << std::setw(5) << "n:"
            << std::scientific << std::setprecision(5)
            << std::setw(15) << n_
            << std::endl;
}

}  // namespace AmanziChemistry
}  // namespace Amanzi
