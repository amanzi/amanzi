/*
  Chemistry 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Class for linear isotherm
*/

#include "sorption_isotherm_linear.hh"

#include <iostream>
#include <iomanip>
#include <string>

#include "VerboseObject.hh"
#include "sorption_isotherm.hh"

namespace Amanzi {
namespace AmanziChemistry {

SorptionIsothermLinear::SorptionIsothermLinear()
    : SorptionIsotherm("linear", SorptionIsotherm::LINEAR),
      KD_(0.0),
      params_(1, 0.0) {
}


SorptionIsothermLinear::SorptionIsothermLinear(const double KD)
    : SorptionIsotherm("linear", SorptionIsotherm::LINEAR),
      KD_(KD),
      params_(1, 0.0) {
}


void SorptionIsothermLinear::Init(const double KD) {
  set_KD(KD);
}


const std::vector<double>& SorptionIsothermLinear::GetParameters(void) {
  params_.at(0) = KD();
  return params_;
}


void SorptionIsothermLinear::SetParameters(const std::vector<double>& params) {
  set_KD(params.at(0));
}


double SorptionIsothermLinear::Evaluate(const Species& primarySpecies ) {
  // Csorb = KD * activity
  // Units:
  // sorbed_concentration [mol/m^3 bulk] = KD [kg water/m^3 bulk] * 
  //   activity [mol/kg water]
  return KD() * primarySpecies.activity();
}


double SorptionIsothermLinear::EvaluateDerivative(const Species& primarySpecies) {
  // Csorb = KD * activity
  // dCsorb/dCaq = KD * activity_coef
  // Units:
  //  KD [kg water/m^3 bulk]
  return KD() * primarySpecies.act_coef();
}


void SorptionIsothermLinear::Display(const Teuchos::RCP<VerboseObject>& vo) const {
  std::stringstream message;
  message << std::setw(5) << "KD:"
          << std::scientific << std::setprecision(5)
          << std::setw(15) << KD() << std::endl;
  vo->Write(Teuchos::VERB_HIGH, message.str());
}

}  // namespace AmanziChemistry
}  // namespace Amanzi
