/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

/*
  Chemistry

  Class for Langmuir isotherm
*/

#include <vector>
#include <iostream>
#include <iomanip>

#include "SorptionIsotherm.hh"
#include "SorptionIsothermLangmuir.hh"

namespace Amanzi {
namespace AmanziChemistry {

SorptionIsothermLangmuir::SorptionIsothermLangmuir()
  : SorptionIsotherm("langmuir", SorptionIsotherm::LANGMUIR), K_(0.0), b_(0.0), params_(2, 0.0)
{}


SorptionIsothermLangmuir::SorptionIsothermLangmuir(double K, double b)
  : SorptionIsotherm("langmuir", SorptionIsotherm::LANGMUIR), K_(K), b_(b), params_(2, 0.0)
{}


void
SorptionIsothermLangmuir::Init(double K, double b)
{
  K_ = K;
  b_ = b;
}


const std::vector<double>&
SorptionIsothermLangmuir::GetParameters()
{
  params_.at(0) = K_;
  params_.at(1) = b_;
  return params_;
}


void
SorptionIsothermLangmuir::SetParameters(const std::vector<double>& params)
{
  K_ = params.at(0);
  b_ = params.at(1);
}


/* *******************************************************************
* Csorb = K * activity * b / (1 + K * activity)
* Units:
* sorbed_concentration [mol/m^3 bulk] =
*   K [kg water/mol] * activity [mol/kg water] * b [mol/m^3 bulk] /
*     (1. + K [kg water/mol] * activity [mol/kg water])
*
* NOTE(bandre): need to be careful with the variable names
* here. Looking at Langmuir (1997), would lead one to expect:
* Csorb = K * b * activity / (1 + b * activity)
******************************************************************* */
double
SorptionIsothermLangmuir::Evaluate(const Species& primary_species)
{
  double K_activity = K_ * primary_species.activity();
  return K_activity * b_ / (1.0 + K_activity);
}


/* *******************************************************************
* Csorb = K * activity * b / (1 + K * activity)
* dCsorb/dCaq = (K * activity_coef * b / (1 + K * activity)) -
*               (K * activity * b / (1 + K * activity)^2 * K * activity_coef)
* Units:
*   KD [kg water/m^3 bulk]
******************************************************************* */
double
SorptionIsothermLangmuir::EvaluateDerivative(const Species& primary_species)
{
  double K_activity = K_ * primary_species.activity();
  double C_sorb = K_activity * b_ / (1.0 + K_activity);
  return C_sorb / primary_species.molality() -
         (C_sorb / (1.0 + K_activity) * K_activity / primary_species.molality());
}

} // namespace AmanziChemistry
} // namespace Amanzi
