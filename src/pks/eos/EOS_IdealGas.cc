/*
  EOS
   
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@lanl.gov)

  EOS for an ideal gas. It does not implement viscosity at this point!
*/

#include "EOS_IdealGas.hh"

namespace Amanzi {
namespace EOS {

EOS_IdealGas::EOS_IdealGas(Teuchos::ParameterList& eos_plist) : eos_plist_(eos_plist) {
  InitializeFromPlist_();
};


double EOS_IdealGas::MolarDensity(double T, double p) {
  return p / (R_ * T);
};


double EOS_IdealGas::DMolarDensityDT(double T, double p) {
  return -p / (R_ * T * T);
};


double EOS_IdealGas::DMolarDensityDp(double T, double p) {
  return 1.0 / (R_ * T);
};


void EOS_IdealGas::InitializeFromPlist_()
{
  R_ = eos_plist_.get<double>("ideal gas constant [J/mol-K]", 8.3144621);

  if (eos_plist_.isParameter("molar mass of gas [kg/mol]")) {
    M_ = eos_plist_.get<double>("molar mass of gas [kg/mol]");
  } else {
    M_ = eos_plist_.get<double>("molar mass of gas [g/mol]", 28.956)*1e-3;
  }
};
 
}  // namespace EOS
}  // namespace Amanzi
