/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  ATS

  EOS for an ideal gas (does not implement viscosity at this point!)

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "eos_vapor_in_gas.hh"

namespace Amanzi {
namespace Flow {
namespace FlowRelations {

double EOSVaporInGas::MassDensity(double T, double p, double mol_frac_vapor) {
  return MassDensity(T,p) * (mol_frac_vapor*Mv_ + (1.-mol_frac_vapor)*Mg_);
};

double EOSVaporInGas::DMassDensityDT(double T, double p, double mol_frac_vapor) {
  return DMolarDensityDT(T,p) * (mol_frac_vapor*Mv_ + (1.-mol_frac_vapor)*Mg_);
};

double EOSVaporInGas::DMassDensityDp(double T, double p, double mol_frac_vapor) {
  return DMolarDensityDp(T,p) * (mol_frac_vapor*Mv_ + (1.-mol_frac_vapor)*Mg_);
};

double EOSVaporInGas::DMassDensityDmol_frac(double T, double p, double mol_frac_vapor) {
  return MolarDensity(T,p)*(Mv_ - Mg_);
};

void EOSVaporInGas::InitializeFromPlist_() {
  R_ = eos_plist_.get<double>("Ideal gas constant [J/mol-K]", 8.3144621);

  if (eos_plist_.isParameter("Molar mass of gas [kg/mol]")) {
    Mg_ = eos_plist_.get<double>("Molar mass of gas [kg/mol]");
  } else {
    Mg_ = eos_plist_.get<double>("Molar mass of gas [g/mol]", 28.956)*1e-3;
  }

  if (eos_plist_.isParameter("Molar mass of vapor [kg/mol]")) {
    Mv_ = eos_plist_.get<double>("Molar mass of vapor [kg/mol]");
  } else {
    Mv_ = eos_plist_.get<double>("Molar mass of vapor [g/mol]", 18.0153)*1e-3;
  }
};

} // namespace
} // namespace
} // namespace
