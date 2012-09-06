/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  ATS

  EOS for an ideal gas (does not implement viscosity at this point!)

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "eos_ideal_gas.hh"

namespace Amanzi {
namespace Relations {

// registry of method
Utils::RegisteredFactory<EOS,EOSIdealGas> EOSIdealGas::factory_("ideal gas");

EOSIdealGas::EOSIdealGas(Teuchos::ParameterList& eos_plist) : eos_plist_(eos_plist) {
  InitializeFromPlist_();
};

double EOSIdealGas::MolarDensity(double T, double p) {
  return p / (R_*T);
};

double EOSIdealGas::DMolarDensityDT(double T, double p) {
  return -p / (R_*T*T);
};

double EOSIdealGas::DMolarDensityDp(double T, double p) {
  return 1.0 / (R_*T);
};


void EOSIdealGas::InitializeFromPlist_() {
  R_ = eos_plist_.get<double>("Ideal gas constant [J/mol-K]", 8.3144621);

  if (eos_plist_.isParameter("Molar mass of gas [kg/mol]")) {
    M_ = eos_plist_.get<double>("Molar mass of gas [kg/mol]");
  } else {
    M_ = eos_plist_.get<double>("Molar mass of gas [g/mol]", 28.956)*1e-3;
  }
};

} // namespace
} // namespace
