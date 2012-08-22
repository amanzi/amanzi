/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  ATS

  EOS for an ideal gas (does not implement viscosity at this point!)

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "eos_factory.hh"
#include "eos_vapor_in_gas.hh"

namespace Amanzi {
namespace Relations {

Utils::RegisteredFactory<EOS,EOSVaporInGas> EOSVaporInGas::factory_("vapor in gas");

EOSVaporInGas::EOSVaporInGas(Teuchos::ParameterList& eos_plist) : eos_plist_(eos_plist) {
  InitializeFromPlist_();
}

double EOSVaporInGas::MolarDensity(double T, double p) {
  return gas_eos_->MolarDensity(T,p);
};

double EOSVaporInGas::DMolarDensityDT(double T, double p) {
  return gas_eos_->DMolarDensityDT(T,p);
};

double EOSVaporInGas::DMolarDensityDp(double T, double p) {
  return gas_eos_->DMolarDensityDp(T,p);
};



void EOSVaporInGas::InitializeFromPlist_() {
  Teuchos::ParameterList gas_eos_plist = eos_plist_.sublist("gas EOS parameters");
  EOSFactory eos_factory;
  gas_eos_ = eos_factory.createEOS(gas_eos_plist);
};

} // namespace
} // namespace
