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


EOSVaporInGas::EOSVaporInGas(Teuchos::ParameterList& eos_plist) : eos_plist_(eos_plist) {
  InitializeFromPlist_();
}

double EOSVaporInGas::MolarDensity(double T, double p) {
  return eos_->MolarDensity(T,p);
};

double EOSVaporInGas::DMolarDensityDT(double T, double p, double mol_frac_vapor) {
  return eos_->DMolarDensityDT(T,p);
};

double EOSVaporInGas::DMolarDensityDp(double T, double p, double mol_frac_vapor) {
  return eos_->DMolarDensityDp(T,p);
};



void EOSVaporInGas::InitializeFromPlist_() {
  Teuchos::ParameterList gas_eos_plist = eos_plist_.sublist("Gas EOS");
  EOSFactory eos_factory;
  gas_eos_ = eos_factory.createEOS(gas_eos_plist);
};

} // namespace
} // namespace
