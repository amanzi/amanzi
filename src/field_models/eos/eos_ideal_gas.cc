/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  ATS

  EOS for an ideal gas

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "eos_ideal_gas.hh"

namespace Amanzi {
namespace Flow {
namespace FlowRelations {

// registry of method
Utils::RegisteredFactoryWithState<EOS,EOSIdealGas> EOSIdealGas::factory_("ideal gas");

EOSIdealGas::EOSIdealGas(Teuchos::ParameterList& eos_plist, const Teuchos::Ptr<State>& S) :
    EOS(eos_plist, S) {
  InitializeFromPlist_();
};

EOSIdealGas::EOSIdealGas(const EOSIdealGas& other) :
    EOS(other),
    R_(other.R_),
    M_(other.M_) {}


// ---------------------------------------------------------------------------
// Virtual copy constructor.
// ---------------------------------------------------------------------------
Teuchos::RCP<FieldModel> EOSIdealGas::Clone() const {
  return Teuchos::rcp(new EOSIdealGas(*this));
}


double EOSIdealGas::Density(double T, double p) {
  return p / (R_*T);
};

double EOSIdealGas::DDensityDT(double T, double p) {
  return -p / (R_*T*T);
};

double EOSIdealGas::DDensityDp(double T, double p) {
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
} // namespace
