/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  ATS

  EOS for an ideal gas (does not implement viscosity at this point!)

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "vapor_pressure_model_factory.hh"
#include "eos_factory.hh"
#include "eos_vapor_in_gas.hh"


namespace Amanzi {
namespace Flow {
namespace FlowRelations {
  EOSVaporInGas::EOSVaporInGas(Teuchos::ParameterList& eos_plist) : eos_plist_(eos_plist) {
  InitializeFromPlist_();
}

double EOSVaporInGas::MassDensity(double T, double p, double mol_frac_vapor) {
  return MolarDensity(T,p) * molar_mass(mol_frac_vapor);
};

double EOSVaporInGas::DMassDensityDT(double T, double p, double mol_frac_vapor) {
  return DMolarDensityDT(T,p) * molar_mass(mol_frac_vapor);
};

double EOSVaporInGas::DMassDensityDp(double T, double p, double mol_frac_vapor) {
  return DMolarDensityDp(T,p) * molar_mass(mol_frac_vapor);
};

double EOSVaporInGas::DMassDensityDmol_frac(double T, double p, double mol_frac_vapor) {
  return MolarDensity(T,p)*(Mv_ - gas_eos_->molar_mass());
};

double EOSVaporInGas::SaturatedVaporPressure(double T) {
  return sat_vapor_model_->SaturatedVaporPressure(T);
}

void EOSVaporInGas::InitializeFromPlist_() {
  Teuchos::ParameterList gas_eos_plist = eos_plist_.sublist("Gas EOS");
  EOSFactory eos_factory;
  gas_eos_ = eos_factory.createEOS(gas_eos_plist);

  Teuchos::ParameterList vappres_plist = eos_plist_.sublist("Vapor Pressure Model");
  VaporPressureModelFactory vappres_factory;
  sat_vapor_model_ = vappres_factory.createVaporPressureModel(vappres_plist);

  // eos variables
  if (eos_plist_.isParameter("Molar mass of vapor [kg/mol]")) {
    Mv_ = eos_plist_.get<double>("Molar mass of vapor [kg/mol]");
  } else {
    Mv_ = eos_plist_.get<double>("Molar mass of vapor [g/mol]", 18.0153)*1e-3;
  }


};

} // namespace
} // namespace
} // namespace
