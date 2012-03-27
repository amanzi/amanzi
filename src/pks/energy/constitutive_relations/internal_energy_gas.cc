/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Internal energy model for gas, relative to water @237.15K

See ATS process model documentation's permafrost physical properties
documentation for details.

UNITS: J/mol
------------------------------------------------------------------------- */

#include "internal_energy_gas.hh"

namespace Amanzi {
namespace Energy {
namespace EnergyRelations {

InternalEnergyGas::InternalEnergyGas(Teuchos::ParameterList& plist) :
    plist_(plist) {
  InitializeFromPlist_();
};

double InternalEnergyGas::InternalEnergy(double temp, double mol_frac_gas) {
  return (1.0 + 0.622*mol_frac_gas) * Cv_air_ * (temp - 273.15) + mol_frac_gas*heat_vaporization_;
};

double InternalEnergyGas::DInternalEnergyDT(double temp, double mol_frac_gas) {
  // evaluated at constant mol_frac gas for now?
  return (1.0 + 0.622*mol_frac_gas) * Cv_air_;
};

void InternalEnergyGas::InitializeFromPlist_() {
  molar_basis_ = plist_.get<bool>("molar-basis (otherwise, mass-basis)", true);
  Cv_air_ = plist_.get<double>("heat capacity of air [J/(mol-K)]", 13.0);
  heat_vaporization_ = plist_.get<double>("heat of vaporization of water [J/mol]", 4.065e4);
};

} // namespace
} // namespace
} // namespace
