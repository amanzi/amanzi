/*
  This is the energy component of the ATS and Amanzi codes. 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon

  Internal energy model for water vapor, relative to water @237.15K
  See ATS process model documentation's permafrost physical properties
  documentation for details.

  UNITS: J/mol
*/

#include "iem_water_vapor.hh"

namespace Amanzi {
namespace Energy {

IEMWaterVapor::IEMWaterVapor(Teuchos::ParameterList& plist) :
    plist_(plist) {
  InitializeFromPlist_();
};


double IEMWaterVapor::InternalEnergy(double temp, double mol_frac_gas) {
  return (1.0 + 0.622*mol_frac_gas) * Cv_air_ * (temp - 273.15) + mol_frac_gas*heat_vaporization_;
};


double IEMWaterVapor::DInternalEnergyDT(double temp, double mol_frac_gas) {
  // evaluated at constant mol_frac gas for now?
  return (1.0 + 0.622*mol_frac_gas) * Cv_air_;
};


double IEMWaterVapor::DInternalEnergyDomega(double temp, double mol_frac_gas) {
  // evaluated at constant mol_frac gas for now?
  return heat_vaporization_ + 0.622 * Cv_air_ * (temp - 273.15);
};


void IEMWaterVapor::InitializeFromPlist_() {
  molar_basis_ = plist_.get<bool>("molar-basis (otherwise, mass-basis)", true);
  Cv_air_ = plist_.get<double>("heat capacity of air [J/(mol-K)]", 13.0);
  heat_vaporization_ = plist_.get<double>("heat of vaporization of water [J/mol]", 4.065e4);
};

}  // namespace Energy
}  // namespace Amanzi

