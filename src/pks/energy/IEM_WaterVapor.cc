/*
  Energy

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon

  Internal energy model for water vapor, relative to water @237.15K
  UNITS: [J/mol]
*/

#include "IEM_WaterVapor.hh"

namespace Amanzi {
namespace Energy {

IEM_WaterVapor::IEM_WaterVapor(Teuchos::ParameterList& plist) :
    plist_(plist) {
  InitializeFromPlist_();
};


double IEM_WaterVapor::InternalEnergy(double temp, double mol_frac_gas) {
  return (1.0 + 0.622 * mol_frac_gas) * Cv_air_ * (temp - 273.15) + mol_frac_gas*heat_vaporization_;
};


double IEM_WaterVapor::DInternalEnergyDT(double temp, double mol_frac_gas) {
  // evaluated at constant mol_frac gas for now?
  return (1.0 + 0.622 * mol_frac_gas) * Cv_air_;
};


double IEM_WaterVapor::DInternalEnergyDomega(double temp, double mol_frac_gas) {
  // evaluated at constant mol_frac gas for now?
  return heat_vaporization_ + 0.622 * Cv_air_ * (temp - 273.15);
};


void IEM_WaterVapor::InitializeFromPlist_() {
  Cv_air_ = plist_.get<double>("molar heat capacity of air", 13.0);
  heat_vaporization_ = plist_.get<double>("heat of vaporization of water [J/mol]", 4.065e4);
};

}  // namespace Energy
}  // namespace Amanzi

