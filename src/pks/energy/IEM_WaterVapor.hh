/*
  Energy

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon

  Internal energy model for water_vapor, relative to water @237.15K
  UNITS: [J/mol]
*/

#ifndef AMANZI_ENERGY_IEM_WATER_VAPOR_HH_
#define AMANZI_ENERGY_IEM_WATER_VAPOR_HH_

#include "Teuchos_ParameterList.hpp"

namespace Amanzi {
namespace Energy {

class IEM_WaterVapor {
 public:
  IEM_WaterVapor(Teuchos::ParameterList& plist);

  double InternalEnergy(double temp, double mol_frac_gas);
  double DInternalEnergyDT(double temp, double mol_frac_gas);
  double DInternalEnergyDomega(double temp, double mol_frac_gas);

 private:
  void InitializeFromPlist_();

  Teuchos::ParameterList plist_;

  double Cv_air_;            // units: J/(mol-K)
  double heat_vaporization_; // units: J/mol
};

} // namespace Energy
} // namespace Amanzi

#endif
