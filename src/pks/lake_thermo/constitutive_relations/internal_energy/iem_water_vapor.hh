/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Internal energy model for water_vapor, relative to water @237.15K

See ATS process model documentation's permafrost physical properties
documentation for details.

UNITS: MJ/mol
------------------------------------------------------------------------- */

#ifndef AMANZI_ENERGY_RELATIONS_IE_WATER_VAPOR_
#define AMANZI_ENERGY_RELATIONS_IE_WATER_VAPOR_

#include "Teuchos_ParameterList.hpp"

namespace Amanzi {
namespace Energy {

class IEMWaterVapor {

public:
  IEMWaterVapor(Teuchos::ParameterList& plist);

  bool IsMolarBasis() { return molar_basis_; }

  double InternalEnergy(double temp, double mol_frac_gas);
  double DInternalEnergyDT(double temp, double mol_frac_gas);
  double DInternalEnergyDomega(double temp, double mol_frac_gas);

private:
  void InitializeFromPlist_();

  Teuchos::ParameterList plist_;

  double Cv_air_; // units: MJ/(mol-K)
  double heat_vaporization_; // units: MJ/mol
  bool molar_basis_;
};

} //namespace
} //namespace

#endif
