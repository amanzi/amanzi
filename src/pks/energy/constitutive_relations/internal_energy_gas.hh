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

#ifndef INTERNAL_ENERGY_GAS_
#define INTERNAL_ENERGY_GAS_

#include "Teuchos_ParameterList.hpp"

namespace Amanzi {
namespace Energy {
namespace EnergyRelations {

class InternalEnergyGas {

public:
  InternalEnergyGas(Teuchos::ParameterList& plist);

  double CalculateInternalEnergy(double temp, double mol_frac_gas);

private:
  void InitializeFromPlist_();

  Teuchos::ParameterList plist_;

  double Cv_air_; // units: J/(mol-K)
  double heat_vaporization_; // units: J/mol
};

}
}
}

#endif
