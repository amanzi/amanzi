/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Base of a Thermal Conductivity relation.

UNITS: ????
------------------------------------------------------------------------- */

#ifndef PK_ENERGY_RELATIONS_THERMAL_CONDUCTIVITY_HH_
#define PK_ENERGY_RELATIONS_THERMAL_CONDUCTIVITY_HH_

#include "Teuchos_ParameterList.hpp"

namespace Amanzi {
namespace Energy {
namespace EnergyRelations {

class TwophaseThermalConductivity {

public:
  TwophaseThermalConductivity(Teuchos::ParameterList& plist);

  double CalculateConductivity(double porosity, double sat_liq);

private:
  void InitializeFromPlist_();

  Teuchos::ParameterList plist_;

  double eps_;
  double alpha_;
  double k_rock_;
  double k_liquid_;
  double k_gas_;
  double d_;
  double rho_rock_;
};

}
}
}

#endif
