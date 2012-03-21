/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Simple model of two-phase thermal conductivity, based upon:

- Interpolation between saturated and dry conductivities via a Kersten number.
- Power-law Kersten number.
- Emperical fit for dry conductivity.

See ATS process model documentation's permafrost model for details.

UNITS: ????
------------------------------------------------------------------------- */

#ifndef TWOPHASE_THERMAL_CONDUCTIVITY_HH_
#define TWOPHASE_THERMAL_CONDUCTIVITY_HH_

#include "Teuchos_ParameterList.hpp"

namespace Amanzi {
namespace Energy {
namespace EnergyRelations {

class TwophaseThermalConductivity {

public:
  TwophaseThermalConductivity(Teuchos::ParameterList& plist);

  double CalculateConductivity(double porosity, double sat_liq, double dens_rock);

private:
  void InitializeFromPlist_();

  Teuchos::ParameterList plist_;

  double eps_;
  double alpha_;
  double k_rock_;
  double k_liquid_;
  double k_gas_;
  double d_;
};

}
}
}

#endif
