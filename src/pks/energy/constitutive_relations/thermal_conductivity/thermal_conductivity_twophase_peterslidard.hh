/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Simple model of two-phase thermal conductivity, based upon:

- Interpolation between saturated and dry conductivities via a Kersten number.
- Power-law Kersten number.
- Emperical fit for dry conductivity from Peters-Lidard et al '98.

See ATS process model documentation's permafrost model for details.

Usage:

  <ParameterList name="Thermal Conductivity Model">
    <Parameter name="Thermal Conductivity Type" type="string" value="two-phase Peters-Lidard"/>
    <Parameter name="thermal conductivity of rock" type="double" value=""/>
    <Parameter name="thermal conductivity of liquid" type="double" value=""/>
    <Parameter name="thermal conductivity of gas" type="double" value=""/>

    <Parameter name="unsaturated alpha" type="double" value="1.0"/>
    <Parameter name="epsilon" type="double" value="1.e-10"/>
  </ParameterList>

Units: ????
------------------------------------------------------------------------- */

#ifndef PK_ENERGY_RELATIONS_THERMAL_CONDUCTIVITY_TWOPHASE_PETERSLIDARD_HH_
#define PK_ENERGY_RELATIONS_THERMAL_CONDUCTIVITY_TWOPHASE_PETERSLIDARD_HH_

#include "Teuchos_ParameterList.hpp"

#include "factory.hh"
#include "thermal_conductivity_twophase.hh"

namespace Amanzi {
namespace Energy {
namespace EnergyRelations {

class ThermalConductivityTwoPhasePetersLidard : public ThermalConductivityTwoPhase {

public:
  ThermalConductivityTwoPhasePetersLidard(Teuchos::ParameterList& plist);

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

private:
  static Utils::RegisteredFactory<ThermalConductivityTwoPhase,
                                  ThermalConductivityTwoPhasePetersLidard> factory_;

};

}
}
}

#endif
