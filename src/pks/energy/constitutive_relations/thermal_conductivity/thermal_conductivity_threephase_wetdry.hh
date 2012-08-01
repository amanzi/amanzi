/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Simple model of three-phase thermal conductivity, based upon:

- Interpolation between saturated and dry conductivities via a Kersten number.
- Power-law Kersten number.

See ATS process model documentation's permafrost model for details.

Usage:

  <ParameterList name="Thermal Conductivity Model">
    <Parameter name="Thermal Conductivity Type" type="string" value="three-phase wet/dry"/>

    <Parameter name="thermal conductivity, wet" type="double" value=""/>
    <Parameter name="thermal conductivity, dry" type="double" value=""/>

    <Parameter name="epsilon" type="double" value="1.e-10"/>
    <Parameter name="unsaturated alpha" type="double" value="1.0"/>
  </ParameterList>

Units: ????
------------------------------------------------------------------------- */

#ifndef PK_ENERGY_RELATIONS_THERMAL_CONDUCTIVITY_THREEPHASE_WETDRY_HH_
#define PK_ENERGY_RELATIONS_THERMAL_CONDUCTIVITY_THREEPHASE_WETDRY_HH_

#include "Teuchos_ParameterList.hpp"

#include "factory.hh"
#include "thermal_conductivity_threephase.hh"

namespace Amanzi {
namespace Energy {
namespace EnergyRelations {

class ThermalConductivityThreePhaseWetDry : public ThermalConductivityThreePhase {

public:
  ThermalConductivityThreePhaseWetDry(Teuchos::ParameterList& plist);

  double CalculateConductivity(double porosity, double sat_liq, double sat_ice);

private:
  void InitializeFromPlist_();

  Teuchos::ParameterList plist_;

  double eps_;
  double alpha_f_;
  double alpha_u_;
  double k_sat_f_;
  double k_sat_u_;
  double k_dry_;

private:
  static Utils::RegisteredFactory<ThermalConductivityThreePhase,
                                  ThermalConductivityThreePhaseWetDry> factory_;

};

}
}
}

#endif
