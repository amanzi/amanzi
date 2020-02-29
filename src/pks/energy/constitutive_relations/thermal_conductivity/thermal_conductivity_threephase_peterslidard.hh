/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Simple model of three-phase thermal conductivity, based upon:

- Interpolation between saturated and dry conductivities via a Kersten number.
- Power-law Kersten number.
- Emperical fit for dry conductivity from Peters-Lidard et al '98.

See ATS process model documentation's permafrost model for details.

Usage:

  <ParameterList name="Thermal Conductivity Model">
    <Parameter name="Thermal Conductivity Type" type="string" value="three-phase Peters-Lidard"/>
    <Parameter name="thermal conductivity of soil" type="double" value=""/>
    <Parameter name="thermal conductivity of liquid" type="double" value=""/>
    <Parameter name="thermal conductivity of gas" type="double" value=""/>

    <Parameter name="unsaturated alpha" type="double" value="1.0"/>
    <Parameter name="epsilon" type="double" value="1.e-10"/>
  </ParameterList>

Units: ????
------------------------------------------------------------------------- */

#ifndef PK_ENERGY_RELATIONS_THERMAL_CONDUCTIVITY_THREEPHASE_PETERSLIDARD_HH_
#define PK_ENERGY_RELATIONS_THERMAL_CONDUCTIVITY_THREEPHASE_PETERSLIDARD_HH_

#include "Teuchos_ParameterList.hpp"

#include "Factory.hh"
#include "thermal_conductivity_threephase.hh"

namespace Amanzi {
namespace Energy {

class ThermalConductivityThreePhasePetersLidard : public ThermalConductivityThreePhase {

public:
  ThermalConductivityThreePhasePetersLidard(Teuchos::ParameterList& plist);

  double ThermalConductivity(double porosity, double sat_liq, double sat_ice, double temp);

private:
  void InitializeFromPlist_();

  Teuchos::ParameterList plist_;

  double eps_;
  double alpha_f_;
  double alpha_u_;
  double k_soil_;
  double k_ice_;
  double k_liquid_;
  double k_gas_;
  double d_;

private:
  static Utils::RegisteredFactory<ThermalConductivityThreePhase,
                                  ThermalConductivityThreePhasePetersLidard> factory_;

};

}
}

#endif
