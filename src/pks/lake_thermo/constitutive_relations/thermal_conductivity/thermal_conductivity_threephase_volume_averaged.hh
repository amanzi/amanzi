/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Simple model of three-phase thermal conductivity, based upon volume-averaging
of four consitutive components.

Usage:

  <ParameterList name="Thermal Conductivity Model">
    <Parameter name="Thermal Conductivity Type" type="string" value="three-phase volume averaged"/>
    <Parameter name="thermal conductivity of soil" type="double" value=""/>
    <Parameter name="thermal conductivity of liquid" type="double" value=""/>
    <Parameter name="thermal conductivity of gas" type="double" value=""/>
    <Parameter name="thermal conductivity of ice" type="double" value=""/>
  </ParameterList>

Units: ????
------------------------------------------------------------------------- */

#ifndef PK_ENERGY_RELATIONS_THERMAL_CONDUCTIVITY_THREEPHASE_VOLUME_AVERAGED_HH_
#define PK_ENERGY_RELATIONS_THERMAL_CONDUCTIVITY_THREEPHASE_VOLUME_AVERAGED_HH_

#include "Teuchos_ParameterList.hpp"

#include "Factory.hh"
#include "thermal_conductivity_threephase.hh"

namespace Amanzi {
namespace Energy {

class ThermalConductivityThreePhaseVolumeAveraged : public ThermalConductivityThreePhase {

public:
  ThermalConductivityThreePhaseVolumeAveraged(Teuchos::ParameterList& plist);

  double ThermalConductivity(double porosity, double sat_liq, double sat_ice, double temp);

private:
  void InitializeFromPlist_();

  Teuchos::ParameterList plist_;

  double k_soil_;
  double k_ice_;
  double k_liquid_;
  double k_gas_;

private:
  static Utils::RegisteredFactory<ThermalConductivityThreePhase,
                                  ThermalConductivityThreePhaseVolumeAveraged> factory_;

};

}
}

#endif
