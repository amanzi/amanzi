/*
  This is the energy component of the ATS / Amanzi code.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon

  Simple model of two-phase thermal conductivity, based upon:
   - Interpolation between saturated and dry conductivities via a Kersten number.
   - Power-law Kersten number.
   - Emperical fit for dry conductivity from Peters-Lidard et al '98.

  See ATS process model documentation's permafrost model for details.
  Units: ????

  Usage:

  <ParameterList name="thermal conductivity parameters">
    <Parameter name="thermal conductivity type" type="string" value="two-phase Peters-Lidard"/>
    <Parameter name="thermal conductivity of rock" type="double" value=""/>
    <Parameter name="thermal conductivity of liquid" type="double" value=""/>
    <Parameter name="thermal conductivity of gas" type="double" value=""/>

    <Parameter name="unsaturated alpha" type="double" value="1.0"/>
    <Parameter name="epsilon" type="double" value="1.e-10"/>
  </ParameterList>
*/

#ifndef AMANZI_ENERGY_THERMAL_CONDUCTIVITY_TWOPHASE_PETERSLIDARD_HH_
#define AMANZI_ENERGY_THERMAL_CONDUCTIVITY_TWOPHASE_PETERSLIDARD_HH_

#include "Teuchos_ParameterList.hpp"

#include "factory.hh"
#include "twophase_thermal_conductivity.hh"

namespace Amanzi {
namespace Energy {

class ThermalConductivityTwoPhasePetersLidard : public ThermalConductivityTwoPhase {
 public:
  ThermalConductivityTwoPhasePetersLidard(Teuchos::ParameterList& plist);

  double ThermalConductivity(double porosity, double sat_liq);

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

}  // namespace Energy
}  // namespace Amanzi

#endif
