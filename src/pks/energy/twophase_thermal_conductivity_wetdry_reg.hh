/*
  This is the energy component of the ATS and Amanzi codes. 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon

  Simple model of two-phase thermal conductivity, based upon:
    - Interpolation between saturated and dry conductivities via a Kersten number.
    - Power-law Kersten number.
  See ATS process model documentation's permafrost model for details.
*/

#include "twophase_thermal_conductivity_wetdry.hh"

namespace Amanzi {
namespace Energy {

// registry of method
Utils::RegisteredFactory<ThermalConductivityTwoPhase,ThermalConductivityTwoPhaseWetDry>
    ThermalConductivityTwoPhaseWetDry::factory_("two-phase wet/dry");

}  // namespace Energy
}  // namespace Amanzi
