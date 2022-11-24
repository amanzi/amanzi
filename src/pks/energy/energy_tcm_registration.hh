/*
  Energy

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon

  Self-registering factory for TCM implementations.
*/

#include "TCMFactory_TwoPhase.hh"
#include "TCM_PetersLidard_TwoPhase.hh"
#include "TCM_WetDry_TwoPhase.hh"

// explicity instantitate the static data of Factory<TCM_TwoPhase>
namespace Amanzi {
namespace Utils {

template <>
Factory<Energy::TCM_TwoPhase>::map_type* Factory<Energy::TCM_TwoPhase>::map_;

} // namespace Utils
} // namespace Amanzi


namespace Amanzi {
namespace Energy {

// linear interpolant of thermal conductivity.
Utils::RegisteredFactory<TCM_TwoPhase, TCM_PetersLidard_TwoPhase>
  TCM_PetersLidard_TwoPhase::factory_("two-phase Peters-Lidard");

// simple model of two-phase thermal conductivity, based upon:
// - Interpolation between saturated and dry conductivities via a Kersten number.
// - Power-law Kersten number.
Utils::RegisteredFactory<TCM_TwoPhase, TCM_WetDry_TwoPhase>
  TCM_WetDry_TwoPhase::factory_("two-phase wet/dry");

} // namespace Energy
} // namespace Amanzi
