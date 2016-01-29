/*
  Energy

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

#include "TCM_WetDry_TwoPhase.hh"

namespace Amanzi {
namespace Energy {

Utils::RegisteredFactory<TCM_TwoPhase,TCM_WetDry_TwoPhase> TCM_WetDry_TwoPhase::factory_("two-phase wet/dry");

}  // namespace Energy
}  // namespace Amanzi
