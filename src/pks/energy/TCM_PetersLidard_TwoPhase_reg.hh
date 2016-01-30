/*
  Energy

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon

  Linear interpolant of thermal conductivity.
*/

#include "TCM_PetersLidard_TwoPhase.hh"

namespace Amanzi {
namespace Energy {

// registry of method
Utils::RegisteredFactory<TCM_TwoPhase,TCM_PetersLidard_TwoPhase>
    TCM_PetersLidard_TwoPhase::factory_("two-phase Peters-Lidard");

}  // namespace Energy
}  // namespace Amanzi
