/*
  Energy

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon

  Process kernel for energy equation for Richard's flow.
*/

#include "EnergyTwoPhase_PK.hh"

namespace Amanzi {
namespace Energy {

RegisteredPKFactory<EnergyTwoPhase_PK> EnergyTwoPhase_PK::reg_("two-phase energy");

}  // namespace Energy
}  // namespace Amanzi
