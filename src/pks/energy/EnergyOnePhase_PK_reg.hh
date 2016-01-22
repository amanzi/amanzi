/*
  Energy

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon

  Registration of process kernel for single-phase energy equation.
*/

#include "EnergyOnePhase_PK.hh"

namespace Amanzi {
namespace Energy {

RegisteredPKFactory<EnergyOnePhase_PK> EnergyOnePhase_PK::reg_("one-phase energy");

}  // namespace Energy
}  // namespace Amanzi
