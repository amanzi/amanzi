/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

/*
  Energy

  Self-registering factory of models for energy PK.
*/

#include "EnergyPressureEnthalpy_PK.hh"
#include "EnergyPressureTemperature_PK.hh"
#include "EnergyTwoPhase_PK.hh"

namespace Amanzi {
namespace Energy {

RegisteredPKFactory<EnergyPressureTemperature_PK> EnergyPressureTemperature_PK::reg_("pt energy");
RegisteredPKFactory<EnergyPressureEnthalpy_PK> EnergyPressureEnthalpy_PK::reg_("ph energy");
RegisteredPKFactory<EnergyTwoPhase_PK> EnergyTwoPhase_PK::reg_("two-phase energy");

} // namespace Energy
} // namespace Amanzi
