/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

/*
  Mechanics PKs

  Self-registering factory for mechanics PKs.
*/

#include "MechanicsElasticity_PK.hh"
#include "MechanicsSmallStrain_PK.hh"
#include "MechanicsFracturedMatrix_PK.hh"

namespace Amanzi {
namespace Mechanics {

RegisteredPKFactory<MechanicsElasticity_PK> MechanicsElasticity_PK::reg_("elastic");
RegisteredPKFactory<MechanicsSmallStrain_PK> MechanicsSmallStrain_PK::reg_("small strain");
RegisteredPKFactory<MechanicsFracturedMatrix_PK>
  MechanicsFracturedMatrix_PK::reg_("mechanics matrix fracture");

} // namespace Mechanics
} // namespace Amanzi
