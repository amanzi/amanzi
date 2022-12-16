/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

/*
  Flow PK

  Self-registering factory for flow PKs.
*/

#include "Darcy_PK.hh"
#include "Richards_PK.hh"
#include "WaterStorage.hh"

namespace Amanzi {
namespace Flow {

RegisteredPKFactory<Darcy_PK> Darcy_PK::reg_("darcy");
RegisteredPKFactory<Richards_PK> Richards_PK::reg_("richards");
Utils::RegisteredFactory<Evaluator, WaterStorage> WaterStorage::reg_("water storage");

} // namespace Flow
} // namespace Amanzi
