/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov
*/

/*
  Shallow Water PK

*/

#include "NumericalFluxFactory.hh"
#include "NumericalFlux_CentralUpwind.hh"
#include "NumericalFlux_Rusanov.hh"

namespace Amanzi {
namespace Utils {

// explicity instantitate the static data of factory
template <>
Factory<ShallowWater::NumericalFlux>::map_type* Factory<ShallowWater::NumericalFlux>::map_;

} // namespace Utils
} // namespace Amanzi


namespace Amanzi {
namespace ShallowWater {

Utils::RegisteredFactory<NumericalFlux, NumericalFlux_Rusanov>
  NumericalFlux_Rusanov::factory_("Rusanov");
Utils::RegisteredFactory<NumericalFlux, NumericalFlux_CentralUpwind>
  NumericalFlux_CentralUpwind::factory_("central upwind");

} // namespace ShallowWater
} // namespace Amanzi
