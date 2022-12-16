/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

/*
  Flow PK

  Self-registering factory for WRMmp implementations.
*/

#include "WRMmp_BrooksCorey.hh"
#include "WRMmp_Corey.hh"
#include "WRMmp_Simple.hh"
#include "WRMmp_vanGenuchten.hh"

namespace Amanzi {
namespace Utils {

// explicity instantitate the static data of factory
template <>
Factory<Multiphase::WRMmp>::map_type* Factory<Multiphase::WRMmp>::map_;

} // namespace Utils
} // namespace Amanzi


namespace Amanzi {
namespace Multiphase {

Utils::RegisteredFactory<WRMmp, WRMmp_BrooksCorey> WRMmp_BrooksCorey::factory_("Brooks Corey");
Utils::RegisteredFactory<WRMmp, WRMmp_Corey> WRMmp_Corey::factory_("Corey");
Utils::RegisteredFactory<WRMmp, WRMmp_vanGenuchten> WRMmp_vanGenuchten::factory_("van Genuchten");
Utils::RegisteredFactory<WRMmp, WRMmp_Simple> WRMmp_Simple::factory_("Simple");

} // namespace Multiphase
} // namespace Amanzi
