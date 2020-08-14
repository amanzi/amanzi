/*
  Flow PK 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon

  Self-registering factory for WRM implementations.
*/

#include "WRMFactory.hh"
#include "WRM_BrooksCorey.hh"
#include "WRM_fake.hh"
#include "WRM_saturated.hh"
#include "WRM_vanGenuchten.hh"

namespace Amanzi {
namespace Utils {

// explicity instantitate the static data of factory
template<> Factory<Flow::WRM>::map_type* Factory<Flow::WRM>::map_;

}  // namespace Utils
}  // namespace Amanzi


namespace Amanzi {
namespace Flow {

Utils::RegisteredFactory<WRM, WRM_BrooksCorey> WRM_BrooksCorey::factory_("Brooks Corey");
Utils::RegisteredFactory<WRM, WRM_vanGenuchten> WRM_vanGenuchten::factory_("van Genuchten");
Utils::RegisteredFactory<WRM, WRM_fake> WRM_fake::factory_("fake");
Utils::RegisteredFactory<WRM, WRM_saturated> WRM_saturated::factory_("saturated");

}  // namespace Flow
}  // namespace Amanzi


