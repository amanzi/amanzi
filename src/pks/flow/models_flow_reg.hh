/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

/*
  Flow PK

  Self-registering factory for water retention and multiscale porosity models.
*/

#include "SpecificStorage_Constant.hh"
#include "SpecificStorage_Standard.hh"

#include "WRMFactory.hh"
#include "WRM_BrooksCorey.hh"
#include "WRM_saturated.hh"
#include "WRM_vanGenuchten.hh"

#include "MultiscaleFlowPorosityFactory.hh"
#include "MultiscaleFlowPorosity_DPM.hh"
#include "MultiscaleFlowPorosity_GDPM.hh"

namespace Amanzi {
namespace Utils {

// explicity instantitate the static data of factory
template <>
Factory<Flow::SpecificStorage>::map_type* Factory<Flow::SpecificStorage>::map_;

template <>
Factory<Flow::WRM>::map_type* Factory<Flow::WRM>::map_;

template <>
Factory<Flow::MultiscaleFlowPorosity>::map_type* Factory<Flow::MultiscaleFlowPorosity>::map_;

} // namespace Utils
} // namespace Amanzi


namespace Amanzi {
namespace Flow {

Utils::RegisteredFactory<SpecificStorage, SpecificStorage_Constant> SpecificStorage_Constant::factory_("constant");
Utils::RegisteredFactory<SpecificStorage, SpecificStorage_Standard> SpecificStorage_Standard::factory_("standard");

Utils::RegisteredFactory<WRM, WRM_BrooksCorey> WRM_BrooksCorey::factory_("Brooks Corey");
Utils::RegisteredFactory<WRM, WRM_vanGenuchten> WRM_vanGenuchten::factory_("van Genuchten");
Utils::RegisteredFactory<WRM, WRM_saturated> WRM_saturated::factory_("saturated");

Utils::RegisteredFactory<MultiscaleFlowPorosity, MultiscaleFlowPorosity_DPM>
  MultiscaleFlowPorosity_DPM::factory_("dual porosity");
Utils::RegisteredFactory<MultiscaleFlowPorosity, MultiscaleFlowPorosity_GDPM>
  MultiscaleFlowPorosity_GDPM::factory_("generalized dual porosity");

} // namespace Flow
} // namespace Amanzi
