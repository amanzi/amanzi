/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov
*/

/*
  Transport PK

  Self-registering factory for MDM and multiscale porosity models.
*/

#include "MDMFactory.hh"
#include "MDM_Bear.hh"
#include "MDM_BurnettFrind.hh"
#include "MDM_Isotropic.hh"
#include "MDM_LichtnerKelkarRobinson.hh"
#include "EvaluatorMDM.hh"

#include "MultiscaleTransportPorosityFactory.hh"
#include "MultiscaleTransportPorosity_GDPM.hh"
#include "MultiscaleTransportPorosity_DPM.hh"

// explicity instantitate the static data of factory
namespace Amanzi {
namespace Utils {

template<>
Factory<Transport::MDM>::map_type* Factory<Transport::MDM>::map_;

template<>
Factory<Transport::MultiscaleTransportPorosity>::map_type*
  Factory<Transport::MultiscaleTransportPorosity>::map_;

} // namespace Utils
} // namespace Amanzi


namespace Amanzi {
namespace Transport {

Utils::RegisteredFactory<MDM, MDM_Bear> MDM_Bear::reg_("Bear");
Utils::RegisteredFactory<MDM, MDM_BurnettFrind> MDM_BurnettFrind::reg_("Burnett-Frind");
Utils::RegisteredFactory<MDM, MDM_Isotropic> MDM_Isotropic::reg_("isotropic");
Utils::RegisteredFactory<MDM, MDM_Isotropic> MDM_Isotropic::reg2_("scalar");
Utils::RegisteredFactory<MDM, MDM_LichtnerKelkarRobinson> MDM_LichtnerKelkarRobinson::reg_(
  "Lichtner-Kelkar-Robinson");

Utils::RegisteredFactory<MultiscaleTransportPorosity, MultiscaleTransportPorosity_DPM>
  MultiscaleTransportPorosity_DPM::reg_("dual porosity");

Utils::RegisteredFactory<MultiscaleTransportPorosity, MultiscaleTransportPorosity_GDPM>
  MultiscaleTransportPorosity_GDPM::reg_("generalized dual porosity");

Utils::RegisteredFactory<Evaluator, EvaluatorMDM> EvaluatorMDM::reg_("dispersion tensor");

} // namespace Transport
} // namespace Amanzi
