/*
  Transport PK 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov

  Self-registering factory for MDM implementations.
*/

#include "MDMFactory.hh"
#include "MDM_Bear.hh"
#include "MDM_BurnettFrind.hh"
#include "MDM_Isotropic.hh"
#include "MDM_LichtnerKelkarRobinson.hh"

// explicity instantitate the static data of factory
namespace Amanzi {
namespace Utils {

template<> Factory<Transport::MDM>::map_type* Factory<Transport::MDM>::map_;

}  // namespace Utils
}  // namespace Amanzi


namespace Amanzi {
namespace Transport {

Utils::RegisteredFactory<MDM, MDM_Bear> MDM_Bear::factory_("Bear");
Utils::RegisteredFactory<MDM, MDM_BurnettFrind> MDM_BurnettFrind::factory_("Burnett-Frind");
Utils::RegisteredFactory<MDM, MDM_Isotropic> MDM_Isotropic::factory_("scalar");
Utils::RegisteredFactory<MDM, MDM_LichtnerKelkarRobinson> 
    MDM_LichtnerKelkarRobinson::factory_("Lichtner-Kelkar-Robinson");

}  // namespace Transport
}  // namespace Amanzi

