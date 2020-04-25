/* 
  EOS
   
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon

  Self-registering factory for viscosity models.
*/

#include <string>

#include "ViscosityBaseFactory.hh"

// explicity instantitate the static data of Factory<Viscosity>
template<> 
Amanzi::Utils::Factory<Amanzi::AmanziEOS::Viscosity_Base>::map_type* 
    Amanzi::Utils::Factory<Amanzi::AmanziEOS::Viscosity_Base>::map_;

namespace Amanzi {
namespace AmanziEOS {

// method for instantiating Viscosity implementations
Teuchos::RCP<Viscosity_Base> ViscosityBaseFactory::CreateViscosity(Teuchos::ParameterList& plist) {
  std::string visc_typename = plist.get<std::string>("viscosity relation type");
  return Teuchos::rcp(CreateInstance(visc_typename, plist));
};

}  // namespace AmanziEOS
}  // namespace Amanzi

