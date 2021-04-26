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
Amanzi::Utils::Factory<Amanzi::AmanziEOS::ViscosityBase>::map_type* 
    Amanzi::Utils::Factory<Amanzi::AmanziEOS::ViscosityBase>::map_;

namespace Amanzi {
namespace AmanziEOS {

// method for instantiating Viscosity implementations
Teuchos::RCP<ViscosityBase> ViscosityBaseFactory::CreateViscosity(Teuchos::ParameterList& plist) {
  std::string type = plist.get<std::string>("viscosity relation type");
  return Teuchos::rcp(CreateInstance(type, plist));
};

}  // namespace AmanziEOS
}  // namespace Amanzi

