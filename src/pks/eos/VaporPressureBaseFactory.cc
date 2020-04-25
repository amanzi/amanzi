/*
  EOS
   
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon

  Self-registering factory for vapor pressure implementations.
*/

#include <string>
#include "VaporPressureBaseFactory.hh"

namespace Amanzi {
namespace AmanziEOS {

// explicity instantitate the static data of Factory<VaporPressure>
// template<> 
// Amanzi::Utils::Factory<Amanzi::AmanziEOS::VaporPressureRelation>::map_type* 
//    Amanzi::Utils::Factory<Amanzi::AmanziEOS::VaporPressureRelation>::map_;

// method for instantiating VaporPressure implementations
Teuchos::RCP<VaporPressure_Base> VaporPressureBaseFactory::CreateVaporPressure(Teuchos::ParameterList& plist) {
  std::string eos_typename = plist.get<std::string>("vapor pressure model type");
  return Teuchos::rcp(CreateInstance(eos_typename, plist));
};

}  // namespace AmanziEOS
}  // namespace Amanzi

