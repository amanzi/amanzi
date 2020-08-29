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

#include "SaturatedVaporPressureFactory.hh"

namespace Amanzi {
namespace AmanziEOS {

// explicity instantitate the static data of Factory<VaporPressure>
// template<> 
// Amanzi::Utils::Factory<Amanzi::AmanziEOS::VaporPressureRelation>::map_type* 
//    Amanzi::Utils::Factory<Amanzi::AmanziEOS::VaporPressureRelation>::map_;

// method for instantiating vapor pressure implementations
Teuchos::RCP<SaturatedVaporPressure> SaturatedVaporPressureFactory::CreateVaporPressure(Teuchos::ParameterList& plist) {
  std::string model = plist.get<std::string>("vapor pressure model type");
  return Teuchos::rcp(CreateInstance(model, plist));
};

}  // namespace AmanziEOS
}  // namespace Amanzi

