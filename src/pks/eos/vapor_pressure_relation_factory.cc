/*
  This is the EOS component of the ATS and Amanzi codes.
   
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon

  Self-registering factory for Vapor Pressure implementations.
*/

#include <string>
#include "vapor_pressure_relation_factory.hh"

namespace Amanzi {
namespace Relations {

// explicity instantitate the static data of Factory<VaporPressure>
// template<> 
//Amanzi::Utils::Factory<Amanzi::Relations::VaporPressureRelation>::map_type* 
//    Amanzi::Utils::Factory<Amanzi::Relations::VaporPressureRelation>::map_;

// method for instantiating VaporPressure implementations
Teuchos::RCP<VaporPressureRelation> VaporPressureRelationFactory::createVaporPressure(Teuchos::ParameterList& plist) {
  std::string eos_typename = plist.get<std::string>("vapor pressure model type");
  return Teuchos::rcp(CreateInstance(eos_typename, plist));
};

}  // namespace Relations
}  // namespace Amanzi

