/* 
  This is the EOS component of the ATS and Amanzi codes.
   
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon

  Self-registering factory for Viscosity implementations.
*/

#include <string>
#include "viscosity_relation_factory.hh"

// explicity instantitate the static data of Factory<Viscosity>
template<> 
Amanzi::Utils::Factory<Amanzi::Relations::ViscosityRelation>::map_type* 
Amanzi::Utils::Factory<Amanzi::Relations::ViscosityRelation>::map_;

namespace Amanzi {
namespace Relations {

// method for instantiating Viscosity implementations
Teuchos::RCP<ViscosityRelation> ViscosityRelationFactory::createViscosity(Teuchos::ParameterList& plist) {
  std::string visc_typename = plist.get<std::string>("viscosity relation type");
  return Teuchos::rcp(CreateInstance(visc_typename, plist));
};

}  // namespace Relations
}  // namespace Amanzi

