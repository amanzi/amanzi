/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------

   ATS
   Author: Ethan Coon

   Self-registering factory for Viscosity implementations.
   ------------------------------------------------------------------------- */

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

} // namespace
} // namespace

