/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------

   ATS
   Author: Ethan Coon

   Self-registering factory for IEM implementations.
   ------------------------------------------------------------------------- */

#include <string>
#include "iem_factory.hh"

namespace Amanzi {
namespace Energy {
namespace EnergyRelations {

// explicity instantitate the static data of Factory<IEM>
//template<> Factory<IEM>::map_type* Factory<IEM>::map_;
template<> Utils::Factory<InternalEnergyModel>::map_type* Utils::Factory<InternalEnergyModel>::map_;

// method for instantiating IEM implementations
Teuchos::RCP<InternalEnergyModel> IEMFactory::createIEM(Teuchos::ParameterList& plist) {
  std::string iem_typename = plist.get<std::string>("Internal Energy Model Type");
  return Teuchos::rcp(CreateInstance(iem_typename, plist));
};

} // namespace
} // namespace
} // namespace

