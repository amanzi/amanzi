/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------

   ATS
   Author: Ethan Coon

   Self-registering factory for WRM implementations.
   ------------------------------------------------------------------------- */

#ifndef AMANZI_ENERGYRELATIONS_IEM_FACTORY_
#define AMANZI_ENERGYRELATIONS_IEM_FACTORY_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "internal_energy_model.hh"
#include "factory.hh"

namespace Amanzi {
namespace Energy {
namespace EnergyRelations {

class IEMFactory : public Utils::Factory<InternalEnergyModel> {

public:
  Teuchos::RCP<InternalEnergyModel> createIEM(Teuchos::ParameterList& plist);
};

} // namespace
} // namespace
} // namespace

#endif
