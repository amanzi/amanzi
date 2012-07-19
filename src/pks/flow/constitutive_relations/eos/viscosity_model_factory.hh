/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------

   ATS
   Author: Ethan Coon

   Self-registering factory for ViscosityModel implementations.
   ------------------------------------------------------------------------- */

#ifndef PK_FLOW_VISCOSITY_MODEL_FACTORY_HH_
#define PK_FLOW_VISCOSITY_MODEL_FACTORY_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "viscosity_model.hh"
#include "factory_with_state.hh"

namespace Amanzi {
namespace Flow {
namespace FlowRelations {

class ViscosityModelFactory : public Utils::FactoryWithState<ViscosityModel> {

public:
  Teuchos::RCP<ViscosityModel> createViscosityModel(Teuchos::ParameterList& plist, const Teuchos::Ptr<State>& S);
};

} // namespace
} // namespace
} // namespace

#endif
