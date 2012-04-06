/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------

   ATS
   Author: Ethan Coon

   Self-registering factory for EOS implementations.
   ------------------------------------------------------------------------- */

#ifndef _PK_FLOW_EOS_VAPOR_PRESSURE_MODEL_FACTORY_HH_
#define _PK_FLOW_EOS_VAPOR_PRESSURE_MODEL_FACTORY_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "vapor_pressure_model.hh"
#include "factory.hh"

namespace Amanzi {
namespace Flow {
namespace FlowRelations {

class VaporPressureModelFactory : public Utils::Factory<VaporPressureModel> {

public:
  Teuchos::RCP<VaporPressureModel> createVaporPressureModel(Teuchos::ParameterList& plist);
};

} // namespace
} // namespace
} // namespace

#endif
