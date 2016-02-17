/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------

   ATS
   Author: Ethan Coon

   Self-registering factory for WRM_PERMAFROST implementations.
   ------------------------------------------------------------------------- */

#ifndef _PK_FLOW_WRM_PERMAFROST_FACTORY_HH_
#define _PK_FLOW_WRM_PERMAFROST_FACTORY_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "wrm_permafrost_model.hh"
#include "factory.hh"

namespace Amanzi {
namespace Flow {
namespace FlowRelations {

class WRMPermafrostFactory : public Utils::Factory<WRMPermafrostModel> {

public:
  Teuchos::RCP<WRMPermafrostModel> createWRMPermafrostModel(Teuchos::ParameterList& plist,
          const Teuchos::RCP<WRM>& wrm);
};

} // namespace
} // namespace
} // namespace

#endif
