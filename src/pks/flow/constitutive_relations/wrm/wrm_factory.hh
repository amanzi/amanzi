/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------

   ATS
   Author: Ethan Coon

   Self-registering factory for EOS implementations.
   ------------------------------------------------------------------------- */

#ifndef _PK_FLOW_WRM_FACTORY_HH_
#define _PK_FLOW_WRM_FACTORY_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "wrm.hh"
#include "factory.hh"

namespace Amanzi {
namespace Flow {
namespace FlowRelations {

class WRMFactory : public Utils::Factory<WRM> {

public:
  Teuchos::RCP<WRM> createWRM(Teuchos::ParameterList& plist);
};

} // namespace
} // namespace
} // namespace

#endif
