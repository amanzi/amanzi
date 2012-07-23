/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------

   ATS
   Author: Ethan Coon

   Self-registering factory for EOS implementations.
   ------------------------------------------------------------------------- */

#ifndef _PK_FLOW_EOS_FACTORY_HH_
#define _PK_FLOW_EOS_FACTORY_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "eos.hh"
#include "factory_with_state.hh"

namespace Amanzi {
namespace Flow {
namespace FlowRelations {

class EOSFactory : public Utils::FactoryWithState<EOS> {

public:
  Teuchos::RCP<EOS> createEOS(Teuchos::ParameterList& plist, const Teuchos::Ptr<State>& S);
};

} // namespace
} // namespace
} // namespace

#endif
