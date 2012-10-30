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
#include "factory.hh"

namespace Amanzi {
namespace Relations {

class EOSFactory : public Utils::Factory<EOS> {

public:
  Teuchos::RCP<EOS> createEOS(Teuchos::ParameterList& plist);
};

} // namespace
} // namespace

#endif
