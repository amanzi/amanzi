/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

#ifndef _EOS_FACTORY_HH_
#define _EOS_FACTORY_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "EOS.hh"

namespace Amanzi
{

  class PK_Factory {

  public:
    Teuchos::RCP<PK> create_pk(Teuchos::ParameterList, Teuchos::RCP<State>&);
  };
}

#endif
