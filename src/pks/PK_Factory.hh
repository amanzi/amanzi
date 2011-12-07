/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

#ifndef _PK_FACTORY_HH_
#define _PK_FACTORY_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "State.hh"
#include "PK.hh"

namespace Amanzi
{

  class PK_Factory {

  public:
    Teuchos::RCP<PK> create_pk(Teuchos::ParameterList, Teuchos::RCP<State>);
  };
}

#endif
