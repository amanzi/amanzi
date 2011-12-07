/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "errors.hh"
#include "State.hh"
#include "PK.hh"
#include "PK_Factory.hh"
#include "Darcy_PK.hh"

namespace Amanzi {

Teuchos::RCP<PK> PK_Factory::create_pk(Teuchos::ParameterList parameter_list,
					 Teuchos::RCP<State> S) {
  std::string pk_type = parameter_list.get<string>("PK model");

  // this could be much fancier and smarter, but for now...
  if (pk_type == "Darcy") {
    return Teuchos::rcp(new Darcy_PK(parameter_list, S));
  } else {
    Errors::Message message("PK_Factory: unknown PK model: "+pk_type);
    Exceptions::amanzi_throw(message);
  }
};

} // namespace
