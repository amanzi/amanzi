/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "errors.hh"
#include "State.hh"
#include "PK.hh"
#include "PK_Factory.hh"
#include "WeakMPC.hh"
// #include "Darcy_PK.hh"
// #include "Richards_PK.hh"
// #include "Permafrost_PK.hh"
#include "NullEnergyPK.hh"

namespace Amanzi {

Teuchos::RCP<PK> PK_Factory::create_pk(Teuchos::ParameterList parameter_list,
					 Teuchos::RCP<State> &S) {
  std::string pk_type = parameter_list.get<string>("PK model");

  // this could be much fancier and smarter, but for now...
  if (pk_type == "Weak MPC") {
    return Teuchos::rcp(new WeakMPC(parameter_list, S));
  // } else if (pk_type == "Strong MPC") {
  //   return Teuchos::rcp(new StrongMPC(parameter_list, S));
  // } else if (pk_type == "Darcy") {
  //   return Teuchos::rcp(new Darcy_PK(parameter_list, S));
  // } else if (pk_type == "Richards") {
  //   return Teuchos::rcp(new Richards_PK(parameter_list, S));
  // } else if (pk_type == "Permafrost") {
  //   return Teuchos::rcp(new Permafrost_PK(parameter_list, S));
  } else if (pk_type == "Constant Temperature") {
    return Teuchos::rcp(new NullEnergyPK(parameter_list, S));
  } else {
    Errors::Message message("PK_Factory: unknown PK model: "+pk_type);
    Exceptions::amanzi_throw(message);
  }
};

} // namespace
