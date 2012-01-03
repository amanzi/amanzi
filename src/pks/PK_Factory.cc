/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Implementation of the process kernel factory.  Currently this must be manually
edited each time a new PK needs to be registered with the factory... this
should eventually get smarter.
------------------------------------------------------------------------- */

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "errors.hh"
#include "State.hh"
#include "PK.hh"
#include "PK_Factory.hh"
#include "WeakMPC.hh"
// #include "DarcyPK.hh"
// #include "RichardsPK.hh"
// #include "PermafrostPK.hh"
#include "NullEnergyPK.hh"

namespace Amanzi {

Teuchos::RCP<PK> PK_Factory::create_pk(Teuchos::ParameterList plist,
        Teuchos::RCP<State>& S, Teuchos::RCP<TreeVector>& soln) {
  std::string pk_type = plist.get<string>("PK model");

  // this could be much fancier and smarter, but for now...
  if (pk_type == "Weak MPC") {
    return Teuchos::rcp(new WeakMPC(plist, S, soln));
  // } else if (pk_type == "Strong MPC") {
  //   return Teuchos::rcp(new StrongMPC(plist, S, soln));
  // } else if (pk_type == "Darcy") {
  //   return Teuchos::rcp(new DarcyPK(plist, S, soln));
  // } else if (pk_type == "Richards") {
  //   return Teuchos::rcp(new RichardsPK(plist, S, soln));
  // } else if (pk_type == "Permafrost") {
  //   return Teuchos::rcp(new PermafrostPK(plist, S, soln));
  } else if (pk_type == "Constant Temperature") {
    return Teuchos::rcp(new NullEnergyPK(plist, S, soln));
  } else {
    Errors::Message message("PK_Factory: unknown PK model: "+pk_type);
    Exceptions::amanzi_throw(message);
  }
};
} // namespace
