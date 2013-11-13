/*
  Factory for timestep control.

  License: BSD
  Authors: Ethan Coon (ecoon@lanl.gov)

*/

#include "TimestepControllerFactory.hh"

#include "TimestepControllerFixed.hh"
#include "TimestepControllerStandard.hh"
#include "TimestepControllerSmarter.hh"

namespace Amanzi {

Teuchos::RCP<TimestepController>
TimestepControllerFactory::Create(
    const Teuchos::ParameterList& slist)
{
  if (slist.isParameter("timestep controller type")) {
    std::string type = slist.get<std::string>("timestep controller type");

    if (type == "fixed") {
      if (!slist.isSublist("timestep controller fixed parameters")) {
        Errors::Message msg("TimestepControllerFactory: missing sublist \"timestep controller fixed parameters\"");
        Exceptions::amanzi_throw(msg);
      }
      Teuchos::ParameterList sublist = slist.sublist("timestep controller fixed parameters");
      return Teuchos::rcp(new TimestepControllerFixed(sublist));

    } else if (type == "standard") {
      if (!slist.isSublist("timestep controller standard parameters")) {
        Errors::Message msg("TimestepControllerFactory: missing sublist \"timestep controller standard parameters\"");
        Exceptions::amanzi_throw(msg);
      }
      Teuchos::ParameterList sublist = slist.sublist("timestep controller standard parameters");
      return Teuchos::rcp(new TimestepControllerStandard(sublist));

    } else if (type == "smarter") {
      if (!slist.isSublist("timestep controller smarter parameters")) {
        Errors::Message msg("TimestepControllerFactory: missing sublist \"timestep controller smarter parameters\"");
        Exceptions::amanzi_throw(msg);
      }
      Teuchos::ParameterList sublist = slist.sublist("timestep controller smarter parameters");
      return Teuchos::rcp(new TimestepControllerSmarter(sublist));

    } else {
      Errors::Message msg("TimestepControllerFactory: invalid value of parameter `\"timestep controller type`\"");
      Exceptions::amanzi_throw(msg);
    }
  } else {
    Errors::Message msg("TimestepControllerFactory: parameter `\"timestep controller type`\" is missing");
    Exceptions::amanzi_throw(msg);
  }
  return Teuchos::null;
}


}  // namespace
