/*
  This is the operators component of the Amanzi code.

  License: BSD
  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#ifndef AMANZI_OPERATOR_UPWIND_FACTORY_HH_
#define AMANZI_OPERATOR_UPWIND_FACTORY_HH_

#include <string>
#include <vector>

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "Upwind.hh"
#include "UpwindStandard.hh"
#include "UpwindMFD.hh"

namespace Amanzi {
namespace Operators {

template<class Model>
class UpwindFactory {
 public:
  UpwindFactory() {};
  ~UpwindFactory() {};

  Teuchos::RCP<Upwind<Model> > Create(Teuchos::RCP<const AmanziMesh::Mesh> mesh,
                               Teuchos::RCP<const Model> model,
                               Teuchos::ParameterList& plist);
};

/* ******************************************************************
* The base class for all upwind methods. 
****************************************************************** */ 
template<class Model>
Teuchos::RCP<Upwind<Model> > UpwindFactory<Model>::Create(
    Teuchos::RCP<const AmanziMesh::Mesh> mesh,
    Teuchos::RCP<const Model> model, Teuchos::ParameterList& plist)
{
  if (!plist.isParameter("upwind method")) {
    Errors::Message msg("UpwindFactory: parameter \"upwind method\" is missing");
    Exceptions::amanzi_throw(msg);
  }

  std::string name = plist.get<std::string>("upwind method");
  if (name == "standard") {
    Teuchos::ParameterList sublist = plist.sublist("upwind standard parameters");
    Teuchos::RCP<UpwindStandard<Model> > upwind = Teuchos::rcp(new UpwindStandard<Model>(mesh, model));
    upwind->Init(sublist);
    return upwind;
  } else if (name == "mfd") {
    Teuchos::ParameterList sublist = plist.sublist("upwind mfd parameters");
    Teuchos::RCP<UpwindMFD<Model> > upwind = Teuchos::rcp(new UpwindMFD<Model>(mesh, model));
    upwind->Init(sublist);
    return upwind;
  } else {
    std::stringstream msgstream;
    msgstream << "UpwindFactory: upwind method \"" << name << "\" is not supported.";
    Errors::Message msg(msgstream.str());
    Exceptions::amanzi_throw(msg);
  }
}

}  // namespace Operators
}  // namespace Amanzi

#endif
