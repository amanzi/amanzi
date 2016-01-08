/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#ifndef AMANZI_OPERATOR_UPWIND_FACTORY_HH_
#define AMANZI_OPERATOR_UPWIND_FACTORY_HH_

#include <string>
#include <vector>

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "Upwind.hh"
#include "UpwindDivK.hh"
#include "UpwindFlux.hh"
#include "UpwindGravity.hh"
#include "UpwindSecondOrder.hh"
#include "UpwindStandard.hh"

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
  Teuchos::ParameterList sublist = plist.sublist("upwind " + name + " standard parameters");
  if (name == "standard") {
    Teuchos::RCP<UpwindFlux<Model> > upwind = Teuchos::rcp(new UpwindFlux<Model>(mesh, model));
    upwind->Init(sublist);
    return upwind;
  } else if (name == "flux") {
    Teuchos::RCP<UpwindFlux<Model> > upwind = Teuchos::rcp(new UpwindFlux<Model>(mesh, model));
    upwind->Init(sublist);
    return upwind;
  } else if (name == "gravity") {
    Teuchos::RCP<UpwindGravity<Model> > upwind = Teuchos::rcp(new UpwindGravity<Model>(mesh, model));
    upwind->Init(sublist);
    return upwind;
  } else if (name == "divk") {
    Teuchos::RCP<UpwindDivK<Model> > upwind = Teuchos::rcp(new UpwindDivK<Model>(mesh, model));
    upwind->Init(sublist);
    return upwind;
  } else if (name == "second-order") {
    Teuchos::RCP<UpwindSecondOrder<Model> > upwind = Teuchos::rcp(new UpwindSecondOrder<Model>(mesh, model));
    upwind->Init(sublist);
    return upwind;
  } else {
    std::stringstream msgstream;
    msgstream << "UpwindFactory: upwind method \"" << name << "\" is not supported.";
    Errors::Message msg(msgstream.str());
    Exceptions::amanzi_throw(msg);
  }

  return Teuchos::null;
}

}  // namespace Operators
}  // namespace Amanzi

#endif
