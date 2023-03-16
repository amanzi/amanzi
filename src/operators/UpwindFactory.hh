/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Operators

*/

#ifndef AMANZI_OPERATOR_UPWIND_FACTORY_HH_
#define AMANZI_OPERATOR_UPWIND_FACTORY_HH_

#include <string>
#include <vector>

// TPLs
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

// Operators
#include "Upwind.hh"
#include "UpwindArithmeticAverage.hh"
#include "UpwindDivK.hh"
#include "UpwindFlux.hh"
#include "UpwindFluxAndGravity.hh"
#include "UpwindFluxManifolds.hh"
#include "UpwindGravity.hh"
#include "UpwindSecondOrder.hh"

namespace Amanzi {
namespace Operators {

class UpwindFactory {
 public:
  Teuchos::RCP<Upwind>
  Create(Teuchos::RCP<const AmanziMesh::Mesh> mesh, Teuchos::ParameterList& plist);
};


/* ******************************************************************
* The base class for all upwind methods.
****************************************************************** */
inline Teuchos::RCP<Upwind>
UpwindFactory::Create(Teuchos::RCP<const AmanziMesh::Mesh> mesh, Teuchos::ParameterList& plist)
{
  if (!plist.isParameter("upwind method")) {
    Errors::Message msg("UpwindFactory: parameter \"upwind method\" is missing");
    Exceptions::amanzi_throw(msg);
  }
  std::string name = plist.get<std::string>("upwind method");
  Teuchos::ParameterList sublist = plist.sublist("upwind parameters");

  bool manifolds =
    (mesh->space_dimension() != mesh->manifold_dimension()) && sublist.get<bool>("manifolds", true);

  if (name == "upwind: darcy velocity" && !manifolds) {
    Teuchos::RCP<UpwindFlux> upwind = Teuchos::rcp(new UpwindFlux(mesh));
    upwind->Init(sublist);
    return upwind;
  } else if (name == "upwind: darcy velocity" && manifolds) {
    Teuchos::RCP<UpwindFluxManifolds> upwind = Teuchos::rcp(new UpwindFluxManifolds(mesh));
    upwind->Init(sublist);
    return upwind;
  } else if (name == "upwind: gravity") {
    Teuchos::RCP<UpwindGravity> upwind = Teuchos::rcp(new UpwindGravity(mesh));
    upwind->Init(sublist);
    return upwind;
  } else if (name == "upwind: amanzi") {
    Teuchos::RCP<UpwindDivK> upwind = Teuchos::rcp(new UpwindDivK(mesh));
    upwind->Init(sublist);
    return upwind;
  } else if (name == "upwind: amanzi new") {
    Teuchos::RCP<UpwindFluxAndGravity> upwind = Teuchos::rcp(new UpwindFluxAndGravity(mesh));
    upwind->Init(sublist);
    return upwind;
  } else if (name == "upwind: second-order") {
    Teuchos::RCP<UpwindSecondOrder> upwind = Teuchos::rcp(new UpwindSecondOrder(mesh));
    upwind->Init(sublist);
    return upwind;
  } else if (name == "other: arithmetic average") {
    Teuchos::RCP<UpwindArithmeticAverage> upwind = Teuchos::rcp(new UpwindArithmeticAverage(mesh));
    upwind->Init(sublist);
    return upwind;
  } else {
    std::stringstream msgstream;
    msgstream << "UpwindFactory: method \"" << name << "\" is not supported.";
    Errors::Message msg(msgstream.str());
    Exceptions::amanzi_throw(msg);
  }

  return Teuchos::null;
}

} // namespace Operators
} // namespace Amanzi

#endif
