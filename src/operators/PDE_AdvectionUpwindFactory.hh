/*
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Simple factory of advection schemes.
*/

#ifndef AMANZI_OPERATOR_PDE_ADVECTION_UPWIND_FACTORY_HH_
#define AMANZI_OPERATOR_PDE_ADVECTION_UPWIND_FACTORY_HH_

#include <vector>

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "PDE_AdvectionUpwind.hh"
#include "PDE_AdvectionUpwindDFN.hh"
#include "PDE_AdvectionUpwindFracturedMatrix.hh"

namespace Amanzi {
namespace Operators {

struct PDE_AdvectionUpwindFactory {
  Teuchos::RCP<PDE_AdvectionUpwind>
  Create(Teuchos::ParameterList& oplist,
         const Teuchos::RCP<const AmanziMesh::Mesh>& mesh)
  {
    if (oplist.isParameter("fracture")) {
      auto op = Teuchos::rcp(new PDE_AdvectionUpwindFracturedMatrix(oplist, mesh));
      return op;
    } else if (oplist.isParameter("single domain")) {
      auto op = Teuchos::rcp(new PDE_AdvectionUpwind(oplist, mesh));
      return op;
    } else {
      auto op = Teuchos::rcp(new PDE_AdvectionUpwindDFN(oplist, mesh));
      return op;
    }
  }

  Teuchos::RCP<PDE_AdvectionUpwind>
  Create(Teuchos::ParameterList& oplist,
         const Teuchos::RCP<Operator>& global_op)
  {
    if (oplist.isParameter("fracture")) {
      auto op = Teuchos::rcp(new PDE_AdvectionUpwindFracturedMatrix(oplist, global_op));
      return op;
    } else if (oplist.isParameter("single domain")) {
      auto op = Teuchos::rcp(new PDE_AdvectionUpwind(oplist, global_op));
      return op;
    } else {
      auto op = Teuchos::rcp(new PDE_AdvectionUpwindDFN(oplist, global_op));
      return op;
    }
  }
};

}  // namespace Operators
}  // namespace Amanzi

#endif
