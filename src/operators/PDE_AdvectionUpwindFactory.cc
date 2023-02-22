/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Operators

  Base factory for advection operators.
*/

#include "errors.hh"

#include "BCs.hh"
#include "OperatorDefs.hh"
#include "PDE_AdvectionUpwindFactory.hh"
#include "PDE_AdvectionUpwind.hh"

namespace Amanzi {
namespace Operators {

/* ******************************************************************
* Create operator
****************************************************************** */
Teuchos::RCP<PDE_AdvectionUpwind>
PDE_AdvectionUpwindFactory::Create(const Teuchos::RCP<Operator>& global_op)
{
  std::string name = oplist_.get<std::string>("discretization primary");
  bool fractured_matrix = oplist_.isParameter("fracture");

  bool manifolds = false;
  if (oplist_.isParameter("manifolds")) manifolds = oplist_.get<bool>("manifolds");

  Teuchos::RCP<PDE_AdvectionUpwind> op;

  if (global_op == Teuchos::null) {
    if (manifolds) {
      op = Teuchos::rcp(new PDE_AdvectionUpwindDFN(oplist_, mesh_));
    } else if (fractured_matrix) {
      oplist_.set<std::string>("name", "AdvectionFracturedMatrix: FACE_CELL");
      op = Teuchos::rcp(new PDE_AdvectionUpwindFracturedMatrix(oplist_, mesh_));
    } else {
      op = Teuchos::rcp(new PDE_AdvectionUpwind(oplist_, mesh_));
    }
    
  } else {
    if (manifolds) {
      op = Teuchos::rcp(new PDE_AdvectionUpwindDFN(oplist_, global_op));
    } else if (fractured_matrix) {
      oplist_.set<std::string>("name", "AdvectionFracturedMatrix: FACE_CELL");
      op = Teuchos::rcp(new PDE_AdvectionUpwindFracturedMatrix(oplist_, global_op));
    } else {
      op = Teuchos::rcp(new PDE_AdvectionUpwind(oplist_, global_op));
    }
  }

  return op;
}


/* ******************************************************************
* Create operator
****************************************************************** */
Teuchos::RCP<PDE_AdvectionUpwind> PDE_AdvectionUpwindFactory::Create(
    Teuchos::ParameterList& oplist, const Teuchos::RCP<const AmanziMesh::Mesh>& mesh)
{
  Teuchos::RCP<PDE_AdvectionUpwind> op;

  if (oplist.isParameter("fracture")) {
    oplist.set<std::string>("name", "AdvectionFracturedMatrix: FACE_CELL");
    op = Teuchos::rcp(new PDE_AdvectionUpwindFracturedMatrix(oplist, mesh));
  } else if (oplist.isParameter("single domain")) {
    op = Teuchos::rcp(new PDE_AdvectionUpwind(oplist, mesh));
  } else {
    op = Teuchos::rcp(new PDE_AdvectionUpwindDFN(oplist, mesh));
  }
  return op;
}


/* ******************************************************************
* Create operator
****************************************************************** */
Teuchos::RCP<PDE_AdvectionUpwind> PDE_AdvectionUpwindFactory::Create(
    Teuchos::ParameterList& oplist, const Teuchos::RCP<Operator>& global_op)
{
  Teuchos::RCP<PDE_AdvectionUpwind> op;

  if (oplist.isParameter("fracture")) {
    oplist.set<std::string>("name", "AdvectionFracturedMatrix: FACE_CELL");
    op = Teuchos::rcp(new PDE_AdvectionUpwindFracturedMatrix(oplist, global_op));
  } else if (oplist.isParameter("single domain")) {
    op = Teuchos::rcp(new PDE_AdvectionUpwind(oplist, global_op));
  } else {
    op = Teuchos::rcp(new PDE_AdvectionUpwindDFN(oplist, global_op));
  }
  return op;
}

} // namespace Operators
} // namespace Amanzi
