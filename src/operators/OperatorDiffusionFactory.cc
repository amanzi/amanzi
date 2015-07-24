/*
  This is the operators component of the Amanzi code.

  Copyright 2010-2013 held jointly by LANS/LANL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
           Ethan Coon (ecoon@lanl.gov)

  Base factory for diffusion operators.
*/

#include "errors.hh"

#include "BCs.hh"
#include "OperatorDefs.hh"
#include "OperatorDiffusionFactory.hh"
#include "OperatorDiffusionMFD.hh"
#include "OperatorDiffusionFV.hh"
#include "OperatorDiffusionMFDwithGravity.hh"
#include "OperatorDiffusionFVwithGravity.hh"

namespace Amanzi {
namespace Operators {

/* ******************************************************************
 * Initialization of the diffusion operators.
 ****************************************************************** */
Teuchos::RCP<OperatorDiffusion> OperatorDiffusionFactory::Create(
    Teuchos::RCP<const AmanziMesh::Mesh> mesh, 
    Teuchos::RCP<BCs> bc, 
    Teuchos::ParameterList& oplist,
    const AmanziGeometry::Point& g,
    int upwind_method)
{
  Teuchos::RCP<OperatorDiffusion> op = Teuchos::null;

  std::string name = oplist.get<std::string>("discretization primary");
  bool flag = oplist.get<bool>("gravity", false);

  // FV methods
  if (name == "fv: default" && !flag) {
    if (!oplist.isParameter("schema")) {
      Teuchos::Array<std::string> schema(1);
      schema[0] = "cell";
      oplist.set("schema",  schema);
    }
    op = Teuchos::rcp(new OperatorDiffusionFV(oplist, mesh));
    op->SetBCs(bc, bc);

  } else if (name == "fv: default" && flag) {
    if (!oplist.isParameter("schema")) {
      Teuchos::Array<std::string> schema(1);
      schema[0] = "cell";
      oplist.set("schema",  schema);
    }
    op = Teuchos::rcp(new OperatorDiffusionFVwithGravity(oplist, mesh));
    op->SetGravity(g);
    op->SetBCs(bc, bc);

  // MFD methods
  } else if (!flag) {
    if (!oplist.isParameter("schema")) {
      Teuchos::Array<std::string> schema(2);
      schema[0] = "cell"; schema[1] = "face";
      oplist.set("schema",  schema);
    }
    op = Teuchos::rcp(new OperatorDiffusionMFD(oplist, mesh));
    op->SetBCs(bc, bc);

  } else {
    if (!oplist.isParameter("schema")) {
      Teuchos::Array<std::string> schema(2);
      schema[0] = "cell"; schema[1] = "face";
      oplist.set("schema",  schema);
    }
    op = Teuchos::rcp(new OperatorDiffusionMFDwithGravity(oplist, mesh));
    op->SetGravity(g);
    op->SetBCs(bc, bc);
  }
  return op;
}


Teuchos::RCP<OperatorDiffusion> OperatorDiffusionFactory::Create(
    Teuchos::RCP<const AmanziMesh::Mesh> mesh,
    Teuchos::RCP<BCs> bc,
    Teuchos::ParameterList& oplist)
{
  Teuchos::RCP<OperatorDiffusion> op = Teuchos::null;

  std::string name = oplist.get<std::string>("discretization primary");

  if (name == "fv: default") {
    if (!oplist.isParameter("schema")) {
      Teuchos::Array<std::string> schema(1);
      schema[0] = "cell";
      oplist.set("schema",  schema);
    }
    op = Teuchos::rcp(new OperatorDiffusionFV(oplist, mesh));
    op->SetBCs(bc, bc);

  } else {
    if (!oplist.isParameter("schema")) {
      Teuchos::Array<std::string> schema(2);
      schema[0] = "cell"; schema[1] = "face";
      oplist.set("schema",  schema);
    }
    op = Teuchos::rcp(new OperatorDiffusionMFD(oplist, mesh));
    op->SetBCs(bc, bc);
  }
  return op;
}
  
Teuchos::RCP<OperatorDiffusion> OperatorDiffusionFactory::Create(Teuchos::ParameterList& oplist,
                                       Teuchos::RCP<const AmanziMesh::Mesh> mesh) {
  Teuchos::RCP<OperatorDiffusion> op = Teuchos::null;
  std::string name = oplist.get<std::string>("discretization primary");
  bool flag = oplist.get<bool>("gravity", false);
  
  // FV methods
  if (name == "fv: default" && !flag) {
    if (!oplist.isParameter("schema")) {
      Teuchos::Array<std::string> schema(1);
      schema[0] = "cell";
      oplist.set("schema",  schema);
    }
    op = Teuchos::rcp(new OperatorDiffusionFV(oplist, mesh));

  } else if (name == "fv: default" && flag) {
    if (!oplist.isParameter("schema")) {
      Teuchos::Array<std::string> schema(1);
      schema[0] = "cell";
      oplist.set("schema",  schema);
    }
    op = Teuchos::rcp(new OperatorDiffusionFVwithGravity(oplist, mesh));
    
  // MFD methods
  } else if (!flag) {
    op = Teuchos::rcp(new OperatorDiffusionMFD(oplist, mesh));
  } else {
    op = Teuchos::rcp(new OperatorDiffusionMFDwithGravity(oplist, mesh));
  }
  return op;
}
  
Teuchos::RCP<OperatorDiffusion> OperatorDiffusionFactory::Create(Teuchos::ParameterList& oplist,
                                       const Teuchos::RCP<Operator>& global_op) {
  Teuchos::RCP<OperatorDiffusion> op = Teuchos::null;
  std::string name = oplist.get<std::string>("discretization primary");
  bool flag = oplist.get<bool>("gravity", false);
  
  // FV methods
  if (name == "fv: default" && !flag) {
    op = Teuchos::rcp(new OperatorDiffusionFV(oplist, global_op));
  } else if (name == "fv: default" && flag) {
    op = Teuchos::rcp(new OperatorDiffusionFVwithGravity(oplist, global_op));
    
  // MFD methods
  } else if (!flag) {
    op = Teuchos::rcp(new OperatorDiffusionMFD(oplist, global_op));
  } else {
    op = Teuchos::rcp(new OperatorDiffusionMFDwithGravity(oplist, global_op));
  }
  return op;
}

  
}  // namespace Operators
}  // namespace Amanzi


