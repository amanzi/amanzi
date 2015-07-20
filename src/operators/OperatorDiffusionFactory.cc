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
  std::string name = oplist.get<std::string>("discretization primary");
  bool flag = oplist.get<bool>("gravity", false);

  // FV methods
  if (name == "fv: default" && !flag) {
    Teuchos::RCP<OperatorDiffusionFV> op = Teuchos::rcp(new OperatorDiffusionFV(oplist, mesh));
    op->SetBCs(bc, bc);
    return op;
  } else if (name == "fv: default" && flag) {
    Teuchos::RCP<OperatorDiffusionFVwithGravity> op = Teuchos::rcp(new OperatorDiffusionFVwithGravity(oplist, mesh));
    op->SetGravity(g);
    op->SetBCs(bc, bc);
    return op;

  // MFD methods
  } else if (!flag) {
    Teuchos::RCP<OperatorDiffusionMFD> op = Teuchos::rcp(new OperatorDiffusionMFD(oplist, mesh));
    op->SetBCs(bc, bc);
    return op;
  } else {
    Teuchos::RCP<OperatorDiffusionMFDwithGravity> op = Teuchos::rcp(new OperatorDiffusionMFDwithGravity(oplist, mesh));
    op->SetGravity(g);
    op->SetBCs(bc, bc);
    return op;
  }
}


Teuchos::RCP<OperatorDiffusion> OperatorDiffusionFactory::Create(
    Teuchos::RCP<const AmanziMesh::Mesh> mesh,
    Teuchos::RCP<BCs> bc,
    Teuchos::ParameterList& oplist)
{
  std::string name = oplist.get<std::string>("discretization primary");

  if (name == "fv: default") {
    Teuchos::RCP<OperatorDiffusionFV> op = Teuchos::rcp(new OperatorDiffusionFV(oplist, mesh));
    op->SetBCs(bc, bc);
    return op;
  } else {
    Teuchos::RCP<OperatorDiffusionMFD> op = Teuchos::rcp(new OperatorDiffusionMFD(oplist, mesh));
    op->SetBCs(bc, bc);
    return op;
  }
}
  
Teuchos::RCP<OperatorDiffusion> OperatorDiffusionFactory::Create(Teuchos::ParameterList& oplist,
                                       Teuchos::RCP<const AmanziMesh::Mesh> mesh) {
  std::string name = oplist.get<std::string>("discretization primary");
  bool flag = oplist.get<bool>("gravity", false);
  
  // FV methods
  if (name == "fv: default" && !flag) {
    Teuchos::RCP<OperatorDiffusionFV> op = Teuchos::rcp(new OperatorDiffusionFV(oplist, mesh));
    return op;
  } else if (name == "fv: default" && flag) {
    Teuchos::RCP<OperatorDiffusionFVwithGravity> op = Teuchos::rcp(new OperatorDiffusionFVwithGravity(oplist, mesh));
    return op;
    
  // MFD methods
  } else if (!flag) {
    Teuchos::RCP<OperatorDiffusionMFD> op = Teuchos::rcp(new OperatorDiffusionMFD(oplist, mesh));
    return op;
  } else {
    Teuchos::RCP<OperatorDiffusionMFDwithGravity> op = Teuchos::rcp(new OperatorDiffusionMFDwithGravity(oplist, mesh));
    return op;
  }
}
  
Teuchos::RCP<OperatorDiffusion> OperatorDiffusionFactory::Create(Teuchos::ParameterList& oplist,
                                       const Teuchos::RCP<Operator>& global_op) {
  std::string name = oplist.get<std::string>("discretization primary");
  bool flag = oplist.get<bool>("gravity", false);
  
  // FV methods
  if (name == "fv: default" && !flag) {
    Teuchos::RCP<OperatorDiffusionFV> op = Teuchos::rcp(new OperatorDiffusionFV(oplist, global_op));
    return op;
  } else if (name == "fv: default" && flag) {
    Teuchos::RCP<OperatorDiffusionFVwithGravity> op = Teuchos::rcp(new OperatorDiffusionFVwithGravity(oplist, global_op));
    return op;
    
  // MFD methods
  } else if (!flag) {
    Teuchos::RCP<OperatorDiffusionMFD> op = Teuchos::rcp(new OperatorDiffusionMFD(oplist, global_op));
    return op;
  } else {
    Teuchos::RCP<OperatorDiffusionMFDwithGravity> op = Teuchos::rcp(new OperatorDiffusionMFDwithGravity(oplist, global_op));
    return op;
  }
}

  
}  // namespace Operators
}  // namespace Amanzi


