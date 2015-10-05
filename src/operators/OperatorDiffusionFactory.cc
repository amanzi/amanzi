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
* Initialization of diffusion operator: method 1 with optional gravity
****************************************************************** */
Teuchos::RCP<OperatorDiffusion> OperatorDiffusionFactory::Create(
    Teuchos::RCP<const AmanziMesh::Mesh> mesh, 
    Teuchos::RCP<BCs> bc, 
    Teuchos::ParameterList& oplist,
    double rho, const AmanziGeometry::Point& g)
{
  Teuchos::RCP<OperatorDiffusion> op = Teuchos::null;

  std::string name = oplist.get<std::string>("discretization primary");
  bool flag = oplist.get<bool>("gravity", false);

  // FV methods
  if (name == "fv: default" && !flag) {
    SetCellSchema_(oplist);
    op = Teuchos::rcp(new OperatorDiffusionFV(oplist, mesh));
    op->SetBCs(bc, bc);

  } else if (name == "fv: default" && flag) {
    SetCellSchema_(oplist);
    op = Teuchos::rcp(new OperatorDiffusionFVwithGravity(oplist, mesh, rho, g));
    op->SetBCs(bc, bc);

  // MFD methods
  } else if (!flag) {
    SetCellFaceSchema_(oplist);
    op = Teuchos::rcp(new OperatorDiffusionMFD(oplist, mesh));
    op->SetBCs(bc, bc);

  } else {
    SetCellFaceSchema_(oplist);
    op = Teuchos::rcp(new OperatorDiffusionMFDwithGravity(oplist, mesh, rho, g));
    op->SetBCs(bc, bc);
  }
  return op;
}


/* ******************************************************************
* Initialization of diffusion operator: method 2 with optional gravity
****************************************************************** */
Teuchos::RCP<OperatorDiffusion> OperatorDiffusionFactory::Create(
    Teuchos::RCP<const AmanziMesh::Mesh> mesh, 
    Teuchos::RCP<BCs> bc, 
    Teuchos::ParameterList& oplist,
    Teuchos::RCP<const CompositeVector> rho,
    const AmanziGeometry::Point& g)
{
  Teuchos::RCP<OperatorDiffusion> op = Teuchos::null;

  std::string name = oplist.get<std::string>("discretization primary");
  bool flag = oplist.get<bool>("gravity", false);

  // FV methods
  if (name == "fv: default" && !flag) {
    SetCellSchema_(oplist);
    op = Teuchos::rcp(new OperatorDiffusionFV(oplist, mesh));
    op->SetBCs(bc, bc);

  } else if (name == "fv: default" && flag) {
    SetCellSchema_(oplist);
    op = Teuchos::rcp(new OperatorDiffusionFVwithGravity(oplist, mesh, rho, g));
    op->SetBCs(bc, bc);

  // MFD methods
  } else if (!flag) {
    SetCellFaceSchema_(oplist);
    op = Teuchos::rcp(new OperatorDiffusionMFD(oplist, mesh));
    op->SetBCs(bc, bc);

  } else {
    SetCellFaceSchema_(oplist);
    op = Teuchos::rcp(new OperatorDiffusionMFDwithGravity(oplist, mesh, rho, g));
    op->SetBCs(bc, bc);
  }
  return op;
}


/* ******************************************************************
* Initialization of diffusion operator: method 1 without gravity.
****************************************************************** */
Teuchos::RCP<OperatorDiffusion> OperatorDiffusionFactory::Create(
    Teuchos::RCP<const AmanziMesh::Mesh> mesh,
    Teuchos::RCP<BCs> bc,
    Teuchos::ParameterList& oplist)
{
  Teuchos::RCP<OperatorDiffusion> op = Teuchos::null;
  std::string name = oplist.get<std::string>("discretization primary");
  
  // FV methods
  if (name == "fv: default") {
    SetCellSchema_(oplist);
    op = Teuchos::rcp(new OperatorDiffusionFV(oplist, mesh));
    op->SetBCs(bc, bc);

  // MFD methods
  } else {
    SetCellFaceSchema_(oplist);
    op = Teuchos::rcp(new OperatorDiffusionMFD(oplist, mesh));
    op->SetBCs(bc, bc);
  }
  return op;
}
  

/* ******************************************************************
* Initialization of diffusion operator: method 2 without gravity.
****************************************************************** */
Teuchos::RCP<OperatorDiffusion> OperatorDiffusionFactory::Create(
    Teuchos::ParameterList& oplist,
    Teuchos::RCP<const AmanziMesh::Mesh> mesh)
{
  Teuchos::RCP<OperatorDiffusion> op = Teuchos::null;
  std::string name = oplist.get<std::string>("discretization primary");
  
  if (name == "fv: default") {
    SetCellSchema_(oplist);
    op = Teuchos::rcp(new OperatorDiffusionFV(oplist, mesh));
  } else {
    op = Teuchos::rcp(new OperatorDiffusionMFD(oplist, mesh));
  }
  return op;
}
  

/* ******************************************************************
* Initialization of diffusion operator: method 3 without gravity.
****************************************************************** */
Teuchos::RCP<OperatorDiffusion> OperatorDiffusionFactory::Create(
    Teuchos::ParameterList& oplist,
    const Teuchos::RCP<Operator>& global_op)
{
  Teuchos::RCP<OperatorDiffusion> op = Teuchos::null;
  std::string name = oplist.get<std::string>("discretization primary");
  
  if (name == "fv: default") {
    op = Teuchos::rcp(new OperatorDiffusionFV(oplist, global_op));
  } else {
    op = Teuchos::rcp(new OperatorDiffusionMFD(oplist, global_op));
  }
  return op;
}


/* ******************************************************************
* Add missing schemas to the operator parameter list.
****************************************************************** */
void OperatorDiffusionFactory::SetCellSchema_(Teuchos::ParameterList& oplist) 
{
  if (!oplist.isParameter("schema")) {
    Teuchos::Array<std::string> schema(1);
    schema[0] = "cell";
    oplist.set("schema",  schema);
  }
}


void OperatorDiffusionFactory::SetCellFaceSchema_(Teuchos::ParameterList& oplist) 
{
  if (!oplist.isParameter("schema")) {
    Teuchos::Array<std::string> schema(2);
    schema[0] = "cell";
    schema[1] = "face";
    oplist.set("schema",  schema);
  }
}

}  // namespace Operators
}  // namespace Amanzi


