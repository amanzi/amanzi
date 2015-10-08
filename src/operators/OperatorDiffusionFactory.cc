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
*
* This is the constructor used by Amanzi
****************************************************************** */
Teuchos::RCP<OperatorDiffusion> OperatorDiffusionFactory::Create(
    Teuchos::ParameterList& oplist,
    const Teuchos::RCP<const AmanziMesh::Mesh>& mesh, 
    const Teuchos::RCP<BCs>& bc, 
    double rho,
    const AmanziGeometry::Point& g)
{
  Teuchos::RCP<OperatorDiffusion> op = Teuchos::null;

  std::string name = oplist.get<std::string>("discretization primary");
  bool flag = oplist.get<bool>("gravity", false);

  // FV methods
  if (name == "fv: default" && !flag) {
    op = Teuchos::rcp(new OperatorDiffusionFV(oplist, mesh));
    op->SetBCs(bc, bc);

  } else if (name == "fv: default" && flag) {
    op = Teuchos::rcp(new OperatorDiffusionFVwithGravity(oplist, mesh, rho, g));
    op->SetBCs(bc, bc);

  // MFD methods
  } else if (!flag) {
    op = Teuchos::rcp(new OperatorDiffusionMFD(oplist, mesh));
    op->SetBCs(bc, bc);

  } else {
    op = Teuchos::rcp(new OperatorDiffusionMFDwithGravity(oplist, mesh, rho, g));
    op->SetBCs(bc, bc);
  }
  return op;
}


/* ******************************************************************
* Initialization of diffusion operator: method 2 with optional gravity
*
* This is the factory used by Amanzi, though it makes life difficult
* for time-varying density.
****************************************************************** */
Teuchos::RCP<OperatorDiffusion> OperatorDiffusionFactory::Create(
    Teuchos::ParameterList& oplist,
    const Teuchos::RCP<const AmanziMesh::Mesh>& mesh, 
    const Teuchos::RCP<BCs>& bc, 
    const Teuchos::RCP<const CompositeVector>& rho,
    const AmanziGeometry::Point& g)
{
  Teuchos::RCP<OperatorDiffusion> op = Teuchos::null;

  std::string name = oplist.get<std::string>("discretization primary");
  bool flag = oplist.get<bool>("gravity", false);

  // FV methods
  if (name == "fv: default" && !flag) {
    op = Teuchos::rcp(new OperatorDiffusionFV(oplist, mesh));
    op->SetBCs(bc, bc);

  } else if (name == "fv: default" && flag) {
    Teuchos::RCP<OperatorDiffusionFVwithGravity> op_g =
      Teuchos::rcp(new OperatorDiffusionFVwithGravity(oplist, mesh, g));
    op_g->SetBCs(bc, bc);
    op_g->SetDensity(rho);
    op = op_g;

  // MFD methods
  } else if (!flag) {
    op = Teuchos::rcp(new OperatorDiffusionMFD(oplist, mesh));
    op->SetBCs(bc, bc);

  } else {
    Teuchos::RCP<OperatorDiffusionMFDwithGravity> op_g =
      Teuchos::rcp(new OperatorDiffusionMFDwithGravity(oplist, mesh, g));
    op_g->SetBCs(bc, bc);
    op_g->SetDensity(rho);
    op = op_g;
  }
  return op;
}


/* ******************************************************************
* Initialization of diffusion operator: method 1 without gravity.
*
* No gravity included, straight diffusion operator.
****************************************************************** */
Teuchos::RCP<OperatorDiffusion> OperatorDiffusionFactory::Create(
    Teuchos::ParameterList& oplist,
    const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
    const Teuchos::RCP<BCs>& bc)
{
  Teuchos::RCP<OperatorDiffusion> op = Teuchos::null;
  std::string name = oplist.get<std::string>("discretization primary");
  
  // FV methods
  if (name == "fv: default") {
    op = Teuchos::rcp(new OperatorDiffusionFV(oplist, mesh));
    op->SetBCs(bc, bc);

  // MFD methods
  } else {
    op = Teuchos::rcp(new OperatorDiffusionMFD(oplist, mesh));
    op->SetBCs(bc, bc);
  }
  return op;
}
  

/* ******************************************************************
* Initialization of diffusion operator: method 2 without gravity.
*
* No gravity included, straight diffusion operator.
****************************************************************** */
Teuchos::RCP<OperatorDiffusion> OperatorDiffusionFactory::Create(
    Teuchos::ParameterList& oplist,
    const Teuchos::RCP<const AmanziMesh::Mesh>& mesh)
{
  Teuchos::RCP<OperatorDiffusion> op = Teuchos::null;
  std::string name = oplist.get<std::string>("discretization primary");
  
  if (name == "fv: default") {
    op = Teuchos::rcp(new OperatorDiffusionFV(oplist, mesh));
  } else {
    op = Teuchos::rcp(new OperatorDiffusionMFD(oplist, mesh));
  }
  return op;
}
  

/* ******************************************************************
* Initialization of diffusion operator: method 3 without gravity.
*
* No gravity included, straight diffusion operator.
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
* Initialization of diffusion operator with gravity: method 1
*
* With gravity, assumed vector, temporally varying density.
* Used by ATS.
****************************************************************** */
Teuchos::RCP<OperatorDiffusionWithGravity>
OperatorDiffusionFactory::CreateWithGravity(
    Teuchos::ParameterList& oplist,
    const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
    const Teuchos::RCP<BCs>& bc)
{
  Teuchos::RCP<OperatorDiffusionWithGravity> op = Teuchos::null;
  std::string name = oplist.get<std::string>("discretization primary");
  
  if (name == "fv: default") {
    op = Teuchos::rcp(new OperatorDiffusionFVwithGravity(oplist, mesh));
  } else {
    op = Teuchos::rcp(new OperatorDiffusionMFDwithGravity(oplist, mesh));
  }
  op->SetBCs(bc, bc);
  return op;
}

/* ******************************************************************
* Initialization of diffusion operator with gravity: method 2
*
* With gravity, assumed vector, temporally varying density.
* Used by ATS.
****************************************************************** */
Teuchos::RCP<OperatorDiffusionWithGravity>
OperatorDiffusionFactory::CreateWithGravity(
    Teuchos::ParameterList& oplist,
    const Teuchos::RCP<Operator>& global_op,
    const Teuchos::RCP<BCs>& bc)
{
  Teuchos::RCP<OperatorDiffusionWithGravity> op = Teuchos::null;
  std::string name = oplist.get<std::string>("discretization primary");
  
  if (name == "fv: default") {
    op = Teuchos::rcp(new OperatorDiffusionFVwithGravity(oplist, global_op));
  } else {
    op = Teuchos::rcp(new OperatorDiffusionMFDwithGravity(oplist, global_op));
  }
  op->SetBCs(bc, bc);
  return op;
}


/* ******************************************************************
* Initialization of diffusion operator with gravity: method 1
*
* With gravity, assumed vector, temporally varying density.
* Used by ATS.
****************************************************************** */
Teuchos::RCP<OperatorDiffusionWithGravity>
OperatorDiffusionFactory::CreateWithGravity(
    Teuchos::ParameterList& oplist,
    const Teuchos::RCP<const AmanziMesh::Mesh>& mesh)
{
  Teuchos::RCP<OperatorDiffusionWithGravity> op = Teuchos::null;
  std::string name = oplist.get<std::string>("discretization primary");
  
  if (name == "fv: default") {
    op = Teuchos::rcp(new OperatorDiffusionFVwithGravity(oplist, mesh));
  } else {
    op = Teuchos::rcp(new OperatorDiffusionMFDwithGravity(oplist, mesh));
  }
  return op;
}

/* ******************************************************************
* Initialization of diffusion operator with gravity: method 2
*
* With gravity, assumed vector, temporally varying density.
* Used by ATS.
****************************************************************** */
Teuchos::RCP<OperatorDiffusionWithGravity>
OperatorDiffusionFactory::CreateWithGravity(
    Teuchos::ParameterList& oplist,
    const Teuchos::RCP<Operator>& global_op)
{
  Teuchos::RCP<OperatorDiffusionWithGravity> op = Teuchos::null;
  std::string name = oplist.get<std::string>("discretization primary");
  
  if (name == "fv: default") {
    op = Teuchos::rcp(new OperatorDiffusionFVwithGravity(oplist, global_op));
  } else {
    op = Teuchos::rcp(new OperatorDiffusionMFDwithGravity(oplist, global_op));
  }
  return op;
}
  

}  // namespace Operators
}  // namespace Amanzi


