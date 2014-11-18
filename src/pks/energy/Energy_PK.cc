/*
  This is the energy component of the Amanzi code. 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <string>
#include <vector>

#include "Teuchos_ParameterList.hpp"

#include "GMVMesh.hh"
#include "Mesh.hh"
#include "mfd3d.hh"
#include "State.hh"

#include "Energy_PK.hh"

namespace Amanzi {
namespace Energy {

/* ******************************************************************
* Default constructor for Energy PK.
****************************************************************** */
Energy_PK::Energy_PK(Teuchos::ParameterList& glist, Teuchos::RCP<State> S)
    : vo_(NULL), passwd_("state")
{
  S_ = S;
  mesh_ = S->GetMesh();
  dim = mesh_->space_dimension();

  plist_ = glist.sublist("Energy");

  ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  ncells_wghost = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::USED);

  nfaces_owned = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::OWNED);
  nfaces_wghost = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::USED);

  Initialize();
}


/* ******************************************************************
* Construction of PK global variables.
****************************************************************** */
void Energy_PK::Initialize()
{
  // require state variables
  if (!S_->HasField("porosity")) {
    S_->RequireField("porosity", passwd_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
  }
  if (!S_->HasField("fluid_density")) {
    S_->RequireScalar("fluid_density", passwd_);
  }
  if (!S_->HasField("darcy_flux")) {
    S_->RequireField("darcy_flux", passwd_)->SetMesh(mesh_)->SetGhosted(true)
       ->SetComponent("face", AmanziMesh::FACE, 1);

    // Teuchos::ParameterList elist;
    // elist.set<std::string>("evaluator name", "darcy_flux");
    // darcy_flux_eval = Teuchos::rcp(new PrimaryVariableFieldEvaluator(elist));
    // S->SetFieldEvaluator("darcy_flux", darcy_flux_eval);
  }

  // for creating fields
  std::vector<std::string> names(2);
  names[0] = "cell";
  names[1] = "face";
 
  std::vector<AmanziMesh::Entity_kind> locations(2);
  locations[0] = AmanziMesh::CELL;
  locations[1] = AmanziMesh::FACE;
 
  std::vector<int> ndofs(2, 1);
  
  if (!S_->HasField("temperature")) {
    S_->RequireField("temperature", passwd_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponents(names, locations, ndofs);
  }

  // create verbosity object
  vo_ = new VerboseObject("EnergyPK", plist_);

  // Process Native XML.
  // ProcessParameterList(plist_);
  K.resize(ncells_wghost);
  for (int c = 0; c < ncells_wghost; c++) {
    K[c].Init(dim, 1);
    K[c](0, 0) = 1.0;
  }

  // Select a proper matrix class. 
  Teuchos::ParameterList& tmp_list = plist_.sublist("operators").sublist("diffusion operator");
  Teuchos::ParameterList oplist_matrix = tmp_list.sublist("matrix");
  Teuchos::ParameterList oplist_pc = tmp_list.sublist("preconditioner");

  op_bc_ = Teuchos::rcp(new Operators:: BCs(Operators::OPERATOR_BC_TYPE_FACE, bc_model_, bc_value_));
  AmanziGeometry::Point g(dim);

  Operators::OperatorDiffusionFactory opfactory;
  op_matrix_ = opfactory.Create(mesh_, op_bc_, oplist_matrix, g, 0);

  op_matrix_->Init();
  op_matrix_->InitOperator(K, Teuchos::null, Teuchos::null, 1.0, 1.0);
}

}  // namespace Energy
}  // namespace Amanzi

