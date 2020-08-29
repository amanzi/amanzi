/*
  Flow PK

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

// TPLs
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "UnitTest++.h"

// Amanzi
#include "GMVMesh.hh"
#include "MeshFactory.hh"
#include "State.hh"

// Flow
#include "Darcy_PK.hh"

/* *********************************************************************
* 2D transient test.
********************************************************************* */
TEST(FLOW_2D_TRANSIENT_DARCY) {
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Flow;

  Comm_ptr_type comm = Amanzi::getDefaultComm();
  int MyPID = comm->MyPID();
  if (MyPID == 0) std::cout << "Test: 2D transient Darcy, 2-layer model" << std::endl;

  // read parameter list
  std::string xmlFileName = "test/flow_darcy_transient_2D.xml";
  Teuchos::RCP<ParameterList> plist = Teuchos::getParametersFromXmlFile(xmlFileName);

  // create SIMPLE mesh framework
  ParameterList regions_list = plist->get<Teuchos::ParameterList>("regions");
  Teuchos::RCP<Amanzi::AmanziGeometry::GeometricModel> gm =
      Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(2, regions_list, *comm));

  Preference pref;
  pref.clear();
  pref.push_back(Framework::MSTK);
  pref.push_back(Framework::STK);

  MeshFactory meshfactory(comm,gm);
  meshfactory.set_preference(pref);
  RCP<const Mesh> mesh = meshfactory.create(0.0, -2.0, 1.0, 0.0, 18, 18);

  // create a simple state and populate it
  Amanzi::VerboseObject::global_hide_line_prefix = true;
  Amanzi::VerboseObject::global_default_level = Teuchos::VERB_EXTREME;

  Teuchos::ParameterList state_list = plist->sublist("state");
  RCP<State> S = rcp(new State(state_list));
  S->RegisterDomainMesh(rcp_const_cast<Mesh>(mesh));

  Teuchos::RCP<TreeVector> soln = Teuchos::rcp(new TreeVector());
  Teuchos::RCP<Darcy_PK> DPK = Teuchos::rcp(new Darcy_PK(plist, "flow", S, soln));
  DPK->Setup(S.ptr());
  std::cout << "Owner of " << S->GetField("permeability")->fieldname() 
            << " is " << S->GetField("permeability")->owner() << "\n";

  S->Setup();
  S->InitializeFields();
  S->InitializeEvaluators();

  // modify the default state for the problem at hand
  // -- permeability
  std::string passwd("flow"); 
  Epetra_MultiVector& K = *S->GetFieldData("permeability", passwd)->ViewComponent("cell", false);

  AmanziMesh::Entity_ID_List block;

  mesh->get_set_entities("Material 1", AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED, &block);
  for (int i = 0; i != block.size(); ++i) {
    int c = block[i];
    K[0][c] = 0.1;
    K[1][c] = 2.0;
  }

  mesh->get_set_entities("Material 2", AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED, &block);
  for (int i = 0; i != block.size(); ++i) {
    int c = block[i];
    K[0][c] = 0.5;
    K[1][c] = 0.5;
  } 
  S->GetField("permeability", "flow")->set_initialized();

  // -- fluid density and viscosity
  *S->GetScalarData("const_fluid_density", passwd) = 1.0;
  S->GetField("const_fluid_density", "flow")->set_initialized();

  *S->GetScalarData("const_fluid_viscosity", passwd) = 1.0;
  S->GetField("const_fluid_viscosity", "flow")->set_initialized();

  // -- gravity
  Epetra_Vector& gravity = *S->GetConstantVectorData("gravity", "state");
  gravity[1] = -1.0;
  S->GetField("gravity", "state")->set_initialized();

  // -- storativity
  S->GetFieldData("specific_storage", passwd)->PutScalar(2.0);
  S->GetField("specific_storage", "flow")->set_initialized();

  // create the initial pressure function
  Epetra_MultiVector& p = *S->GetFieldData("pressure", passwd)->ViewComponent("cell", false);

  for (int c = 0; c < p.MyLength(); c++) {
    const Point& xc = mesh->cell_centroid(c);
    p[0][c] = xc[1] * (xc[1] + 2.0);
  }
  S->GetField("pressure", "flow")->set_initialized();

  // initialize the Darcy process kernel
  DPK->Initialize(S.ptr());
  S->CheckAllFieldsInitialized();

  // transient solution
  double t_old(0.0), t_new, dt(0.1);
  for (int n = 0; n < 10; n++) {
    t_new = t_old + dt;

    DPK->AdvanceStep(t_old, t_new);
    DPK->CommitStep(t_old, t_new, S);

    t_old = t_new;

    if (MyPID == 0 && n > 5) {
      GMV::open_data_file(*mesh, (std::string)"flow.gmv");
      GMV::start_data();
      GMV::write_cell_data(p, 0, "pressure");
      GMV::close_data_file();
    }
  }
}


/* **************************************************************** */
TEST(FLOW_3D_TRANSIENT_DARCY) {
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Flow;

  Comm_ptr_type comm = Amanzi::getDefaultComm();
  int MyPID = comm->MyPID();
  if (MyPID == 0) std::cout << "\nTest: 3D transient Darcy, 3-layer model" << std::endl;

  // read parameter list
  std::string xmlFileName = "test/flow_darcy_transient_3D.xml";
  Teuchos::RCP<ParameterList> plist = Teuchos::getParametersFromXmlFile(xmlFileName);

  /* create an SIMPLE mesh framework */
  ParameterList regions_list = plist->get<Teuchos::ParameterList>("regions");
  Teuchos::RCP<Amanzi::AmanziGeometry::GeometricModel> gm =
      Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(3, regions_list, *comm));


  Preference pref;
  pref.clear();
  pref.push_back(Framework::MSTK);

  MeshFactory meshfactory(comm,gm);
  meshfactory.set_preference(pref);
  RCP<const Mesh> mesh = meshfactory.create(0.0, 0.0, -2.0, 1.0, 2.0, 0.0, 18, 18, 18);

  /* create and populate flow state */
  Amanzi::VerboseObject::global_hide_line_prefix = true;

  Teuchos::ParameterList state_list = plist->sublist("state");
  RCP<State> S = rcp(new State(state_list));
  S->RegisterDomainMesh(rcp_const_cast<Mesh>(mesh));

  Teuchos::RCP<TreeVector> soln = Teuchos::rcp(new TreeVector());
  Teuchos::RCP<Darcy_PK> DPK = Teuchos::rcp(new Darcy_PK(plist, "flow", S, soln));
  DPK->Setup(S.ptr());
  S->Setup();
  S->InitializeFields();
  S->InitializeEvaluators();

  /* modify the default state for the problem at hand */
  std::string passwd("flow"); 
  Epetra_MultiVector& K = *S->GetFieldData("permeability", passwd)->ViewComponent("cell", false);
  
  AmanziMesh::Entity_ID_List block;
  mesh->get_set_entities("Material 1", AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED, &block);
  for (int i = 0; i != block.size(); ++i) {
    int c = block[i];
    K[0][c] = 0.1;
    K[1][c] = 0.1;
    K[2][c] = 2.0;
  }

  mesh->get_set_entities("Material 2", AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED, &block);
  for (int i = 0; i != block.size(); ++i) {
    int c = block[i];
    K[0][c] = 0.5;
    K[1][c] = 0.5;
    K[2][c] = 0.5;
  }

  *S->GetScalarData("const_fluid_density", passwd) = 1.0;
  *S->GetScalarData("const_fluid_viscosity", passwd) = 1.0;
  Epetra_Vector& gravity = *S->GetConstantVectorData("gravity", "state");
  gravity[2] = -1.0;

  S->GetFieldData("specific_storage", passwd)->PutScalar(1.0);
  S->GetFieldData("specific_yield", passwd)->PutScalar(0.0);

  /* create the initial pressure function */
  Epetra_MultiVector& p = *S->GetFieldData("pressure", passwd)->ViewComponent("cell", false);

  for (int c = 0; c < p.MyLength(); c++) {
    const Point& xc = mesh->cell_centroid(c);
    p[0][c] = xc[2] * (xc[2] + 2.0) * (xc[1] + 1.0) * (xc[0] - 1.0);
  }

  // initialize the Darcy process kernel
  DPK->Initialize(S.ptr());
  S->CheckAllFieldsInitialized();

  // transient solution
  double t_old(0.0), t_new, dt(0.1);
  for (int n = 0; n < 5; n++) {
    t_new = t_old + dt;

    DPK->AdvanceStep(t_old, t_new);
    DPK->CommitStep(t_old, t_new, S);

    if (MyPID == 0) {
      GMV::open_data_file(*mesh, (std::string)"flow.gmv");
      GMV::start_data();
      GMV::write_cell_data(p, 0, "pressure");
      GMV::close_data_file();
    }
  }
}
