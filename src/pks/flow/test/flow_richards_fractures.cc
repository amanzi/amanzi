/*
  Transport PK 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  License: BSD
  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <vector>

// TPLs
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "UnitTest++.h"

// Amanzi
#include "MeshFactory.hh"
#include "Mesh_MSTK.hh"
#include "State.hh"

// Flow
#include "Richards_PK.hh"

/* **************************************************************** */
TEST(RICHARDS_TWO_FRACTURES) {
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::Flow;
  using namespace Amanzi::AmanziGeometry;

  std::cout << "Test: Richards PK in two fractures" << std::endl;
  Comm_ptr_type comm = Amanzi::getDefaultComm();
  int MyPID = comm->MyPID();

  // read parameter list
  std::string xmlFileName = "test/flow_richards_fractures.xml";
  Teuchos::RCP<Teuchos::ParameterList> plist = Teuchos::getParametersFromXmlFile(xmlFileName);

  // create a mesh framework
  ParameterList region_list = plist->get<Teuchos::ParameterList>("regions");
  auto gm = Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(3, region_list, *comm));

  MeshFactory meshfactory(comm, gm);
  meshfactory.set_preference(Preference({Framework::MSTK}));
  RCP<const Mesh> mesh3D = meshfactory.create(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 10, 10, 10);

  // extract fractures mesh
  std::vector<std::string> setnames;
  setnames.push_back("fracture 1");
  setnames.push_back("fracture 2");

  RCP<const Mesh> mesh = meshfactory.create(mesh3D, setnames, AmanziMesh::FACE);

  int ncells_owned = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  int ncells_wghost = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::ALL);

  std::cout << "pid=" << MyPID << " cells: " << ncells_owned << " " << ncells_wghost << std::endl;

  Teuchos::ParameterList state_list = plist->sublist("state");
  RCP<State> S = rcp(new State(state_list));
  S->RegisterDomainMesh(rcp_const_cast<Mesh>(mesh));

  Teuchos::RCP<TreeVector> soln = Teuchos::rcp(new TreeVector());
  Teuchos::RCP<Richards_PK> RPK = Teuchos::rcp(new Richards_PK(plist, "flow", S, soln));
  RPK->Setup(S.ptr());

  S->Setup();
  S->InitializeFields();
  S->InitializeEvaluators();

  // initialize the Darcy process kernel
  RPK->Initialize(S.ptr());
  S->CheckAllFieldsInitialized();

  // transient solution
  double t_old(0.0), t_new, dt(0.5);
  for (int n = 0; n < 2; n++) {
    t_new = t_old + dt;

    RPK->AdvanceStep(t_old, t_new);
    RPK->CommitStep(t_old, t_new, S);

    t_old = t_new;
  }

  auto& p = *S->GetFieldData("pressure", "flow")->ViewComponent("cell");
  for (int c = 0; c < p.MyLength(); c++) {
    CHECK(p[0][c] > -1.0 && p[0][c] < 2.0);
  }
}





