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
#include "GMVMesh.hh"
#include "MeshFactory.hh"
#include "Mesh_MSTK.hh"
#include "State.hh"

// Flow
#include "Darcy_PK.hh"

/* **************************************************************** */
TEST(DARCY_TWO_FRACTURES) {
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::Flow;
  using namespace Amanzi::AmanziGeometry;

  std::cout << "Test: Darcy PK in two fractures" << std::endl;
  Comm_ptr_type comm = Amanzi::getDefaultComm();
  int MyPID = comm->MyPID();

  // read parameter list
  std::string xmlFileName = "test/flow_darcy_fractures.xml";
  Teuchos::RCP<Teuchos::ParameterList> plist = Teuchos::getParametersFromXmlFile(xmlFileName);

  // create a mesh framework
  ParameterList region_list = plist->get<Teuchos::ParameterList>("regions");
  auto gm = Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(3, region_list, *comm));

  MeshFactory meshfactory(comm, gm);
  meshfactory.set_preference(Preference({Framework::MSTK, Framework::STK}));
  RCP<const Mesh> mesh3D = meshfactory.create(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 10, 10, 10);

  // extract fractures mesh
  std::vector<std::string> setnames;
  setnames.push_back("fracture 1");
  setnames.push_back("fracture 2");

  RCP<Mesh> mesh = meshfactory.create(mesh3D, setnames, AmanziMesh::FACE);

  int ncells_owned = 0;
  int ncells_wghost = 0;
  if (mesh.get()) {
    ncells_owned = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
    ncells_wghost = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::ALL);
  }

  std::cout << "pid=" << MyPID << " cells: " << ncells_owned << " " << ncells_wghost << std::endl;

  Teuchos::ParameterList state_list = plist->sublist("state");
  RCP<State> S = rcp(new State(state_list));
  S->RegisterDomainMesh(mesh);

  Teuchos::RCP<TreeVector> soln = Teuchos::rcp(new TreeVector());
  Teuchos::RCP<Darcy_PK> DPK = Teuchos::rcp(new Darcy_PK(plist, "flow", S, soln));
  DPK->Setup(S.ptr());

  S->Setup();
  S->InitializeFields();
  S->InitializeEvaluators();

  // modify the default state
  // -- storativity
  S->GetW<CompositeVector>("specific_storage", Tags::DEFAULT, "flow").PutScalar(2.0);
  S->GetRecordW("specific_storage", Tags::DEFAULT, "flow").set_initialized();

  // create the initial pressure function
  auto& p = *S->GetW<CompositeVector>("pressure", Tags::DEFAULT, "flow").ViewComponent("cell");

  for (int c = 0; c < p.MyLength(); c++) {
    const Point& xc = mesh->cell_centroid(c);
    p[0][c] = xc[0] * (xc[1] + 2.0);
  }
  S->GetRecordW("pressure", "flow").set_initialized();

  // initialize the Darcy process kernel
  DPK->Initialize(S.ptr());
  S->CheckAllFieldsInitialized();

  // transient solution
  double t_old(0.0), t_new, dt(0.5);
  for (int n = 0; n < 2; n++) {
    t_new = t_old + dt;

    DPK->AdvanceStep(t_old, t_new);
    DPK->CommitStep(t_old, t_new, S);

    t_old = t_new;
  }

  for (int c = 0; c < p.MyLength(); c++) {
    CHECK(p[0][c] > -1.0 && p[0][c] < 2.0);
  }

  if (MyPID == 0) {
    GMV::open_data_file(*mesh, (std::string)"flow.gmv");
    GMV::start_data();
    GMV::write_cell_data(p, 0, "pressure");
    GMV::close_data_file();
  }

  
}





