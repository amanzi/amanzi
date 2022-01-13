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
#include "Mesh.hh"
#include "MeshFactory.hh"
#include "State.hh"

// Flow
#include "Richards_PK.hh"

/* **************************************************************** */
TEST(FLOW_3D_RICHARDS) {
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Flow;

  Comm_ptr_type comm = Amanzi::getDefaultComm();
  int MyPID = comm->MyPID();
  if (MyPID == 0) std::cout << "Test: 3D Richards, crib model" << std::endl;

  /* read parameter list */
  std::string xmlFileName = "test/flow_richards_bc_cribs.xml";
  Teuchos::RCP<Teuchos::ParameterList> plist = Teuchos::getParametersFromXmlFile(xmlFileName);

  // create a mesh framework
  ParameterList regions_list = plist->get<Teuchos::ParameterList>("regions");
  Teuchos::RCP<Amanzi::AmanziGeometry::GeometricModel> gm =
      Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(3, regions_list, *comm));

  Preference pref;
  pref.clear();
  pref.push_back(Framework::MSTK);
  pref.push_back(Framework::STK);

  MeshFactory meshfactory(comm,gm);
  meshfactory.set_preference(pref);  
  // ParameterList mesh_list = plist->get<ParameterList>("mesh").get<ParameterList>("unstructured");
  // ParameterList factory_list = mesh_list.get<ParameterList>("generate mesh");
  // Teuchos::RCP<Mesh> mesh(factory(factory_list, gm));
  Teuchos::RCP<Mesh> mesh = meshfactory.create(0.,0.,0.,64.5,1.0,103.2,1,1,516);

  /* create a simple state and populate it */
  Amanzi::VerboseObject::global_hide_line_prefix = true;

  ParameterList state_list = plist->get<ParameterList>("state");
  RCP<State> S = rcp(new State(state_list));
  S->RegisterDomainMesh(rcp_const_cast<Mesh>(mesh));

  Teuchos::RCP<TreeVector> soln = Teuchos::rcp(new TreeVector());
  Richards_PK* RPK = new Richards_PK(plist, "flow", S, soln);

  RPK->Setup();
  S->Setup();
  S->InitializeFields();
  S->InitializeEvaluators();

  // initialize the Richards process kernel
  RPK->Initialize();
  S->CheckAllFieldsInitialized();

  RPK->CommitStep(0.0, 1.0, Tags::DEFAULT);  // dummay times

  // derive dependent variable
  const auto& p = *S->Get<CompositeVector>("pressure").ViewComponent("cell");
  const auto& ws = *S->Get<CompositeVector>("saturation_liquid").ViewComponent("cell");
  const auto& K = *S->Get<CompositeVector>("permeability").ViewComponent("cell");

  GMV::open_data_file(*mesh, (std::string)"flow.gmv");
  GMV::start_data();
  GMV::write_cell_data(p, 0, "pressure");
  GMV::write_cell_data(ws, 0, "saturation");
  GMV::write_cell_data(*K(0), 0, "permeability_x");
  GMV::write_cell_data(*K(1), 0, "permeability_y");
  GMV::write_cell_data(*K(2), 0, "permeability_z");
  GMV::close_data_file();

  /* check the pressure profile */
  int ncells = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  //for (int c = 0; c < ncells; c++) std::cout << (mesh->cell_centroid(c))[2] << " " << pressure[c] << std::endl;
  for (int c = 0; c < ncells; c++) CHECK(p[0][c] > 4500.0 && p[0][c] < 101325.0);

  delete RPK;
}
