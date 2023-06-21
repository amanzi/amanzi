/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Flow PK

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
#include "IO.hh"
#include "MeshAudit.hh"
#include "MeshFactory.hh"
#include "State.hh"

// Flow
#include "Richards_PK.hh"
#include "Richards_SteadyState.hh"

/* **************************************************************** */
TEST(FLOW_POROSITY_MODELS)
{
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Flow;

  Comm_ptr_type comm = Amanzi::getDefaultComm();
  int MyPID = comm->MyPID();
  if (MyPID == 0) std::cout << "Test: Richards, compressible porosity model" << std::endl;

  // read parameter list
  std::string xmlFileName = "test/flow_richards_porosity.xml";
  Teuchos::RCP<Teuchos::ParameterList> plist = Teuchos::getParametersFromXmlFile(xmlFileName);

  // create a mesh framework
  Teuchos::ParameterList regions_list = plist->get<Teuchos::ParameterList>("regions");
  Teuchos::RCP<Amanzi::AmanziGeometry::GeometricModel> gm =
    Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(2, regions_list, *comm));


  Preference pref;
  pref.clear();
  pref.push_back(Framework::MSTK);
  pref.push_back(Framework::STK);

  MeshFactory meshfactory(comm, gm);
  meshfactory.set_preference(pref);
  Teuchos::RCP<const Mesh> mesh = meshfactory.create(0.0, -2.0, 1.0, 0.0, 18, 18);

  // create a simple state and populate it
  Teuchos::ParameterList state_list = plist->sublist("state");
  Teuchos::RCP<State> S = Teuchos::rcp(new State(state_list));
  S->RegisterDomainMesh(Teuchos::rcp_const_cast<Mesh>(mesh));

  Teuchos::RCP<TreeVector> soln = Teuchos::rcp(new TreeVector());
  Teuchos::RCP<Richards_PK> RPK = Teuchos::rcp(new Richards_PK(plist, "flow", S, soln));

  RPK->Setup();
  S->Setup();
  S->InitializeFields();
  S->InitializeEvaluators();

  // modify the default state for the problem at hand
  std::string passwd("");
  auto& K = *S->GetW<CompositeVector>("permeability", "permeability").ViewComponent("cell");

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
  S->GetRecordW("permeability", "permeability").set_initialized();

  // -- fluid density and vicosity
  S->GetW<CompositeVector>("viscosity_liquid", "viscosity_liquid").PutScalar(1.0);

  // create the initial pressure function
  auto& p = *S->GetW<CompositeVector>("pressure", passwd).ViewComponent("cell");

  for (int c = 0; c < p.MyLength(); c++) {
    const Point& xc = mesh->cell_centroid(c);
    p[0][c] = xc[1] * (xc[1] + 2.0);
  }

  // initialize the Richards process kernel
  RPK->Initialize();
  S->CheckAllFieldsInitialized();

  // solve the problem
  TI_Specs ti_specs;
  ti_specs.T0 = 0.0;
  ti_specs.dT0 = 1.0;
  ti_specs.T1 = 100.0;
  ti_specs.max_itrs = 400;

  AdvanceToSteadyState(S, *RPK, ti_specs, soln);
  RPK->CommitStep(0.0, 1.0, Tags::DEFAULT); // dummy times

  // output
  double pmin, pmax;
  auto& phi = *S->GetW<CompositeVector>("porosity", "porosity").ViewComponent("cell");
  phi.MinValue(&pmin);
  phi.MaxValue(&pmax);
  CHECK(pmin + 0.02 < pmax);

  Teuchos::ParameterList& tmp =
    plist->sublist("PKs").sublist("flow").sublist("porosity models").sublist("POM for Material 2");
  std::cout << "Mat2: ref pressure:" << tmp.get<double>("reference pressure") << std::endl;
  std::cout << "      compressibility:" << tmp.get<double>("pore compressibility") << std::endl;
  std::cout << "Porosity min = " << pmin << std::endl;
  std::cout << "Porosity max = " << pmax << std::endl;

  // WriteStateStatistics(*S);

  if (MyPID == 0) {
    GMV::open_data_file(*mesh, (std::string) "flow.gmv");
    GMV::start_data();
    GMV::write_cell_data(p, 0, "pressure");
    GMV::close_data_file();
  }

  // check the pressure
  int ncells = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  for (int c = 0; c < ncells; c++) CHECK(p[0][c] > -4.0 && p[0][c] < 0.01);
}
