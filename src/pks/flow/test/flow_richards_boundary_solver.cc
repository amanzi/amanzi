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
#include "MeshAudit.hh"
#include "MeshFactory.hh"
#include "State.hh"

// Flow
#include "Richards_PK.hh"
#include "Richards_SteadyState.hh"

/* *******************************************************************
* A
******************************************************************* */
TEST(FLOW_BOUNDARY_SOLVER)
{
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Flow;

  Comm_ptr_type comm = Amanzi::getDefaultComm();
  int MyPID = comm->MyPID();
  if (MyPID == 0) std::cout << "Test: Richards, boundary solver" << std::endl;

  // read parameter list
  std::string xmlFileName = "test/flow_richards_boundary_solver.xml";
  Teuchos::RCP<Teuchos::ParameterList> plist = Teuchos::getParametersFromXmlFile(xmlFileName);

  // create a mesh framework
  Teuchos::ParameterList regions_list = plist->get<Teuchos::ParameterList>("regions");
  Teuchos::RCP<Amanzi::AmanziGeometry::GeometricModel> gm =
    Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(3, regions_list, *comm));


  Preference pref;
  pref.clear();
  pref.push_back(Framework::MSTK);

  double bottom = -0.5;
  MeshFactory meshfactory(comm, gm);
  meshfactory.set_preference(pref);
  // Teuchos::RCP<const Mesh> mesh = meshfactory.create(0.0, bottom, 1.0, 0.0, 1, 10);
  Teuchos::RCP<const Mesh> mesh1 = meshfactory.create("test/hex_2x2x1-1.exo");
  Teuchos::RCP<const Mesh> mesh2 = meshfactory.create("test/hex_2x2x1-2.exo");

  // // create a simple state and populate it
  Teuchos::ParameterList state_list = plist->sublist("state");
  Teuchos::RCP<State> S1 = Teuchos::rcp(new State(state_list));
  Teuchos::RCP<State> S2 = Teuchos::rcp(new State(state_list));
  S1->RegisterDomainMesh(Teuchos::rcp_const_cast<Mesh>(mesh1));
  S2->RegisterDomainMesh(Teuchos::rcp_const_cast<Mesh>(mesh2));

  Teuchos::RCP<TreeVector> soln1 = Teuchos::rcp(new TreeVector());
  Teuchos::RCP<Richards_PK> RPK1 = Teuchos::rcp(new Richards_PK(plist, "flow", S1, soln1));
  Teuchos::RCP<TreeVector> soln2 = Teuchos::rcp(new TreeVector());
  Teuchos::RCP<Richards_PK> RPK2 = Teuchos::rcp(new Richards_PK(plist, "flow", S2, soln2));

  RPK1->Setup();
  S1->Setup();
  S1->InitializeFields();
  S1->InitializeEvaluators();

  RPK2->Setup();
  S2->Setup();
  S2->InitializeFields();
  S2->InitializeEvaluators();

  // modify the default state for the problem at hand
  std::string passwd("");
  auto& K1 = *S1->GetW<CompositeVector>("permeability", "permeability").ViewComponent("cell");

  AmanziMesh::Entity_ID_List block;
  mesh1->getSetEntities(
    "All", AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_type::OWNED, &block);
  for (int i = 0; i != block.size(); ++i) {
    int c = block[i];
    K1[0][c] = 1e-9;
    K1[1][c] = 1e-9;
    K1[2][c] = 1e-9;
  }
  S1->GetRecordW("permeability", "permeability").set_initialized();

  auto& K2 = *S2->GetW<CompositeVector>("permeability", "permeability").ViewComponent("cell");
  mesh2->getSetEntities(
    "Material 1", AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_type::OWNED, &block);
  for (int i = 0; i != block.size(); ++i) {
    int c = block[i];
    K2[0][c] = 1e-9;
    K2[1][c] = 1e-9;
    K2[2][c] = 1e-9;
  }
  S2->GetRecordW("permeability", "permeability").set_initialized();

  double atm_pressure = 101325.0;
  double rho = S1->Get<double>("const_fluid_density");

  auto gravity1 = S1->Get<AmanziGeometry::Point>("gravity");
  auto gravity2 = S2->Get<AmanziGeometry::Point>("gravity");
  S1->GetRecordW("gravity", "state").set_initialized();
  S2->GetRecordW("gravity", "state").set_initialized();

  // create the initial pressure function
  auto& p1 = *S1->GetW<CompositeVector>("pressure", passwd).ViewComponent("cell");
  auto& p2 = *S2->GetW<CompositeVector>("pressure", passwd).ViewComponent("cell");

  for (int c = 0; c < p1.MyLength(); c++) {
    const Point& xc = mesh1->getCellCentroid(c);
    p1[0][c] = 0.6 * atm_pressure + rho * gravity1[1] * (xc[1] - bottom);
    p2[0][c] = 0.6 * atm_pressure + rho * gravity2[1] * (xc[1] - bottom);
  }

  // initialize the Richards process kernel
  RPK1->Initialize();
  S1->CheckAllFieldsInitialized();
  RPK2->Initialize();
  S2->CheckAllFieldsInitialized();

  std::cout << p1 << "\n";
  std::cout << p2 << "\n";

  double bnd_val1, bnd_val2;

  std::cout << "MESH1\n";

  int nfaces =
    mesh1->getNumEntities(AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_type::OWNED);
  for (int f = 0; f < nfaces; f++) {
    AmanziMesh::Entity_ID_List cells;
    cells = mesh1->getFaceCells(f, AmanziMesh::Parallel_type::ALL);
    int dir;
    const Point& norm = mesh1->getFaceNormal(f, cells[0], &dir);
    if ((cells.size() == 1) && (norm[2] * dir > 0)) {
      bnd_val1 = RPK1->BoundaryFaceValue(f, S1->GetW<CompositeVector>("pressure", passwd));
      std::cout << ": " << f << " " << bnd_val1 << "\n";
    }
  }

  std::cout << "MESH2\n";

  nfaces = mesh2->getNumEntities(AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_type::OWNED);
  for (int f = 0; f < nfaces; f++) {
    AmanziMesh::Entity_ID_List cells;
    cells = mesh2->getFaceCells(f, AmanziMesh::Parallel_type::ALL);
    int dir;
    const Point& norm = mesh2->getFaceNormal(f, cells[0], &dir);
    if ((cells.size() == 1) && (norm[2] * dir > 0)) {
      bnd_val2 = RPK2->BoundaryFaceValue(f, S2->GetW<CompositeVector>("pressure", passwd));
      std::cout << ": " << f << " " << bnd_val2 << "\n";
    }
  }
}
