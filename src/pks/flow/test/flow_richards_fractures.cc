/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Transport PK

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
#include "IO.hh"
#include "MeshFactory.hh"
#include "Mesh_MSTK.hh"
#include "OutputXDMF.hh"
#include "State.hh"

// Flow
#include "Richards_PK.hh"

/* *******************************************************************
* Compute pressure.
******************************************************************* */
double
MeanPressure(const Amanzi::State& S)
{
  const auto& p = *S.Get<Amanzi::CompositeVector>("pressure").ViewComponent("cell");
  int ncells = p.MyLength();

  double mean(0.0);
  for (int c = 0; c < ncells; c++) { mean += p[0][c]; }
  return mean / ncells;
}


/* *******************************************************************
* Run the model.
******************************************************************* */
TEST(RICHARDS_TWO_FRACTURES)
{
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
  meshfactory.set_preference(Preference({ Framework::MSTK }));
  RCP<const Mesh> mesh3D = meshfactory.create(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 10, 10, 10);

  // extract fractures mesh
  std::vector<std::string> setnames;
  setnames.push_back("fracture 1");
  setnames.push_back("fracture 2");

  RCP<const Mesh> mesh = meshfactory.create(mesh3D, setnames, AmanziMesh::Entity_kind::FACE);

  int ncells_owned =
    mesh->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);
  int ncells_wghost =
    mesh->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::ALL);

  std::cout << "pid=" << MyPID << " cells: " << ncells_owned << " " << ncells_wghost << std::endl;

  Teuchos::ParameterList state_list = plist->sublist("state");
  RCP<State> S = rcp(new State(state_list));
  S->RegisterDomainMesh(rcp_const_cast<Mesh>(mesh));

  Teuchos::ParameterList pk_tree("flow");
  Teuchos::RCP<TreeVector> soln = Teuchos::rcp(new TreeVector());
  Teuchos::RCP<Richards_PK> RPK = Teuchos::rcp(new Richards_PK(pk_tree, plist, S, soln));
  RPK->parseParameterList();
  RPK->Setup();

  S->Setup();
  S->InitializeFields();
  S->InitializeEvaluators();

  // initialize the Darcy process kernel
  RPK->Initialize();
  S->CheckAllFieldsInitialized();
  WriteStateStatistics(*S);

  // initialize io
  Teuchos::ParameterList iolist;
  iolist.get<std::string>("file name base", "plot");
  auto io = Teuchos::rcp(new OutputXDMF(iolist, mesh, true, false));

  // transient solution
  double t_old(0.0), t_new, dt(2.5);
  for (int n = 0; n < 10; n++) {
    t_new = t_old + dt;

    RPK->AdvanceStep(t_old, t_new);
    RPK->CommitStep(t_old, t_new, Tags::DEFAULT);

    t_old = t_new;
    dt *= 1.2;

    // measure of wetting front speed
    double pmean = MeanPressure(*S);
    CHECK(pmean > 80000 + (t_old - 2) * 160);

    io->InitializeCycle(t_old, n, "");
    const auto& u0 = *S->Get<CompositeVector>("pressure").ViewComponent("cell");
    const auto& u1 = *S->Get<CompositeVector>("saturation_liquid").ViewComponent("cell");
    const auto& u2 = *S->Get<CompositeVector>("hydraulic_head").ViewComponent("cell");

    io->WriteVector(*u0(0), "pressure", AmanziMesh::Entity_kind::CELL);
    io->WriteVector(*u1(0), "saturation", AmanziMesh::Entity_kind::CELL);
    io->WriteVector(*u2(0), "hydraulic head", AmanziMesh::Entity_kind::CELL);
    io->FinalizeCycle();
  }
  WriteStateStatistics(*S);

  const auto& p = *S->Get<CompositeVector>("pressure").ViewComponent("cell");
  for (int c = 0; c < p.MyLength(); c++) { CHECK(p[0][c] > 76000.0 && p[0][c] < 101000.0); }
}
