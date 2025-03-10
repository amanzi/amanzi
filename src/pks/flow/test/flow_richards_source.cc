/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Flow PK

  Various sources: verify solution bounds and PK interfaces
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
#include "IO.hh"
#include "MeshFactory.hh"
#include "OutputXDMF.hh"
#include "State.hh"
// #include "state_evaluators_registration.hh"

// Flow
#include "Richards_PK.hh"
#include "Richards_SteadyState.hh"

void
Richards2D_SourceTest(std::string filename, bool deform)
{
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Flow;

  Comm_ptr_type comm = Amanzi::getDefaultComm();
  int MyPID = comm->MyPID();
  if (MyPID == 0) std::cout << "Test: 2D Richards, various sources" << std::endl;

  // create a mesh framework
  Teuchos::RCP<Teuchos::ParameterList> plist = Teuchos::getParametersFromXmlFile(filename);
  ParameterList regions_list = plist->get<Teuchos::ParameterList>("regions");
  auto gm = Teuchos::rcp(new AmanziGeometry::GeometricModel(2, regions_list, *comm));

  MeshFactory meshfactory(comm, gm);
  meshfactory.set_preference(Preference({ Framework::MSTK }));
  RCP<Mesh> mesh = meshfactory.create(0.0, 0.0, 100.0, 50.0, 50, 25);

  if (deform) {
    int nnodes =
      mesh->getNumEntities(AmanziMesh::Entity_kind::NODE, AmanziMesh::Parallel_kind::ALL);
    AmanziMesh::Entity_ID_View nodeids("nodeids", nnodes);
    AmanziMesh::Point_View new_positions("new_positions", nnodes);

    for (int v = 0; v < nnodes; ++v) {
      auto xv = mesh->getNodeCoordinate(v);
      nodeids[v] = v;
      xv[1] *= (xv[0] * 0.4 + (100.0 - xv[0])) / 100.0;
      new_positions[v] = xv;
    }
    AmanziMesh::deform(*mesh, nodeids, new_positions);
  }

  // create State and PK
  Amanzi::VerboseObject::global_hide_line_prefix = true;

  Teuchos::ParameterList state_list = plist->get<Teuchos::ParameterList>("state");
  RCP<State> S = rcp(new State(state_list));
  S->RegisterDomainMesh(rcp_const_cast<Mesh>(mesh));

  Teuchos::ParameterList pk_tree("flow");
  Teuchos::RCP<TreeVector> soln = Teuchos::rcp(new TreeVector());
  Richards_PK* RPK = new Richards_PK(pk_tree, plist, S, soln);

  RPK->parseParameterList();
  RPK->Setup();
  S->Setup();
  S->InitializeFields();
  S->InitializeEvaluators();

  RPK->Initialize();
  S->CheckAllFieldsInitialized();

  // solve the steady-state problem
  TI_Specs ti_specs;
  ti_specs.T0 = 0.0;
  ti_specs.dT0 = 10.0;
  ti_specs.T1 = 4000.0;
  ti_specs.max_itrs = 3000;

  AdvanceToSteadyState(S, *RPK, ti_specs, soln);
  RPK->CommitStep(0.0, 1.0, Tags::DEFAULT); // dummy times for flow

  WriteStateStatistics(*S);

  // output
  const auto& e = *S->Get<CompositeVector>("strain_rate").ViewComponent("cell");
  const auto& s = *S->Get<CompositeVector>("saturation_liquid").ViewComponent("cell");
  const auto& p = *S->Get<CompositeVector>("pressure").ViewComponent("cell");

  Teuchos::ParameterList iolist;
  iolist.get<std::string>("file name base", "plot");
  OutputXDMF io(iolist, mesh, true, false);
  io.InitializeCycle(ti_specs.T1, 1, "");
  io.WriteVector(*p(0), "pressure", AmanziMesh::Entity_kind::CELL);
  io.WriteVector(*s(0), "saturation", AmanziMesh::Entity_kind::CELL);
  io.WriteVector(*e(0), "strain rate", AmanziMesh::Entity_kind::CELL);
  io.FinalizeCycle();

  delete RPK;
}

TEST(FLOW_2D_RICHARDS_SOURCE)
{
  Richards2D_SourceTest("test/flow_richards_source.xml", true);
}
