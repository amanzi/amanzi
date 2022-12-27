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

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <vector>

// TPLs
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "UnitTest++.h"

// Amanzi
#include "GMVMesh.hh"
#include "IO.hh"
#include "MeshFactory.hh"
#include "MeshAudit.hh"
#include "State.hh"

// Transport
#include "TransportExplicit_PK.hh"

/* **************************************************************** */
void
runTest(std::string xmlfile)
{
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::Transport;
  using namespace Amanzi::AmanziGeometry;

  std::cout << "Test: Advance over a dwet/dry soil" << std::endl;

  Comm_ptr_type comm = Amanzi::getDefaultComm();
  Amanzi::VerboseObject::global_hide_line_prefix = true;

  // read parameter list
  Teuchos::RCP<Teuchos::ParameterList> plist = Teuchos::getParametersFromXmlFile(xmlfile);

  // create a mesh framework
  ParameterList region_list = plist->get<Teuchos::ParameterList>("regions");
  auto gm = Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(2, region_list, *comm));

  MeshFactory meshfactory(comm, gm);
  meshfactory.set_preference(Preference({ Framework::MSTK }));
  RCP<Mesh> mesh = meshfactory.create(0.0, 0.0, 1.0, 1.0, 10, 10);

  // create a state and PK
  std::vector<std::string> component_names;
  component_names.push_back("Component 0");

  Teuchos::ParameterList state_list = plist->sublist("state");
  RCP<State> S = rcp(new State(state_list));
  S->RegisterDomainMesh(mesh);

  TransportExplicit_PK TPK(plist, S, "transport", component_names);
  TPK.Setup();
  S->Setup();
  S->InitializeFields();
  S->InitializeEvaluators();
  S->set_time(0.0);
  TPK.Initialize();
  WriteStateStatistics(*S);

  // advance the transport state
  int iter(0);
  double t_old(0.0), t_new(0.0), dt;
  auto tcc =
    S->GetW<CompositeVector>("total_component_concentration", "state").ViewComponent("cell");

  while (t_new < 0.5) {
    dt = TPK.StableTimeStep(-1);
    t_new = t_old + dt;

    S->set_initial_time(t_old);
    S->set_intermediate_time(t_old);
    S->set_final_time(t_new);

    TPK.AdvanceStep(t_old, t_new);
    TPK.CommitStep(t_old, t_new, Tags::DEFAULT);

    t_old = t_new;
    iter++;

    printf("T=%8.4f  C_0(x):", t_new);
    for (int k = 3; k < 100; k += 10) { printf("%8.5f", (*tcc)[0][k]); }
    printf("\n");

    for (int c = 0; c < 10; ++c) { CHECK((*tcc)[0][c] <= 1.0 && (*tcc)[0][c] >= 0.0); }
  }

  WriteStateStatistics(*S);
}


TEST(ADVANCE_2D_MESH)
{
  runTest("test/transport_wet_dry.xml");
}
