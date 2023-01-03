/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <vector>

// TPLs
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "UnitTest++.h"

// Amanzi
#include "MeshFactory.hh"
#include "State.hh"

// Transport
#include "TransportExplicit_PK.hh"


TEST(ADVANCE_WITH_MULTISCALE)
{
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::Transport;
  using namespace Amanzi::AmanziGeometry;

  std::cout << "Test: advance multiscale " << std::endl;
  Comm_ptr_type comm = Amanzi::getDefaultComm();

  // read parameter list
  std::string xmlFileName("test/transport_multiscale.xml");
  Teuchos::RCP<Teuchos::ParameterList> plist = Teuchos::getParametersFromXmlFile(xmlFileName);

  // create a mesh
  ParameterList region_list = plist->get<Teuchos::ParameterList>("regions");
  auto gm = Teuchos::rcp(new AmanziGeometry::GeometricModel(3, region_list, *comm));

  Preference pref;
  pref.clear();
  pref.push_back(Framework::MSTK);

  MeshFactory meshfactory(comm, gm);
  meshfactory.set_preference(pref);
  RCP<const Mesh> mesh = meshfactory.create(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 20, 1, 1);

  // create a simple state and populate it
  Amanzi::VerboseObject::global_hide_line_prefix = false;

  std::vector<std::string> component_names;
  component_names.push_back("Component 0");
  component_names.push_back("Component 1");

  Teuchos::ParameterList state_list = plist->sublist("state");
  RCP<State> S = rcp(new State(state_list));
  S->RegisterDomainMesh(rcp_const_cast<Mesh>(mesh));

  TransportExplicit_PK TPK(plist, S, "transport", component_names);
  TPK.Setup();
  S->Setup();
  S->InitializeFields();
  S->InitializeEvaluators();
  S->set_time(0.0);
  S->set_intermediate_time(0.0);

  // modify the default state for the problem at hand
  std::string passwd("state");
  auto& flux = *S->GetW<CompositeVector>("volumetric_flow_rate", passwd).ViewComponent("face");

  AmanziGeometry::Point velocity(1.0, 0.0, 0.0);
  int nfaces_owned = mesh->getNumEntities(AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_type::OWNED);
  for (int f = 0; f < nfaces_owned; f++) {
    const AmanziGeometry::Point& normal = mesh->getFaceNormal(f);
    flux[0][f] = velocity * normal;
  }

  // initialize a transport process kernel
  TPK.Initialize();

  // advance the state
  double t_old(0.0), t_new, dt;
  dt = TPK.StableTimeStep(-1);
  t_new = t_old + dt;
  TPK.AdvanceStep(t_old, t_new);

  // printing cell concentration
  auto& tcc = *S->Get<CompositeVector>("total_component_concentration").ViewComponent("cell");
  auto& tcc_msp =
    *S->Get<CompositeVector>("total_component_concentration_msp").ViewComponent("cell");

  while (t_new < 0.2) {
    dt = TPK.StableTimeStep(-1);
    t_new = t_old + dt;

    TPK.AdvanceStep(t_old, t_new);
    TPK.CommitStep(t_old, t_new, Tags::DEFAULT);

    t_old = t_new;

    printf("T=%6.3f  C_0(x):", t_new);
    for (int k = 0; k < 9; k++) printf("%9.2g", tcc[0][k]);
    std::cout << std::endl;
    printf("T=%6.3f  C_1(x):", t_new);
    for (int k = 0; k < 9; k++) printf("%9.2g", tcc[1][k]);
    std::cout << std::endl;
    printf("T=%6.3f  c_0(x):", t_new);
    for (int k = 0; k < 9; k++) printf("%9.2g", tcc_msp[0][k]);
    std::cout << std::endl;
  }

  // check that the final state is constant
  for (int k = 0; k < 6; k++) CHECK_CLOSE(tcc[0][k], 1.0, 1e-6);

  for (int k = 0; k < 20; k++) CHECK_CLOSE(tcc[0][k], tcc[1][k], 1e-12);

  for (int k = 0; k < 19; k++) {
    CHECK(tcc[0][k] - tcc[0][k + 1] > -1e-15);
    CHECK(tcc_msp[0][k] - tcc_msp[0][k + 1] > -1e-15);
  }
}
