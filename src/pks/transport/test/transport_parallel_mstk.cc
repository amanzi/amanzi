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
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "UnitTest++.h"

// Amanzi
#include "MeshFactory.hh"
#include "State.hh"

// Transport
#include "TransportExplicit_PK.hh"


double
f_step(const Amanzi::AmanziGeometry::Point& x, double t)
{
  if (x[0] <= t) return 1;
  return 0;
}


TEST(ADVANCE_WITH_MSTK_PARALLEL)
{
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::Transport;
  using namespace Amanzi::AmanziGeometry;

  std::cout << "Test: second-order transport with parallel MSTK framework" << std::endl;
  Comm_ptr_type comm = Amanzi::getDefaultComm();

  // read parameter list
  std::string xmlFileName = "test/transport_parallel_mstk.xml";
  Teuchos::RCP<Teuchos::ParameterList> plist = Teuchos::getParametersFromXmlFile(xmlFileName);

  /* create an MSTK mesh framework */
  ParameterList region_list = plist->get<Teuchos::ParameterList>("regions");
  Teuchos::RCP<Amanzi::AmanziGeometry::GeometricModel> gm =
    Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(3, region_list, *comm));

  MeshFactory meshfactory(comm, gm);
  meshfactory.set_preference(Preference({ Framework::MSTK }));
  RCP<const Mesh> mesh = meshfactory.create("test/hex_3x3x3_ss.exo");

  /* create a simple state and populate it */
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

  /* modify the default state for the problem at hand */
  std::string passwd("state");
  auto& flux = *S->GetW<CompositeVector>("volumetric_flow_rate", passwd).ViewComponent("face");

  AmanziGeometry::Point velocity(1.0, 0.0, 0.0);
  int nfaces_owned =
    mesh->getNumEntities(AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_kind::OWNED);
  for (int f = 0; f < nfaces_owned; f++) {
    const AmanziGeometry::Point& normal = mesh->getFaceNormal(f);
    flux[0][f] = velocity * normal;
  }

  auto& tcc =
    *S->GetW<CompositeVector>("total_component_concentration", passwd).ViewComponent("cell");

  int ncells_owned =
    mesh->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);
  for (int c = 0; c < ncells_owned; c++) {
    const AmanziGeometry::Point& xc = mesh->getCellCentroid(c);
    tcc[0][c] = f_step(xc, 0.0);
  }

  // initialize a transport process kernel
  TPK.Initialize();

  /* advance the state */
  double t_old(0.0), t_new, dt;
  dt = TPK.StableTimeStep(-1);
  t_new = t_old + dt;
  TPK.AdvanceStep(t_old, t_new);

  // print cell concentrations
  int iter = 0;
  while (t_new < 1.0) {
    dt = TPK.StableTimeStep(-1);
    t_new = t_old + dt;

    TPK.AdvanceStep(t_old, t_new);
    TPK.CommitStep(t_old, t_new, Tags::DEFAULT);

    t_old = t_new;
    iter++;

    if (iter < 10 && TPK.MyPID == 3) {
      printf("T=%7.2f  C_0(x):", t_new);
      for (int k = 0; k < 2; k++) printf("%7.4f", tcc[0][k]);
      std::cout << std::endl;
    }

    for (int k = 0; k < ncells_owned; ++k) CHECK(tcc[0][k] >= 0.0 && tcc[0][k] <= 1.0);
  }

  for (int k = 0; k < ncells_owned; ++k) CHECK_CLOSE(tcc[0][k], 1.0, 1e-6);
}
