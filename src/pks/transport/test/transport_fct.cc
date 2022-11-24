/*
  Transport PK
  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
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


TEST(ADVANCE_FCT)
{
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::Transport;
  using namespace Amanzi::AmanziGeometry;

  std::cout << "Test: flux-corrected transport" << std::endl;
  Comm_ptr_type comm = Amanzi::getDefaultComm();

  // read parameter list
  std::string xmlFileName("test/transport_fct.xml");
  Teuchos::RCP<Teuchos::ParameterList> plist = Teuchos::getParametersFromXmlFile(xmlFileName);

  // create a mesh
  ParameterList region_list = plist->get<Teuchos::ParameterList>("regions");
  auto gm = Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(2, region_list, *comm));

  MeshFactory meshfactory(comm, gm);
  meshfactory.set_preference(Preference({ Framework::MSTK }));
  RCP<const Mesh> mesh = meshfactory.create(0.0, 0.0, 1.0, 1.0, 20, 1);

  int ncells_owned = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  int nfaces_wghost = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::ALL);

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
  auto& flux =
    *S->GetW<CompositeVector>("volumetric_flow_rate", passwd).ViewComponent("face", true);

  AmanziGeometry::Point velocity(1.0, 0.0);
  for (int f = 0; f < nfaces_wghost; ++f) {
    const AmanziGeometry::Point& normal = mesh->face_normal(f);
    flux[0][f] = velocity * normal;
  }

  // initialize a transport process kernel
  TPK.Initialize();

  // advance the state
  int p0(1);
  double t_old(0.0), t_new, dt, tol(1e-10);
  dt = TPK.StableTimeStep(-1);
  t_new = t_old + dt;
  TPK.AdvanceStep(t_old, t_new);

  // printing cell concentration
  auto tcc =
    S->GetW<CompositeVector>("total_component_concentration", passwd).ViewComponent("cell");

  while (t_new < 0.2) {
    dt = TPK.StableTimeStep(-1);
    t_new = t_old + dt;

    TPK.AdvanceStep(t_old, t_new);
    TPK.CommitStep(t_old, t_new, Tags::DEFAULT);

    t_old = t_new;
    S->set_intermediate_time(t_old);

    if (comm->MyPID() == p0) {
      printf("T=%7.3f  C_0(x):", t_new);
      for (int k = 0; k < 9; k++) printf("%11.4g", (*tcc)[0][k]);
      std::cout << std::endl;
      for (int k = 0; k < ncells_owned; k++) CHECK((*tcc)[0][k] < 1.0 + tol && (*tcc)[0][k] > -tol);
    }
  }

  if (comm->MyPID() == p0) {
    for (int k = 0; k < 4; k++) CHECK_CLOSE((*tcc)[0][k], 1.0, 1e-6);
  }

  // turn off boundary conditions
  while (t_new < 0.4) {
    dt = TPK.StableTimeStep(-1);
    t_new = t_old + dt;

    TPK.AdvanceStep(t_old, t_new);
    TPK.CommitStep(t_old, t_new, Tags::DEFAULT);

    t_old = t_new;
    S->set_intermediate_time(t_old);

    if (comm->MyPID() == p0) {
      printf("T=%7.3f  C_0(x):", t_new);
      for (int k = 0; k < 9; k++) printf("%11.4g", (*tcc)[0][k]);
      std::cout << std::endl;
      for (int k = 0; k < ncells_owned; k++) CHECK((*tcc)[0][k] < 1.0 + tol && (*tcc)[0][k] > -tol);
    }
  }

  if (comm->MyPID() == p0) {
    for (int k = 0; k < 4; k++) CHECK_CLOSE((*tcc)[0][k], 0.0, 1e-6);
  }
}
