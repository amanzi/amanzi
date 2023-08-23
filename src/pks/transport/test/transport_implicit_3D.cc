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
#include "IO.hh"
#include "MeshFactory.hh"
#include "State.hh"

// Transport
#include "TransportImplicit_PK.hh"


void
runTest(int order, const std::string& linsolver)
{
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::Transport;
  using namespace Amanzi::AmanziGeometry;

  std::string framework_name = "MSTK";
  Framework framework = Framework::MSTK;

  std::cout << "\nTEST: implicit advance, spatial order=" << order << std::endl;

  Comm_ptr_type comm = Amanzi::getDefaultComm();

  // read parameter list
  std::string xmlFileName("test/transport_implicit_3D.xml");
  Teuchos::RCP<Teuchos::ParameterList> plist = Teuchos::getParametersFromXmlFile(xmlFileName);

  plist->sublist("PKs")
    .sublist("transport implicit")
    .set<int>("spatial discretization order", order);
  plist->sublist("PKs")
    .sublist("transport implicit")
    .sublist("time integrator")
    .set<std::string>("linear solver", linsolver);

  // create a mesh
  ParameterList region_list = plist->get<Teuchos::ParameterList>("regions");
  auto gm = Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(3, region_list, *comm));

  Preference pref;
  pref.clear();
  pref.push_back(framework);

  MeshFactory meshfactory(comm, gm);
  meshfactory.set_preference(pref);

  int nx(7);
  // RCP<const Mesh> mesh = meshfactory.create("test/hex_3x3x3_ss.exo");
  RCP<const Mesh> mesh = meshfactory.create(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, nx, nx, nx);

  // create a simple state and populate it
  Amanzi::VerboseObject::global_hide_line_prefix = false;

  std::vector<std::string> component_names;
  component_names.push_back("Component 0");

  Teuchos::ParameterList state_list = plist->sublist("state");
  RCP<State> S = rcp(new State(state_list));
  S->RegisterDomainMesh(rcp_const_cast<Mesh>(mesh));

  Teuchos::ParameterList pk_tree =
    plist->sublist("cycle driver").sublist("pk_tree").sublist("transport implicit");

  // create the global solution vector
  Teuchos::RCP<TreeVector> soln = Teuchos::rcp(new TreeVector());
  TransportImplicit_PK TPK(pk_tree, plist, S, soln);

  TPK.Setup();
  S->Setup();
  S->InitializeFields();
  S->InitializeEvaluators();
  S->set_time(0.0);
  S->set_intermediate_time(0.0);

  // modify the default state for the problem at hand
  std::string passwd("state");
  auto& flux = *S->GetW<CompositeVector>("volumetric_flow_rate", passwd).ViewComponent("face");

  AmanziGeometry::Point velocity(1.0, 1.0, 0.0);
  int nfaces_owned =
    mesh->getNumEntities(AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_type::OWNED);
  for (int f = 0; f < nfaces_owned; f++) {
    const AmanziGeometry::Point& normal = mesh->getFaceNormal(f);
    flux[0][f] = velocity * normal;
  }

  // initialize a transport process kernel
  TPK.Initialize();

  // advance the state
  int nloop(0);
  bool failed;
  double t_old(0.0), t_new(0.0), dt(0.01);
  auto tcc =
    S->GetW<CompositeVector>("total_component_concentration", passwd).ViewComponent("cell");
  tcc->PutScalar(0.0);

  while (t_new < 1.2) {
    t_new = t_old + dt;

    if (comm->MyPID() == 0)
      std::cout << "\nCycle: " << nloop << " T=" << t_old << " dT=" << dt << std::endl;

    failed = TPK.AdvanceStep(t_old, t_new);
    if (failed) {
      dt /= 2;
    } else {
      TPK.CommitStep(t_old, t_new, Tags::DEFAULT);
      dt = std::min(dt * 1.1, 0.01);
      t_old = t_new;
      nloop++;

      if (t_new < 0.4) {
        printf("T=%6.2f  C_0(x):", t_new);
        for (int k = 0; k < nx * nx * nx; k += nx * nx) printf("%7.4f", (*tcc)[0][k]);
        std::cout << std::endl;
      }
    }
  }

  // check that the final state is constant
  for (int k = 0; k < 4; k++) CHECK_CLOSE((*tcc)[0][k], 1.0, 1e-5);

  WriteStateStatistics(*S);
}


TEST(IMPLICIT_TRANSPORT_3D_FIRST_ORDER)
{
  runTest(1, "PCG");
}

TEST(IMPLICIT_TRANSPORT_3D_SECOND_ORDER)
{
  runTest(2, "");
}
