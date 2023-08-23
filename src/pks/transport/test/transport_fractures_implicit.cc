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
#include "UniqueLocalIndex.hh"

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

  std::cout << "\nTest: Implicit transport in fractures, order=" << order << std::endl;
#ifdef HAVE_MPI
  Comm_ptr_type comm = Amanzi::getDefaultComm();
#else
  Epetra_SerialComm* comm = new Epetra_SerialComm();
#endif

  // read parameter list
  std::string xmlFileName = "test/transport_fractures_implicit.xml";
  Teuchos::RCP<Teuchos::ParameterList> plist = Teuchos::getParametersFromXmlFile(xmlFileName);

  // read and modify parameter list
  ParameterList region_list = plist->get<Teuchos::ParameterList>("regions");
  auto gm = Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(3, region_list, *comm));

  plist->sublist("PKs").sublist("transport").set<int>("spatial discretization order", order);
  plist->sublist("PKs")
    .sublist("transport")
    .sublist("time integrator")
    .set<std::string>("linear solver", linsolver);

  // create meshes
  MeshFactory meshfactory(comm, gm);
  meshfactory.set_preference(Preference({ Framework::MSTK, Framework::STK }));
  RCP<const Mesh> mesh3D = meshfactory.create(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 10, 10, 10);

  std::vector<std::string> setnames;
  setnames.push_back("fracture 1");
  setnames.push_back("fracture 2");

  // RCP<const Mesh> mesh = meshfactory.create(mesh3D, setnames, AmanziMesh::Entity_kind::FACE);
  RCP<const Mesh> mesh = meshfactory.create("test/fractures.exo");

  int ncells_owned =
    mesh->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_type::OWNED);
  int ncells_wghost =
    mesh->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_type::ALL);
  int nfaces_owned =
    mesh->getNumEntities(AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_type::OWNED);
  std::cout << "pid=" << comm->MyPID() << " cells: " << ncells_owned << " faces: " << nfaces_owned
            << std::endl;

  // initialize io
  Teuchos::ParameterList iolist;
  iolist.get<std::string>("file name base", "transport");
  OutputXDMF io(iolist, mesh, true, false);

  // create a simple state and populate it
  Amanzi::VerboseObject::global_hide_line_prefix = true;

  std::vector<std::string> component_names;
  component_names.push_back("tracer");

  Teuchos::ParameterList state_list = plist->sublist("state");
  RCP<State> S = rcp(new State(state_list));
  S->RegisterDomainMesh(rcp_const_cast<Mesh>(mesh));

  Teuchos::ParameterList pk_tree =
    plist->sublist("cycle driver").sublist("pk_tree").sublist("transport");
  Teuchos::RCP<TreeVector> soln = Teuchos::rcp(new TreeVector());
  TransportImplicit_PK TPK(pk_tree, plist, S, soln);
  TPK.Setup();
  S->Setup();
  S->InitializeFields();
  S->InitializeEvaluators();
  S->set_time(0.0);
  S->set_intermediate_time(0.0);
  S->set_initial_time(0.0);
  S->set_final_time(0.0);

  // modify the default state
  auto& flux =
    *S->GetW<CompositeVector>("volumetric_flow_rate", "state").ViewComponent("face", true);
  const auto flux_map =
    S->GetW<CompositeVector>("volumetric_flow_rate", "state").Map().Map("face", true);

  int dir;
  AmanziGeometry::Point velocity(1.0, 0.2, -0.2);
  velocity /= 4e+5;

  for (int c = 0; c < ncells_wghost; c++) {
    const auto& faces = mesh->getCellFaces(c);
    int nfaces = faces.size();

    for (int n = 0; n < nfaces; ++n) {
      int f = faces[n];
      int g = flux_map->FirstPointInElement(f);
      int ndofs = flux_map->ElementSize(f);
      if (ndofs > 1) g += Operators::UniqueIndexFaceToCells(*mesh, f, c);

      const AmanziGeometry::Point& normal = mesh->getFaceNormal(f, c, &dir);
      flux[0][g] = (velocity * normal) * dir;
    }
  }
  S->GetRecordW("volumetric_flow_rate", "state").set_initialized();

  // initialize the transport process kernel
  TPK.Initialize();

  // advance the transport state
  int iter(0);
  double tol = (order == 1) ? 1e-8 : 1e-4;
  double t_old(0.0), t_new(0.0), dt(2.0e+3);
  auto& tcc =
    *S->GetW<CompositeVector>("total_component_concentration", "state").ViewComponent("cell");

  while (t_new < 1.0e+5) {
    t_new = t_old + dt;

    if (comm->MyPID() == 0) std::cout << "\nCycle: T=" << t_old << " dT=" << dt << std::endl;

    TPK.AdvanceStep(t_old, t_new);
    TPK.CommitStep(t_old, t_new, Tags::DEFAULT);

    t_old = t_new;
    iter++;

    // verify solution
    for (int c = 0; c < ncells_owned; ++c) CHECK(tcc[0][c] >= -tol && tcc[0][c] <= 1.0 + tol);

    if (iter % 2 == 0) {
      io.InitializeCycle(t_new, iter, "");
      io.WriteVector(*tcc(0), "tcc", AmanziMesh::Entity_kind::CELL);
      io.FinalizeCycle();
    }
  }

  // test the maximum principle
  AmanziMesh::Entity_ID_List block;
  mesh->getSetEntities(
    "fracture 2", AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_type::OWNED, &block);

  // test that solute enter the second fracture
  double tcc_max(0.0);
  for (int n = 0; n < block.size(); ++n) {
    tcc_max = std::max(tcc_max, tcc[0][block[n]]);
    // std::cout << n << " " << tcc[0][block[n]] << " " << tcc_max << std::endl;
  }
  double tmp = tcc_max;
  mesh->getComm()->MaxAll(&tmp, &tcc_max, 1);
  CHECK(tcc_max > 0.25);

  WriteStateStatistics(*S);
}


TEST(IMPLICIT_TRANSPORT_TWO_FRACTURES_1ST)
{
  runTest(1, "PCG");
}


TEST(IMPLICIT_TRANSPORT_TWO_FRACTURES_2ND)
{
  runTest(2, "GMRES");
}
