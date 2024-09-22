/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Meachanics PK

*/

// TPLs
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "UnitTest++.h"

// Amanzi
#include "IO.hh"
#include "Mesh.hh"
#include "MeshFactory.hh"
#include "OutputXDMF.hh"

#include "MechanicsElasticity_PK.hh"

/* **************************************************************** */
TEST(COOK_MEMBRANE_2D)
{
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Mechanics;

  Comm_ptr_type comm = Amanzi::getDefaultComm();
  int MyPID = comm->MyPID();
  if (MyPID == 0) std::cout << "Test: Elastic deformation in 2D" << std::endl;

  // read parameter list
  std::string xmlFileName = "test/mechanics_cook.xml";
  Teuchos::RCP<Teuchos::ParameterList> plist = Teuchos::getParametersFromXmlFile(xmlFileName);

  // create a mesh framework
  Teuchos::ParameterList regions_list = plist->sublist("regions");
  auto gm = Teuchos::rcp(new GeometricModel(2, regions_list, *comm));

  MeshFactory meshfactory(comm, gm);
  meshfactory.set_preference(Preference({ Framework::MSTK }));
  Teuchos::RCP<const Mesh> mesh = meshfactory.create("test/mechanics_cook.exo");

  // create a simple state and populate it
  Teuchos::ParameterList state_list = plist->sublist("state");
  Teuchos::RCP<State> S = Teuchos::rcp(new State(state_list));
  S->RegisterDomainMesh(Teuchos::rcp_const_cast<Mesh>(mesh));

  Teuchos::RCP<TreeVector> soln = Teuchos::rcp(new TreeVector());
  Teuchos::ParameterList pk_tree = plist->sublist("PK tree").sublist("mechanics elasticity");
  auto MPK = Teuchos::rcp(new MechanicsElasticity_PK(pk_tree, plist, S, soln));

  MPK->Setup();
  S->Setup();
  S->InitializeFields();
  S->InitializeEvaluators();

  MPK->Initialize();
  S->CheckAllFieldsInitialized();

  // solve the problem
  int itrs(0), max_itrs(1);
  double T1 = plist->get<double>("end time", 100.0);
  double dT = plist->get<double>("initial timestep", 1.0);
  double T(0.0);

  while (T < T1 && itrs < max_itrs) {
    MPK->AdvanceStep(T, T + dT);
    MPK->CommitStep(T, T + dT, Tags::DEFAULT);

    T += dT;
    dT *= 1.2;
    itrs++;
  }

  // update mesh
  auto mesh_tmp = Teuchos::rcp_const_cast<AmanziMesh::Mesh>(mesh);
  Teuchos::RCP<Mesh> mesh_vis = meshfactory.create("test/mechanics_cook.exo");
  mesh_tmp->setVisMesh(mesh_vis);

  int nnodes = mesh->getNumEntities(Entity_kind::NODE, AmanziMesh::Parallel_kind::OWNED);
  double scale(1.0);
  const auto& u = *S->Get<CompositeVector>("displacement").ViewComponent("node");

  double ymax(0.0);
  for (int n = 0; n < nnodes; ++n) {
    auto xp = mesh->getNodeCoordinate(n);
    for (int k = 0; k < 2; ++k) xp[k] += u[k][n] * scale;
    mesh_vis->setNodeCoordinate(n, xp);
    ymax = std::max(ymax, u[1][n]);
  }
  CHECK_CLOSE(0.032, ymax, 0.0034);

  // initialize I/O
  Teuchos::ParameterList iolist;
  iolist.get<std::string>("file name base", "plot");
  OutputXDMF io(iolist, mesh, true, false);

  io.InitializeCycle(T, 1, "");
  io.WriteVector(*u(0), "displacement-x", AmanziMesh::NODE);
  io.WriteVector(*u(1), "displacement-y", AmanziMesh::NODE);
  io.FinalizeCycle();

  // summary
  WriteStateStatistics(*S);
}
