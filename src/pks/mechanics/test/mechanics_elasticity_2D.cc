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
TEST(ELASTIC_2D)
{
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Mechanics;

  Comm_ptr_type comm = Amanzi::getDefaultComm();
  int MyPID = comm->MyPID();
  if (MyPID == 0) std::cout << "Test: Elastic deformation in 2D" << std::endl;

  // read parameter list
  std::string xmlFileName = "test/mechanics_elasticity_2D.xml";
  Teuchos::RCP<Teuchos::ParameterList> plist = Teuchos::getParametersFromXmlFile(xmlFileName);

  // create a mesh framework
  Teuchos::ParameterList regions_list = plist->sublist("regions");
  auto gm = Teuchos::rcp(new GeometricModel(2, regions_list, *comm));

  MeshFactory meshfactory(comm, gm);
  meshfactory.set_preference(Preference({ Framework::MSTK }));
  Teuchos::RCP<const Mesh> mesh = meshfactory.create(0.0, 0.0, 1.0, 1.0, 20, 20);

  // create a simple state and populate it
  Teuchos::ParameterList state_list = plist->sublist("state");
  Teuchos::RCP<State> S = Teuchos::rcp(new State(state_list));
  S->RegisterDomainMesh(Teuchos::rcp_const_cast<Mesh>(mesh));

  Teuchos::RCP<TreeVector> soln = Teuchos::rcp(new TreeVector());
  Teuchos::ParameterList pk_tree = plist->sublist("PK tree").sublist("mechanics elasticity");
  auto EPK = Teuchos::rcp(new MechanicsElasticity_PK(pk_tree, plist, S, soln));

  EPK->Setup();
  S->Setup();
  S->InitializeFields();
  S->InitializeEvaluators();

  EPK->Initialize();
  S->CheckAllFieldsInitialized();

  // solve the problem
  int itrs(0);
  int max_itrs = plist->get<int>("max iterations", 50);
  double T1 = plist->get<double>("end time", 100.0);
  double dT = plist->get<double>("initial time step", 1.0);
  double T(0.0);

  while (T < T1 && itrs < max_itrs) {
    std::cout << "Cycle #" << itrs << std::endl;
    EPK->AdvanceStep(T, T + dT);
    EPK->CommitStep(T, T + dT, Tags::DEFAULT);
 
    T += dT;
    dT *= 1.2;
    itrs++;
  }

  // initialize I/O
  Teuchos::ParameterList iolist;
  iolist.get<std::string>("file name base", "plot");
  OutputXDMF io(iolist, mesh, true, false);

  io.InitializeCycle(T, 1, "");
  const auto& u = *S->Get<CompositeVector>("displacement").ViewComponent("node");
  const auto& p = *S->Get<CompositeVector>("hydrostatic_stress").ViewComponent("cell");
  io.WriteVector(*u(0), "displacement", AmanziMesh::NODE);
  io.WriteVector(*p(0), "hydrostatic_stress", AmanziMesh::CELL);
  io.FinalizeCycle();

  // summary
  WriteStateStatistics(*S);
  // S->Get<CompositeVector>("displacement", Tags::DEFAULT).Print(std::cout);

  int nnodes = mesh->getNumEntities(Entity_kind::NODE, AmanziMesh::Parallel_kind::OWNED);
  for (int n = 0; n < nnodes; ++n) {
    const auto& xv = mesh->getNodeCoordinate(n);
    CHECK_CLOSE(u[0][n], xv[0], 1e-6);
  }
}
