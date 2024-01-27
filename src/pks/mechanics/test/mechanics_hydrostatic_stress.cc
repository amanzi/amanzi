/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Navier Stokes PK

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
TEST(HYDROSTATIC_STRESS)
{
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Mechanics;

  Comm_ptr_type comm = Amanzi::getDefaultComm();
  int MyPID = comm->MyPID();
  if (MyPID == 0) std::cout << "Test: hydrostatic stress" << std::endl;

  // read parameter list
  std::string xmlFileName = "test/mechanics_hydrostatic_stress.xml";
  Teuchos::RCP<Teuchos::ParameterList> plist = Teuchos::getParametersFromXmlFile(xmlFileName);

  // create a mesh framework
  Teuchos::ParameterList regions_list = plist->sublist("regions");
  auto gm = Teuchos::rcp(new GeometricModel(3, regions_list, *comm));

  MeshFactory meshfactory(comm, gm);
  meshfactory.set_preference(Preference({ Framework::MSTK }));
  Teuchos::RCP<const Mesh> mesh = meshfactory.create(0.0, 0.0, 0.0, 4.0, 1.0, 25.0, 16, 4, 25);

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
  double dT = plist->get<double>("initial time step", 1.0);

  Teuchos::RCP<TreeVector> udot = Teuchos::rcp(new TreeVector(*soln));
  udot->PutScalar(0.0);
  MPK->bdf1_dae()->SetInitialState(0.0, soln, udot);
  MPK->UpdatePreconditioner(0.0, soln, dT);

  MPK->bdf1_dae()->TimeStep(dT, dT, soln);
  MPK->bdf1_dae()->CommitSolution(dT, soln);

  MPK->CommitStep(0.0, dT, Tags::DEFAULT);

  // initialize I/O
  const auto& u = *S->Get<CompositeVector>("displacement").ViewComponent("node");
  const auto& p = *S->Get<CompositeVector>("hydrostatic_stress").ViewComponent("cell");

  Teuchos::ParameterList iolist;
  iolist.get<std::string>("file name base", "plot");
  OutputXDMF io(iolist, mesh, true, false);

  io.InitializeCycle(dT, 1, "");
  io.WriteVector(*u(0), "displacement_x", AmanziMesh::NODE);
  io.WriteVector(*u(1), "displacement_y", AmanziMesh::NODE);
  io.WriteVector(*u(2), "displacement_z", AmanziMesh::NODE);
  io.WriteVector(*p(0), "hydrostatic_stress", AmanziMesh::CELL);
  io.FinalizeCycle();

  // summary
  WriteStateStatistics(*S);

  // check
  int ncells_owned =
    mesh->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);
  for (int c = 0; c < ncells_owned; ++c) {
    const auto& xp = mesh->getCellCentroid(c);
    if (xp[0] < 0.135 && xp[2] < 0.51) CHECK_CLOSE(-p[0][c], 200000.0, 50000.0);
    if (xp[0] > 3.871 && xp[2] < 0.51) CHECK_CLOSE(-p[0][c], 490000.0, 40000.0);
  }
}
