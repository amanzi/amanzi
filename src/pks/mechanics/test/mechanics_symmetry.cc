/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Mechanics PK

  Symmetry tests
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
TEST(SYMMETRY_2D)
{
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Mechanics;

  Comm_ptr_type comm = Amanzi::getDefaultComm();
  int MyPID = comm->MyPID();
  if (MyPID == 0) std::cout << "Test: Symmetric deformation" << std::endl;

  // read parameter list
  std::string xmlFileName = "test/mechanics_symmetry.xml";
  Teuchos::RCP<Teuchos::ParameterList> plist = Teuchos::getParametersFromXmlFile(xmlFileName);

  // create a mesh framework
  Teuchos::ParameterList regions_list = plist->sublist("regions");
  auto gm = Teuchos::rcp(new GeometricModel(2, regions_list, *comm));

  int nx(20), ny(10);
  MeshFactory meshfactory(comm, gm);
  meshfactory.set_preference(Preference({ Framework::MSTK }));
  Teuchos::RCP<const Mesh> mesh = meshfactory.create(0.0, 0.0, 2.0, 0.5, nx, ny);

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
  double dT(1.0);
  EPK->AdvanceStep(0.0, dT);
  EPK->CommitStep(0.0, dT, Tags::DEFAULT);

  // initialize I/O
  Teuchos::ParameterList iolist;
  iolist.get<std::string>("file name base", "plot");
  OutputXDMF io(iolist, mesh, true, false);

  io.InitializeCycle(0.0, 1, "");
  const auto& u = *S->Get<CompositeVector>("displacement").ViewComponent("node");
  const auto& e = *S->Get<CompositeVector>("volumetric_strain").ViewComponent("cell");
  io.WriteVector(*u(0), "displacement-x", AmanziMesh::NODE);
  io.WriteVector(*u(1), "displacement-y", AmanziMesh::NODE);
  io.WriteVector(*e(0), "volumetric_strain", AmanziMesh::CELL);
  io.FinalizeCycle();

  // summary
  WriteStateStatistics(*S);

  int n(0), n0;
  for (int j = 0; j < ny + 1; ++j) {
    n0 = n;
    for (int i = 0; i < nx + 1; ++i) {
      CHECK_CLOSE(u[0][n], u[0][n0], 1e-12);
      CHECK_CLOSE(u[1][n], u[1][n0], 1e-12);
      n++;
    }
  }
}
