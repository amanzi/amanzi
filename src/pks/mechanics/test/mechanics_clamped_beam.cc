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
TEST(CLAMPED_BEAM)
{
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Mechanics;

  Comm_ptr_type comm = Amanzi::getDefaultComm();
  int MyPID = comm->MyPID();
  if (MyPID == 0) std::cout << "Test: clamped beam" << std::endl;

  // read parameter list
  std::string xmlFileName = "test/mechanics_clamped_beam.xml";
  Teuchos::RCP<Teuchos::ParameterList> plist = Teuchos::getParametersFromXmlFile(xmlFileName);

  // create a mesh framework
  Teuchos::ParameterList regions_list = plist->sublist("regions");
  auto gm = Teuchos::rcp(new GeometricModel(3, regions_list, *comm));

  int nx(100), ny(4);
  double L(25.0), W(1.3);
  MeshFactory meshfactory(comm, gm);
  meshfactory.set_preference(Preference({ Framework::MSTK }));
  Teuchos::RCP<const Mesh> mesh = meshfactory.create(0.0, 0.0, 0.0, L, W, W, nx, ny, ny);

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
  double dT = plist->get<double>("initial timestep", 1.0);
  double T(0.0);

  while (T < T1 && itrs < max_itrs) {
    EPK->AdvanceStep(T, T + dT);
    EPK->CommitStep(T, T + dT, Tags::DEFAULT);

    T += dT;
    dT *= 1.2;
    itrs++;
  }

  // update mesh
  auto mesh_tmp = Teuchos::rcp_const_cast<AmanziMesh::Mesh>(mesh);
  Teuchos::RCP<Mesh> mesh_vis = meshfactory.create(0.0, 0.0, 0.0, L, W, W, nx, ny, ny);
  mesh_tmp->setVisMesh(mesh_vis);

  int nnodes = mesh->getNumEntities(Entity_kind::NODE, AmanziMesh::Parallel_kind::OWNED);
  double scale(20.0);
  const auto& u = *S->Get<CompositeVector>("displacement").ViewComponent("node");
  const auto& p = *S->Get<CompositeVector>("hydrostatic_stress").ViewComponent("cell");

  for (int n = 0; n < nnodes; ++n) {
    auto xp = mesh->getNodeCoordinate(n);
    for (int k = 0; k < 3; ++k) xp[k] += u[k][n] * scale;
    mesh_vis->setNodeCoordinate(n, xp);
  }

  // initialize I/O
  Teuchos::ParameterList iolist;
  iolist.get<std::string>("file name base", "plot");
  OutputXDMF io(iolist, mesh, true, false);

  io.InitializeCycle(T, 1, "");
  io.WriteVector(*u(0), "displacement_x", AmanziMesh::NODE);
  io.WriteVector(*u(1), "displacement_y", AmanziMesh::NODE);
  io.WriteVector(*u(2), "displacement_z", AmanziMesh::NODE);
  io.WriteVector(*p(0), "hydrostatic_stress", AmanziMesh::CELL);
  io.FinalizeCycle();

  // summary
  WriteStateStatistics(*S);

  // exact solution is quadratic
  const auto& E = *S->Get<CompositeVector>("young_modulus").ViewComponent("cell");
  const auto& rho = *S->Get<CompositeVector>("particle_density").ViewComponent("cell");
  const auto& gravity = S->Get<AmanziGeometry::Point>("gravity");

  double f0 = rho[0][0] * gravity[2] * W * W;
  double I = std::pow(W, 4) / 12;
  double umax = 1.5 * rho[0][0] * std::fabs(gravity[2]) * std::pow(L, 4) / E[0][0] / W / W;

  double err(0.0);
  for (int n = 0; n < nx + 1; ++n) {
    double x = (mesh->getNodeCoordinate(n))[0];
    double exact = f0 / (24 * E[0][0] * I) * x * x * (x * x + 6 * L * L - 4 * L * x);
    // if (n < 101) std::cout << x << " " << u[2][n] << " " << exact << " " << err << std::endl;
    err = std::max(err, std::fabs(u[2][n] - exact) / umax);
    CHECK_CLOSE(u[2][n], exact, 6e-2 * umax);
  }

  std::cout << "Maximum error " << err * 100 << " %\n";
}
