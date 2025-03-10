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
#include "MeshExtractedManifold.hh"
#include "MeshFactory.hh"
#include "OutputXDMF.hh"

#include "MechanicsFracturedMatrix_PK.hh"

/* **************************************************************** */
TEST(MECHANICS_FRACTURED_MATRIX)
{
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Mechanics;

  Comm_ptr_type comm = Amanzi::getDefaultComm();
  int MyPID = comm->MyPID();
  if (MyPID == 0) std::cout << "Test: Elastic deformation of fracture matrix" << std::endl;

  // read parameter list
  std::string xmlFileName = "test/mechanics_fractured_matrix.xml";
  Teuchos::RCP<Teuchos::ParameterList> plist = Teuchos::getParametersFromXmlFile(xmlFileName);

  // create a mesh framework
  Teuchos::ParameterList regions_list = plist->sublist("regions");
  auto gm = Teuchos::rcp(new GeometricModel(3, regions_list, *comm));

  auto mesh_list = Teuchos::rcp(new Teuchos::ParameterList());
  mesh_list->set<bool>("request edges", true);
  mesh_list->set<bool>("request faces", true);

  MeshFactory factory(comm, gm, mesh_list);
  factory.set_preference(Preference({ Framework::MSTK }));
  auto mesh = factory.create(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 6, 6, 6);

  // derive a fracture mesh
  std::vector<std::string> names;
  names.push_back("fracture");
  auto mesh_fracture_fw = Teuchos::rcp(new AmanziMesh::MeshExtractedManifold(
    mesh, names[0], AmanziMesh::Entity_kind::FACE, mesh->getComm(), gm, mesh_list));
  auto mesh_fracture = Teuchos::rcp(new AmanziMesh::Mesh(
    mesh_fracture_fw, Teuchos::rcp(new AmanziMesh::MeshAlgorithms()), mesh_list));

  // create a simple state and populate it
  Teuchos::ParameterList state_list = plist->sublist("state");
  Teuchos::RCP<State> S = Teuchos::rcp(new State(state_list));
  S->RegisterMesh("domain", mesh);
  S->RegisterMesh("fracture", mesh_fracture);

  Teuchos::RCP<TreeVector> soln = Teuchos::rcp(new TreeVector());
  Teuchos::ParameterList pk_tree = plist->sublist("PK tree").sublist("mechanics matrix fracture");
  auto MPK = Teuchos::rcp(new MechanicsFracturedMatrix_PK(pk_tree, plist, S, soln));
  MPK->parseParameterList();

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

  // initialize I/O
  Teuchos::ParameterList iolist;
  iolist.get<std::string>("file name base", "plot");
  OutputXDMF io(iolist, mesh_fracture, true, false);

  const auto& a = *S->Get<CompositeVector>("fracture-aperture").ViewComponent("cell");
  const auto& u = S->Get<CompositeVector>("displacement");
  const auto& nmap = *u.Map().Map("node", true);
  int nnodes =
    mesh->getNumEntities(AmanziMesh::Entity_kind::NODE, AmanziMesh::Parallel_kind::OWNED);

  Epetra_MultiVector u_long = *u.ViewComponent("node", true);
  Epetra_MultiVector u_short(mesh->getMap(AmanziMesh::Entity_kind::NODE, true), 3);

  for (int v = 0; v < nnodes; ++v) {
    int g = nmap.FirstPointInElement(v);
    for (int k = 0; k < 3; ++k) u_short[k][v] = u_long[k][g];
  }

  io.InitializeCycle(T, 1, "");
  io.WriteVector(*a(0), "fracture-aperture", AmanziMesh::CELL);
  // io.WriteVector(*u_short(0), "displacmenet-x", AmanziMesh::NODE);
  // io.WriteVector(*s(0), "shear_strain", AmanziMesh::CELL);
  io.FinalizeCycle();

  // summary
  WriteStateStatistics(*S);
}
