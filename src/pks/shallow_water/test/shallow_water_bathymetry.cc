/*
  Shallow water PK
 
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.
 
  Author: Svetlana Tokareva (tokareva@lanl.gov)
*/

// TPLs
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "UnitTest++.h"

// Amanzi
#include "IO.hh"
#include "LeastSquare.hh"
#include "Mesh.hh"
#include "MeshFactory.hh"
#include "OutputXDMF.hh"

#include "ShallowWater_PK.hh"


/* **************************************************************** */
Epetra_MultiVector RunTest(int ntest) {
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::ShallowWater;

  Comm_ptr_type comm = Amanzi::getDefaultComm();
  int MyPID = comm->MyPID();
  if (MyPID == 0) std::cout << "Test: 2D shallow water: equilibrium solution with a hump" << std::endl;

  // read parameter list
  std::string xmlFileName = "test/shallow_water_bathymetry.xml";
  Teuchos::RCP<Teuchos::ParameterList> plist = Teuchos::getParametersFromXmlFile(xmlFileName);

  // create a mesh
  ParameterList regions_list = plist->get<Teuchos::ParameterList>("regions");
  auto gm = Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(3, regions_list, *comm));

  MeshFactory meshfactory(comm, gm);
  meshfactory.set_preference(Preference({Framework::MSTK}));
  RCP<Mesh> mesh;
  if (ntest == 1) {
    mesh = meshfactory.create(0.0, 0.0, 100.0, 100.0, 20, 20, true, true);
  } else {
    RCP<Mesh> mesh3D = meshfactory.create(0.0, 0.0, 0.0, 100.0, 100.0, 10.0, 20, 20, 4, true, true);
    mesh = meshfactory.create(mesh3D, { "TopSurface" }, AmanziMesh::FACE, true, true, true);
  }

  // create a state
  Teuchos::ParameterList state_list = plist->sublist("state");
  RCP<State> S = rcp(new State(state_list));
  S->RegisterMesh("surface", mesh);

  Teuchos::RCP<TreeVector> soln = Teuchos::rcp(new TreeVector());
  Teuchos::ParameterList pk_tree = plist->sublist("PK tree").sublist("shallow water");

  // create a shallow water PK
  ShallowWater_PK SWPK(pk_tree,plist,S,soln);
  SWPK.Setup(S.ptr());
  S->Setup();
  S->InitializeFields();
  S->InitializeEvaluators();
  S->set_time(0.0);
  SWPK.Initialize(S.ptr());
  S->CheckAllFieldsInitialized();

  const auto& hh = *S->Get<CompositeVector>("surface-ponded_depth").ViewComponent("cell");
  const auto& ht = *S->Get<CompositeVector>("surface-total_depth").ViewComponent("cell");
  const auto& vel = *S->Get<CompositeVector>("surface-velocity").ViewComponent("cell");
  const auto& B = *S->Get<CompositeVector>("surface-bathymetry").ViewComponent("cell");

  // create screen io
  auto vo = Teuchos::rcp(new Amanzi::VerboseObject("ShallowWater", *plist));
  WriteStateStatistics(*S, *vo);

  // advance in time
  double t_old(0.0), t_new(0.0), dt;

  // initialize io
  Teuchos::ParameterList iolist;
  std::string fname;
  fname = "SW_sol";
  iolist.get<std::string>("file name base", fname);
  OutputXDMF io(iolist, mesh, true, false);

  std::string passwd("state");
  int iter = 0;

  while (t_new < 2.0 && iter < 20) {
    double t_out = t_new;

    if (iter % 5 == 0) {
      io.InitializeCycle(t_out, iter, "");
      io.WriteVector(*hh(0), "depth", AmanziMesh::CELL);
      io.WriteVector(*ht(0), "total_depth", AmanziMesh::CELL);
      io.WriteVector(*vel(0), "vx", AmanziMesh::CELL);
      io.WriteVector(*vel(1), "vy", AmanziMesh::CELL);
      io.WriteVector(*B(0), "B", AmanziMesh::CELL);
      io.FinalizeCycle();
    }

    dt = SWPK.get_dt();
    if (iter < 10) dt = 0.01 * dt;

    t_new = t_old + dt;

    SWPK.AdvanceStep(t_old, t_new);
    SWPK.CommitStep(t_old, t_new, S);

    t_old = t_new;
    iter++;
  }

  WriteStateStatistics(*S, *vo);
  
  return hh;
}

TEST(SHALLOW_WATER_BATHYMETRY) {
  auto fa = RunTest(1);
  auto fb = RunTest(2);

  double vala[1], valb[1];
  fa.MeanValue(vala);
  fb.MeanValue(valb);
  CHECK_CLOSE(vala[0], valb[0], 1e-10); 
}
