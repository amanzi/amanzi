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
#include "LeastSquare.hh"
#include "GMVMesh.hh"
#include "Mesh.hh"
#include "MeshFactory.hh"

#include "ShallowWater_PK.hh"

#include "OutputXDMF.hh"

//--------------------------------------------------------------
// Initial conditions
//--------------------------------------------------------------
void hump_setIC(Teuchos::RCP<const Amanzi::AmanziMesh::Mesh> mesh,
                Teuchos::RCP<Amanzi::State>& S)
{
  double H_inf = 0.5;

  int ncells_owned = mesh->num_entities(Amanzi::AmanziMesh::CELL, Amanzi::AmanziMesh::Parallel_type::OWNED);
  std::string passwd = "state";

  Epetra_MultiVector& B_c = *S->GetFieldData("surface-bathymetry", passwd)->ViewComponent("cell");
  Epetra_MultiVector& h_c = *S->GetFieldData("surface-ponded_depth", passwd)->ViewComponent("cell");
  Epetra_MultiVector& ht_c = *S->GetFieldData("surface-total_depth", passwd)->ViewComponent("cell");
  Epetra_MultiVector& vel_c = *S->GetFieldData("surface-velocity", passwd)->ViewComponent("cell");
  Epetra_MultiVector& q_c = *S->GetFieldData("surface-discharge", "surface-discharge")->ViewComponent("cell");

  for (int c = 0; c < ncells_owned; c++) {
    const Amanzi::AmanziGeometry::Point& xc = mesh->cell_centroid(c);

    B_c[0][c] = std::max(0.0, 0.25 - 5 * ((xc[0] - 0.5) * (xc[0] - 0.5) + (xc[1] - 0.5) * (xc[1] - 0.5)));
    ht_c[0][c] = H_inf;
    h_c[0][c] = H_inf - B_c[0][c];

    vel_c[0][c] = 0.0;
    vel_c[1][c] = 0.0;
    q_c[0][c] = 0.0;
    q_c[1][c] = 0.0;
  }
}


/* **************************************************************** */
TEST(SHALLOW_WATER_BATHYMETRY) {
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
  auto gm = Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(2, regions_list, *comm));

  MeshFactory meshfactory(comm, gm);
  meshfactory.set_preference(Preference({Framework::MSTK}));
  RCP<Mesh> mesh = meshfactory.create(0.0, 0.0, 1.0, 1.0, 40, 40, true, false);

  // create a state
  Teuchos::ParameterList state_list = plist->sublist("state");
  RCP<State> S = rcp(new State(state_list));
  S->RegisterMesh("surface", mesh);
  S->set_time(0.0);

  Teuchos::RCP<TreeVector> soln = Teuchos::rcp(new TreeVector());
  Teuchos::ParameterList pk_tree = plist->sublist("PK tree").sublist("shallow water");

  // create a shallow water PK
  ShallowWater_PK SWPK(pk_tree,plist,S,soln);
  SWPK.Setup(S.ptr());
  S->Setup();
  S->InitializeFields();
  S->InitializeEvaluators();
  SWPK.Initialize(S.ptr());

  hump_setIC(mesh, S);

  const Epetra_MultiVector& hh = *S->GetFieldData("surface-ponded_depth")->ViewComponent("cell");
  const Epetra_MultiVector& ht = *S->GetFieldData("surface-total_depth")->ViewComponent("cell");
  const Epetra_MultiVector& vel = *S->GetFieldData("surface-velocity")->ViewComponent("cell");
  const Epetra_MultiVector& B = *S->GetFieldData("surface-bathymetry")->ViewComponent("cell");

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

  while (t_new < 0.1) {
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
}
