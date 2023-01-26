/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Svetlana Tokareva (tokareva@lanl.gov)
*/

/*
 Shallow water PK

 */

// TPLs
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "UnitTest++.h"

// Amanzi
#include "GMVMesh.hh"
#include "Mesh.hh"
#include "MeshFactory.hh"

#include "ShallowWater_PK.hh"

#include "OutputXDMF.hh"

void
dam_break_1D_exact(double hL, double x0, double t, double x, double& h, double& u)
{
  double g = 9.81;
  double xA = x0 - t * std::sqrt(g * hL);
  double xB = x0 + 2. * t * std::sqrt(g * hL);

  //  std::cout << "xA = " << xA << ", xB = " << xB << std::endl;

  if (0. <= x && x < xA) {
    h = hL;
    u = 0.;
  } else {
    if (xA <= x && x < xB) {
      h = 4. / (9. * g) * (std::sqrt(g * hL) - (x - x0) / (2. * t)) *
          (std::sqrt(g * hL) - (x - x0) / (2. * t));
      u = 2. / 3. * ((x - x0) / t + std::sqrt(g * hL));
    } else {
      h = 0.;
      u = 0.;
    }
  }
}

void
dam_break_1D_exact_field(Teuchos::RCP<const Amanzi::AmanziMesh::Mesh> mesh,
                         Epetra_MultiVector& hh_ex,
                         Epetra_MultiVector& vx_ex,
                         double t)
{
  double hL, x0, x, h, u;

  hL = 10.;
  x0 = 1000.;

  int ncells_owned = mesh->getNumEntities(Amanzi::AmanziMesh::Entity_kind::CELL,
                                          Amanzi::AmanziMesh::Parallel_type::OWNED);

  for (int c = 0; c < ncells_owned; c++) {
    Amanzi::AmanziGeometry::Point xc = mesh->getCellCentroid(c);

    x = xc[0];

    dam_break_1D_exact(hL, x0, t, x, h, u);
    hh_ex[0][c] = h;
    vx_ex[0][c] = u;
  }
}

void
error(Teuchos::RCP<const Amanzi::AmanziMesh::Mesh> mesh,
      Epetra_MultiVector& hh_ex,
      Epetra_MultiVector& vx_ex,
      const Epetra_MultiVector& hh,
      const Epetra_MultiVector& vx)
{
  int ncells_owned = mesh->getNumEntities(Amanzi::AmanziMesh::Entity_kind::CELL,
                                          Amanzi::AmanziMesh::Parallel_type::OWNED);

  double err_max, err_L1;

  err_max = 0.;
  err_L1 = 0.;

  for (int c = 0; c < ncells_owned; c++) {
    err_max = std::max(err_max, std::abs(hh_ex[0][c] - hh[0][c]));
    err_L1 += std::abs(hh_ex[0][c] - hh[0][c]) * mesh->getCellVolume(c);
  }

  std::cout << "err_max = " << err_max << std::endl;
  std::cout << "err_L1  = " << err_L1 << std::endl;
}


/* **************************************************************** */
TEST(SHALLOW_WATER_1D_CONVERGENCE)
{
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::ShallowWater;

  Comm_ptr_type comm = Amanzi::getDefaultComm();
  int MyPID = comm->MyPID();
  if (MyPID == 0) std::cout << "Test: 1D shallow water" << std::endl;

  // read parameter list
  std::string xmlFileName = "test/shallow_water_1D.xml";
  Teuchos::RCP<Teuchos::ParameterList> plist = Teuchos::getParametersFromXmlFile(xmlFileName);

  /* create a mesh framework */
  auto gm = Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(2));
  if (MyPID == 0) std::cout << "Geometric model created." << std::endl;

  // create a mesh
  bool request_faces = true, request_edges = true;
  MeshFactory meshfactory(comm, gm);
  meshfactory.set_preference(Preference({ Framework::MSTK }));
  if (MyPID == 0) std::cout << "Mesh factory created." << std::endl;

  RCP<const Mesh> mesh;
  mesh = meshfactory.create(0.0, 0.0, 10.0, 1.0, 100, 1, request_faces, request_edges);
  // mesh = meshfactory.create("test/median63x64.exo",request_faces,request_edges); // works only with first order, no reconstruction
  if (MyPID == 0) std::cout << "Mesh created." << std::endl;

  // create a state
  RCP<State> S = rcp(new State());
  S->RegisterMesh("surface", rcp_const_cast<Mesh>(mesh));
  S->set_time(0.0);
  if (MyPID == 0) std::cout << "State created." << std::endl;

  Teuchos::RCP<TreeVector> soln = Teuchos::rcp(new TreeVector());

  Teuchos::ParameterList pk_tree = plist->sublist("PK tree").sublist("shallow water");

  // create a shallow water PK
  ShallowWater_PK SWPK(pk_tree, plist, S, soln);
  SWPK.Setup(S.ptr());
  S->Setup();

  S->InitializeFields();
  S->InitializeEvaluators();

  SWPK.Initialize(S.ptr()); // WE NEED A TEMPLATE OR SOMETHING TO INITIALIZE SWPK FROM THE TEST !!!

  if (MyPID == 0) std::cout << "Shallow water PK created." << std::endl;

  // create screen io
  auto vo = Teuchos::rcp(new Amanzi::VerboseObject("ShallowWater", *plist));
  WriteStateStatistics(S.ptr(), vo);

  // advance in time
  double t_old(0.0), t_new(0.0), dt;
  // Teuchos::RCP<Epetra_MultiVector>
  // tcc = S->GetFieldData("total_component_concentration", passwd)->ViewComponent("cell", false);

  // initialize io
  Teuchos::ParameterList iolist;
  std::string fname;
  fname = "SW_sol";
  iolist.get<std::string>("file name base", fname);
  OutputXDMF io(iolist, mesh, true, false);

  std::string passwd("state");

  auto& h_vec = *S->GetFieldData("surface-ponded_depth", passwd)->ViewComponent("cell");
  auto& u_vec = *S->GetFieldData("surface-velocity-x", passwd)->ViewComponent("cell");
  auto& v_vec = *S->GetFieldData("surface-velocity-y", passwd)->ViewComponent("cell");

  int iter = 0;
  bool flag = true;

  while (t_new < 0.1) {
    // cycle 1, time t
    double t_out = t_new;

    const Epetra_MultiVector& hh =
      *S->GetFieldData("surface-ponded_depth", passwd)->ViewComponent("cell");
    const Epetra_MultiVector& ht =
      *S->GetFieldData("surface-total_depth", passwd)->ViewComponent("cell");
    const Epetra_MultiVector& vx =
      *S->GetFieldData("surface-velocity-x", passwd)->ViewComponent("cell");
    const Epetra_MultiVector& vy =
      *S->GetFieldData("surface-velocity-y", passwd)->ViewComponent("cell");
    const Epetra_MultiVector& qx =
      *S->GetFieldData("surface-discharge-x", passwd)->ViewComponent("cell");
    const Epetra_MultiVector& qy =
      *S->GetFieldData("surface-discharge-y", passwd)->ViewComponent("cell");
    const Epetra_MultiVector& B =
      *S->GetFieldData("surface-bathymetry", passwd)->ViewComponent("cell");
    const Epetra_MultiVector& pid = *S->GetFieldData("surface-PID", passwd)->ViewComponent("cell");

    Epetra_MultiVector hh_ex(hh);
    Epetra_MultiVector vx_ex(vx);

    dam_break_1D_exact_field(mesh, hh_ex, vx_ex, t_out);

    io.InitializeCycle(t_out, iter);
    io.WriteVector(*hh(0), "depth", AmanziMesh::Entity_kind::CELL);
    io.WriteVector(*ht(0), "total_depth", AmanziMesh::Entity_kind::CELL);
    io.WriteVector(*vx(0), "vx", AmanziMesh::Entity_kind::CELL);
    io.WriteVector(*vy(0), "vy", AmanziMesh::Entity_kind::CELL);
    io.WriteVector(*qx(0), "qx", AmanziMesh::Entity_kind::CELL);
    io.WriteVector(*qy(0), "qy", AmanziMesh::Entity_kind::CELL);
    io.WriteVector(*B(0), "B", AmanziMesh::Entity_kind::CELL);
    io.WriteVector(*pid(0), "pid", AmanziMesh::Entity_kind::CELL);
    io.WriteVector(*hh_ex(0), "hh_ex", AmanziMesh::Entity_kind::CELL);
    io.WriteVector(*vx_ex(0), "vx_ex", AmanziMesh::Entity_kind::CELL);
    io.FinalizeCycle();

    dt = SWPK.get_dt();

    if (iter < 10) dt = 0.01 * dt;

    t_new = t_old + dt;

    std::cout << "h_vec.MyLength() = " << h_vec.MyLength() << std::endl;


    SWPK.AdvanceStep(t_old, t_new);
    SWPK.CommitStep(t_old, t_new, S);

    t_old = t_new;
    iter++;
  }

  if (MyPID == 0) std::cout << "Time-stepping finished." << std::endl;

  std::cout << "MyPID = " << MyPID << ", iter = " << iter << std::endl;

  // cycle 1, time t
  double t_out = t_new;

  const Epetra_MultiVector& hh =
    *S->GetFieldData("surface-ponded_depth", passwd)->ViewComponent("cell");
  const Epetra_MultiVector& ht =
    *S->GetFieldData("surface-total_depth", passwd)->ViewComponent("cell");
  const Epetra_MultiVector& vx =
    *S->GetFieldData("surface-velocity-x", passwd)->ViewComponent("cell");
  const Epetra_MultiVector& vy =
    *S->GetFieldData("surface-velocity-y", passwd)->ViewComponent("cell");
  const Epetra_MultiVector& qx =
    *S->GetFieldData("surface-discharge-x", passwd)->ViewComponent("cell");
  const Epetra_MultiVector& qy =
    *S->GetFieldData("surface-discharge-y", passwd)->ViewComponent("cell");
  const Epetra_MultiVector& B =
    *S->GetFieldData("surface-bathymetry", passwd)->ViewComponent("cell");
  const Epetra_MultiVector& pid = *S->GetFieldData("surface-PID", passwd)->ViewComponent("cell");

  Epetra_MultiVector hh_ex(hh);
  Epetra_MultiVector vx_ex(vx);

  dam_break_1D_exact_field(mesh, hh_ex, vx_ex, t_out);

  error(mesh, hh_ex, vx_ex, hh, vx);

  io.InitializeCycle(t_out, iter);
  io.WriteVector(*hh(0), "depth", AmanziMesh::Entity_kind::CELL);
  io.WriteVector(*ht(0), "total_depth", AmanziMesh::Entity_kind::CELL);
  io.WriteVector(*vx(0), "vx", AmanziMesh::Entity_kind::CELL);
  io.WriteVector(*vy(0), "vy", AmanziMesh::Entity_kind::CELL);
  io.WriteVector(*qx(0), "qx", AmanziMesh::Entity_kind::CELL);
  io.WriteVector(*qy(0), "qy", AmanziMesh::Entity_kind::CELL);
  io.WriteVector(*B(0), "B", AmanziMesh::Entity_kind::CELL);
  io.WriteVector(*pid(0), "pid", AmanziMesh::Entity_kind::CELL);
  io.WriteVector(*hh_ex(0), "hh_ex", AmanziMesh::Entity_kind::CELL);
  io.WriteVector(*vx_ex(0), "vx_ex", AmanziMesh::Entity_kind::CELL);
  io.FinalizeCycle();
}
