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
#include "IO.hh"
#include "LeastSquare.hh"
#include "Mesh.hh"
#include "MeshFactory.hh"
#include "OutputXDMF.hh"

// Amanzi::ShallowWater
#include "ShallowWater_PK.hh"
#include "ShallowWater_Helper.hh"

using namespace Amanzi;

//--------------------------------------------------------------
// Analytic solution
//--------------------------------------------------------------
double
hf(double x)
{
  return 2. * cos(x) + 2. * x * sin(x) + 1. / 8. * cos(2. * x) + x / 4. * sin(2. * x) +
         12. / 16. * x * x;
}

void
vortex_2D_exact(double t, double x, double y, double& h, double& u, double& v)
{
  double g = 9.81;
  double gmm = 15.;
  double omega = M_PI;
  double u_inf = 6., v_inf = 0.;
  double H_inf = 10.;
  double xc = 5., yc = 5.;

  double xt = x - xc - u_inf * t;
  double yt = y - yc - v_inf * t;

  double rc = std::sqrt(xt * xt + yt * yt); // distance from the vortex core

  // if (omega*rc <= M_PI) {
  if (rc <= 1.) {
    h = H_inf + 1. / g * (gmm / omega) * (gmm / omega) * (hf(omega * rc) - hf(M_PI));
    u = u_inf + gmm * (1. + cos(omega * rc)) * (yc - y);
    v = v_inf + gmm * (1. + cos(omega * rc)) * (x - xc);
  } else {
    h = H_inf;
    u = u_inf;
    v = v_inf;
  }
}


void
vortex_2D_exact_field(Teuchos::RCP<const Amanzi::AmanziMesh::Mesh> mesh,
                      Epetra_MultiVector& hh_ex,
                      Epetra_MultiVector& vel_ex,
                      double t)
{
  double x, y, h, u, v;

  int ncells_owned = mesh->getNumEntities(Amanzi::AmanziMesh::Entity_kind::CELL,
                                          Amanzi::AmanziMesh::Parallel_kind::OWNED);

  for (int c = 0; c < ncells_owned; c++) {
    const Amanzi::AmanziGeometry::Point& xc = mesh->getCellCentroid(c);

    x = xc[0];
    y = xc[1];

    vortex_2D_exact(t, x, y, h, u, v);
    hh_ex[0][c] = h;
    vel_ex[0][c] = u;
    vel_ex[1][c] = v;
  }
}


void
vortex_2D_setIC(Teuchos::RCP<const Amanzi::AmanziMesh::Mesh> mesh, Teuchos::RCP<Amanzi::State>& S)
{
  int ncells_owned = mesh->getNumEntities(Amanzi::AmanziMesh::Entity_kind::CELL,
                                          Amanzi::AmanziMesh::Parallel_kind::OWNED);

  std::string passwd("");

  const auto& B_vec_c = *S->Get<CompositeVector>("surface-bathymetry").ViewComponent("cell");
  auto& h_vec_c =
    *S->GetW<CompositeVector>("surface-ponded_depth", Tags::DEFAULT, passwd).ViewComponent("cell");
  auto& ht_vec_c =
    *S->GetW<CompositeVector>("surface-total_depth", Tags::DEFAULT, passwd).ViewComponent("cell");
  auto& vel_vec_c =
    *S->GetW<CompositeVector>("surface-velocity", Tags::DEFAULT, passwd).ViewComponent("cell");
  auto& q_vec_c = *S->GetW<CompositeVector>("surface-discharge", Tags::DEFAULT, "surface-discharge")
                     .ViewComponent("cell");

  for (int c = 0; c < ncells_owned; c++) {
    const Amanzi::AmanziGeometry::Point& xc = mesh->getCellCentroid(c);
    double h, u, v;
    vortex_2D_exact(0., xc[0], xc[1], h, u, v);
    h_vec_c[0][c] = h;
    ht_vec_c[0][c] = h_vec_c[0][c] + B_vec_c[0][c];
    vel_vec_c[0][c] = u;
    vel_vec_c[1][c] = v;
    q_vec_c[0][c] = u * h;
    q_vec_c[1][c] = v * h;
  }
}

void
error(Teuchos::RCP<const Amanzi::AmanziMesh::Mesh> mesh,
      Epetra_MultiVector& hh_ex,
      Epetra_MultiVector& vel_ex,
      const Epetra_MultiVector& hh,
      const Epetra_MultiVector& vel,
      double& err_max,
      double& err_L1,
      double& hmax)
{
  int ncells_owned = mesh->getNumEntities(Amanzi::AmanziMesh::Entity_kind::CELL,
                                          Amanzi::AmanziMesh::Parallel_kind::OWNED);

  err_max = 0.;
  err_L1 = 0.;
  hmax = 0.;

  for (int c = 0; c < ncells_owned; c++) {
    double tmp = std::abs(hh_ex[0][c] - hh[0][c]);
    err_max = std::max(err_max, tmp);
    err_L1 += tmp * mesh->getCellVolume(c);
    hmax = std::sqrt(mesh->getCellVolume(c));
  }

  double err_max_tmp;
  double err_L1_tmp;

  mesh->getComm()->MaxAll(&err_max, &err_max_tmp, 1);
  mesh->getComm()->SumAll(&err_L1, &err_L1_tmp, 1);

  err_max = err_max_tmp;
  err_L1 = err_L1_tmp;

  std::cout << "err_max = " << err_max << std::endl;
  std::cout << "err_L1  = " << err_L1 << std::endl;
}


/* **************************************************************** */
TEST(SHALLOW_WATER_2D_SMOOTH)
{
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::ShallowWater;

  Comm_ptr_type comm = Amanzi::getDefaultComm();
  int MyPID = comm->MyPID();
  if (MyPID == 0) std::cout << "Test: 2D shallow water: smooth vortex" << std::endl;

  // read parameter list
  std::string xmlFileName = "test/shallow_water_2D_smooth.xml";
  Teuchos::RCP<Teuchos::ParameterList> plist = Teuchos::getParametersFromXmlFile(xmlFileName);

  /* create a mesh framework */
  ParameterList regions_list = plist->get<Teuchos::ParameterList>("regions");
  auto gm = Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(2, regions_list, *comm));

  // create a mesh
  MeshFactory meshfactory(comm, gm);
  meshfactory.set_preference(Preference({ Framework::MSTK }));

  std::vector<double> dx, Linferror, L1error, L2error, dt_val;

  for (int NN = 20; NN <= 80; NN *= 2) {
    RCP<Mesh> mesh = meshfactory.create(0.0, 0.0, 10.0, 10.0, NN, NN);
    // mesh = meshfactory.create("test/median63x64.exo",request_faces,request_edges);
    // works only with first order, no reconstruction

    // create a state
    Teuchos::ParameterList state_list = plist->sublist("state");
    RCP<State> S = rcp(new State(state_list));
    S->RegisterMesh("surface", mesh);

    Teuchos::RCP<TreeVector> soln = Teuchos::rcp(new TreeVector());

    Teuchos::ParameterList pk_tree = plist->sublist("PKs").sublist("shallow water");

    // create a shallow water PK
    ShallowWater_PK SWPK(pk_tree, plist, S, soln);
    SWPK.Setup();
    S->Setup();
    S->InitializeFields();
    S->InitializeEvaluators();
    S->set_time(0.0);
    SWPK.Initialize();

    vortex_2D_setIC(mesh, S);
    S->CheckAllFieldsInitialized();

    const auto& hh = *S->Get<CompositeVector>("surface-ponded_depth").ViewComponent("cell");
    const auto& vel = *S->Get<CompositeVector>("surface-velocity").ViewComponent("cell");

    // create screen io
    auto vo = Teuchos::rcp(new Amanzi::VerboseObject("ShallowWater", pk_tree));
    WriteStateStatistics(*S, *vo);

    // advance in time
    double t_old(0.0), t_new(0.0), dt, dt_max(0.0);

    // initialize io
    Teuchos::ParameterList iolist;
    std::string fname;
    fname = "SW_sol";
    iolist.get<std::string>("file name base", fname);
    OutputXDMF io(iolist, mesh, true, false);

    std::string passwd("state");

    int iter = 0;

    while (t_new < 0.001) {
      // cycle 1, time t
      double t_out = t_new;

      Epetra_MultiVector hh_ex(hh);
      Epetra_MultiVector vel_ex(vel);

      vortex_2D_exact_field(mesh, hh_ex, vel_ex, t_out);

      if (iter % 100 == 0) { IO_Fields(t_out, iter, MyPID, io, *S, &hh_ex, &vel_ex); }

      dt = SWPK.get_dt();

      t_new = t_old + dt;
      dt_max = std::max(dt_max, dt);

      SWPK.AdvanceStep(t_old, t_new);
      SWPK.CommitStep(t_old, t_new, Tags::DEFAULT);

      t_old = t_new;
      iter++;
    }

    if (MyPID == 0) std::cout << "Time-stepping finished." << std::endl;

    // cycle 1, time t
    double t_out = t_new;

    Epetra_MultiVector hh_ex(hh);
    Epetra_MultiVector vel_ex(vel);

    vortex_2D_exact_field(mesh, hh_ex, vel_ex, t_out);

    double err_max, err_L1, hmax;

    error(mesh, hh_ex, vel_ex, hh, vel, err_max, err_L1, hmax);

    dx.push_back(hmax);
    Linferror.push_back(err_max);
    L1error.push_back(err_L1);
    dt_val.push_back(dt_max);

    IO_Fields(t_out, iter, MyPID, io, *S, &hh_ex, &vel_ex);
  } // NN

  double order = Amanzi::Utils::bestLSfit(dx, L1error);

  std::cout << "computed order L1 (dx) = " << order << std::endl;

  double order_Linf_dt = Amanzi::Utils::bestLSfit(dt_val, Linferror);
  std::cout << "computed order Linf (dt) = " << order_Linf_dt << std::endl;

  CHECK(order > 0.9); // first order scheme
}
