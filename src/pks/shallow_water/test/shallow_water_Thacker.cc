/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
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
#include "CompositeVector.hh"
#include "IO.hh"
#include "LeastSquare.hh"
#include "Mesh.hh"
#include "MeshFactory.hh"
#include "OutputXDMF.hh"

// Amanzi::ShallowWater
#include "ShallowWater_PK.hh"
#include "ShallowWater_Helper.hh"

//--------------------------------------------------------------
// Analytic solution (Thacker's solution [Beljadid et. al. 2016]
//--------------------------------------------------------------

using namespace Amanzi;

void
analytical_exact(double t, double x, double y, double& h, double& u, double& v)
{
  // Thacker's time dependent solution [Beljadid et. al. 2016]
  double eta = 2.0, T = 1.436739427831727;
  //  double R0 = T * std::sqrt(2.0*g*eta);
  double R0 = 9.0;

  h = eta * ((T * T) / (T * T + t * t)) *
      (1.0 - ((x * x + y * y) / (R0 * R0)) * ((T * T) / (t * t + T * T)));
  u = (x * t) / (t * t + T * T);
  v = (y * t) / (t * t + T * T);
}


void
analytical_exact_field(Teuchos::RCP<const AmanziMesh::Mesh> mesh,
                       Epetra_MultiVector& hh_ex,
                       Epetra_MultiVector& vel_ex,
                       double t)
{
  double x, y, h, u, v;

  int ncells_owned = mesh->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);

  for (int c = 0; c < ncells_owned; c++) {
    const AmanziGeometry::Point& xc = mesh->getCellCentroid(c);

    x = xc[0];
    y = xc[1];

    analytical_exact(t, x, y, h, u, v);
    hh_ex[0][c] = h;
    vel_ex[0][c] = u;
    vel_ex[1][c] = v;
  }
}


void
analytical_setIC(Teuchos::RCP<const AmanziMesh::Mesh> mesh, Teuchos::RCP<State>& S)
{
  int ncells_owned = mesh->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);

  std::string passwd("");

  auto& B_vec_c =
    *S->GetW<CompositeVector>("surface-bathymetry", Tags::DEFAULT, passwd).ViewComponent("cell");
  auto& h_vec_c =
    *S->GetW<CompositeVector>("surface-ponded_depth", Tags::DEFAULT, passwd).ViewComponent("cell");
  auto& ht_vec_c =
    *S->GetW<CompositeVector>("surface-total_depth", Tags::DEFAULT, passwd).ViewComponent("cell");
  auto& vel_vec_c =
    *S->GetW<CompositeVector>("surface-velocity", Tags::DEFAULT, passwd).ViewComponent("cell");
  auto& q_vec_c = *S->GetW<CompositeVector>("surface-discharge", Tags::DEFAULT, "surface-discharge")
                     .ViewComponent("cell");

  for (int c = 0; c < ncells_owned; c++) {
    const AmanziGeometry::Point& xc = mesh->getCellCentroid(c);
    double h, u, v;
    analytical_exact(0., xc[0], xc[1], h, u, v);
    h_vec_c[0][c] = h;
    B_vec_c[0][c] = 0.0;
    ht_vec_c[0][c] = h_vec_c[0][c] + B_vec_c[0][c];
    vel_vec_c[0][c] = u;
    vel_vec_c[1][c] = v;
    q_vec_c[0][c] = u * h;
    q_vec_c[1][c] = v * h;
  }
}

void
error(Teuchos::RCP<const AmanziMesh::Mesh> mesh,
      Epetra_MultiVector& hh_ex,
      Epetra_MultiVector& vel_ex,
      const Epetra_MultiVector& hh,
      const Epetra_MultiVector& vel,
      double& err_max,
      double& err_L1,
      double& hmax,
      double& err_max_u,
      double& err_L1_u,
      double& err_max_v,
      double& err_L1_v)
{
  int ncells_owned = mesh->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);

  err_max = 0.;
  err_L1 = 0.;
  hmax = 0.;
  err_max_u = 0.;
  err_L1_u = 0.;
  err_max_v = 0.;
  err_L1_v = 0.;

  for (int c = 0; c < ncells_owned; c++) {
    double tmp = std::abs(hh_ex[0][c] - hh[0][c]);
    err_max = std::max(err_max, tmp);
    err_L1 += tmp * mesh->getCellVolume(c);
    hmax = std::sqrt(mesh->getCellVolume(c));

    tmp = std::abs(vel_ex[0][c] - vel[0][c]);
    err_max_u = std::max(err_max_u, tmp);
    err_L1_u += tmp * mesh->getCellVolume(c);

    tmp = std::abs(vel_ex[1][c] - vel[1][c]);
    err_max_v = std::max(err_max_v, tmp);
    err_L1_v += tmp * mesh->getCellVolume(c);
  }

  double err_max_tmp;
  double err_L1_tmp;
  double err_max_tmp_u;
  double err_L1_tmp_u;
  double err_max_tmp_v;
  double err_L1_tmp_v;

  mesh->getComm()->MaxAll(&err_max, &err_max_tmp, 1);
  mesh->getComm()->SumAll(&err_L1, &err_L1_tmp, 1);
  mesh->getComm()->MaxAll(&err_max_u, &err_max_tmp_u, 1);
  mesh->getComm()->SumAll(&err_L1_u, &err_L1_tmp_u, 1);
  mesh->getComm()->MaxAll(&err_max_v, &err_max_tmp_v, 1);
  mesh->getComm()->SumAll(&err_L1_v, &err_L1_tmp_v, 1);

  err_max = err_max_tmp;
  err_L1 = err_L1_tmp;

  err_max_u = err_max_tmp_u;
  err_L1_u = err_L1_tmp_u;

  err_max_v = err_max_tmp_v;
  err_L1_v = err_L1_tmp_v;

  if (mesh->getComm()->MyPID() == 0) {
    std::cout << "err_max = " << err_max << std::endl;
    std::cout << "err_L1  = " << err_L1 << std::endl;

    std::cout << "err_max (u)= " << err_max_u << std::endl;
    std::cout << "err_L1 (u) = " << err_L1_u << std::endl;

    std::cout << "err_max (v)= " << err_max_v << std::endl;
    std::cout << "err_L1 (v) = " << err_L1_v << std::endl;
  }
}


/* **************************************************************** */
void
RunTest(int icase)
{
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::ShallowWater;

  Comm_ptr_type comm = Amanzi::getDefaultComm();
  int MyPID = comm->MyPID();
  if (MyPID == 0) std::cout << "Test: 2D shallow water: Thacker's solution" << std::endl;

  // read parameter list
  std::string xmlFileName = "test/shallow_water_Thacker.xml";
  Teuchos::RCP<Teuchos::ParameterList> plist = Teuchos::getParametersFromXmlFile(xmlFileName);

  // set temporal discretization order
  std::string temporal_order_string;
  if (icase == 1) {
    plist->sublist("PKs").sublist("shallow water").set<int>("temporal discretization order", 1);
    temporal_order_string = "RK1: Explicit_RK_Euler";
  } else if (icase == 2) {
    plist->sublist("PKs").sublist("shallow water").set<int>("temporal discretization order", 2);
    temporal_order_string = "RK2: Explicit_RK_Midpoint";
  } else if (icase == 3) {
    plist->sublist("PKs").sublist("shallow water").set<int>("temporal discretization order", 3);
    temporal_order_string = "RK3: Explicit_TVD_RK3";
  }

  // create a mesh framework
  ParameterList regions_list = plist->get<Teuchos::ParameterList>("regions");
  auto gm = Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(2, regions_list, *comm));

  // create a mesh
  MeshFactory meshfactory(comm, gm);
  meshfactory.set_preference(Preference({ Framework::MSTK }));

  std::vector<double> dx, Linferror, L1error, L2error;

  for (int NN = 40; NN <= 160; NN *= 2) {
    RCP<Mesh> mesh = meshfactory.create(-3.0, -3.0, 3.0, 3.0, NN, NN);
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

    analytical_setIC(mesh, S);
    S->CheckAllFieldsInitialized();

    const auto& hh = *S->Get<CompositeVector>("surface-ponded_depth").ViewComponent("cell");
    const auto& vel = *S->Get<CompositeVector>("surface-velocity").ViewComponent("cell");

    // create screen io
    auto vo = Teuchos::rcp(new Amanzi::VerboseObject("ShallowWater", pk_tree));
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

    while (t_new < 0.01) {
      double t_out = t_new;

      Epetra_MultiVector hh_ex(hh);
      Epetra_MultiVector vel_ex(vel);

      analytical_exact_field(mesh, hh_ex, vel_ex, t_out);

      if (iter % 5 == 0) { IO_Fields(t_out, iter, MyPID, io, *S, &hh_ex, &vel_ex); }

      dt = SWPK.get_dt();

      t_new = t_old + dt;

      SWPK.AdvanceStep(t_old, t_new);
      SWPK.CommitStep(t_old, t_new, Tags::DEFAULT);

      t_old = t_new;
      iter++;
    }

    if (MyPID == 0) std::cout << "Time-stepping finished, T=" << t_new << std::endl;

    // cycle 1, time t
    double t_out = t_new;

    Epetra_MultiVector hh_ex(hh);
    Epetra_MultiVector vel_ex(vel);

    analytical_exact_field(mesh, hh_ex, vel_ex, t_out);

    double err_max, err_L1, hmax, err_max_u, err_L1_u, err_max_v, err_L1_v;

    error(mesh,
          hh_ex,
          vel_ex,
          hh,
          vel,
          err_max,
          err_L1,
          hmax,
          err_max_u,
          err_L1_u,
          err_max_v,
          err_L1_v);

    dx.push_back(hmax);
    Linferror.push_back(err_max);
    L1error.push_back(err_L1);

    IO_Fields(t_out, iter, MyPID, io, *S, &hh_ex, &vel_ex);
  } // NN

  double L1_order = Amanzi::Utils::bestLSfit(dx, L1error);
  double Linf_order = Amanzi::Utils::bestLSfit(dx, Linferror);

  if (MyPID == 0) {
    std::cout << temporal_order_string << std::endl;
    std::cout << "computed order L1  = " << L1_order << std::endl;
    std::cout << "computed order Linf = " << Linf_order << std::endl;
  }

  if (icase == 1) {
    CHECK(L1_order > 0.9); // first order scheme (first order time stepping)
  } else {
    CHECK(L1_order > 1.8); //second order scheme (second/third order time stepping)
  }
}

TEST(SHALLOW_WATER_THACKER)
{
  RunTest(1); // RK1: Forward Euler
  RunTest(2); // RK2: Midpoint
  RunTest(3); // RK3: TVD 3rd Order
}
