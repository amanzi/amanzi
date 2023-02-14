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
#include "IO.hh"
#include "Mesh.hh"
#include "MeshFactory.hh"
#include "ShallowWater_PK.hh"
#include "OutputXDMF.hh"

// Amznzi::ShallowWater
#include "ShallowWater_Helper.hh"

// General
#define _USE_MATH_DEFINES
#include "math.h"

using namespace Amanzi;

double H_inf = 0.5; // Lake at rest total height

// ---- ---- ---- ---- ---- ---- ---- ----
// Exact Solution
// ---- ---- ---- ---- ---- ---- ---- ----
void
lake_at_rest_exact(double t, double x, double y, double& ht, double& u, double& v)
{
  ht = H_inf;
  u = 0;
  v = 0;
}

void
lake_at_rest_exact_field(Teuchos::RCP<const Amanzi::AmanziMesh::Mesh> mesh,
                         Epetra_MultiVector& ht_ex,
                         Epetra_MultiVector& vel_ex,
                         double t)
{
  double x, y, ht, u, v;

  int ncells_owned =
    mesh->getNumEntities(Amanzi::AmanziMesh::Entity_kind::CELL, Amanzi::AmanziMesh::Parallel_kind::OWNED);

  for (int c = 0; c < ncells_owned; ++c) {
    const Amanzi::AmanziGeometry::Point& xc = mesh->getCellCentroid(c);

    x = xc[0];
    y = xc[1];

    lake_at_rest_exact(t, x, y, ht, u, v);
    ht_ex[0][c] = ht;
    vel_ex[0][c] = u;
    vel_ex[1][c] = v;
  }
}


// ---- ---- ---- ---- ---- ---- ---- ----
// Inital Conditions
// ---- ---- ---- ---- ---- ---- ---- ----
void
lake_at_rest_setIC(Teuchos::RCP<const Amanzi::AmanziMesh::Mesh> mesh,
                   Teuchos::RCP<Amanzi::State>& S,
                   int icase)
{
  double pi = M_PI;
  const double rho = S->Get<double>("const_fluid_density");
  const double patm = S->Get<double>("atmospheric_pressure");

  double g = norm(S->Get<AmanziGeometry::Point>("gravity"));

  int ncells_wghost =
    mesh->getNumEntities(Amanzi::AmanziMesh::Entity_kind::CELL, Amanzi::AmanziMesh::Parallel_kind::ALL);
  int nnodes_wghost =
    mesh->getNumEntities(Amanzi::AmanziMesh::Entity_kind::NODE, Amanzi::AmanziMesh::Parallel_kind::ALL);
  std::string passwd("");

  auto& B_c =
    *S->GetW<CompositeVector>("surface-bathymetry", Tags::DEFAULT, passwd).ViewComponent("cell");
  auto& B_n =
    *S->GetW<CompositeVector>("surface-bathymetry", Tags::DEFAULT, passwd).ViewComponent("node");
  auto& h_c =
    *S->GetW<CompositeVector>("surface-ponded_depth", Tags::DEFAULT, passwd).ViewComponent("cell");
  auto& ht_c =
    *S->GetW<CompositeVector>("surface-total_depth", Tags::DEFAULT, passwd).ViewComponent("cell");
  auto& vel_c =
    *S->GetW<CompositeVector>("surface-velocity", Tags::DEFAULT, passwd).ViewComponent("cell");
  auto& q_c = *S->GetW<CompositeVector>("surface-discharge", Tags::DEFAULT, "surface-discharge")
                 .ViewComponent("cell");
  auto& p_c =
    *S->GetW<CompositeVector>("surface-ponded_pressure", Tags::DEFAULT, "surface-ponded_pressure")
       .ViewComponent("cell");

  // Define bathymetry at the cell vertices (Bn)
  for (int n = 0; n < nnodes_wghost; ++n) {
    Amanzi::AmanziGeometry::Point node_crd;

    node_crd = mesh->getNodeCoordinate(n); // Coordinate of current node

    double x = node_crd[0], y = node_crd[1];

    if (icase == 1) {
      B_n[0][n] = std::max(
        0.0, 0.25 - 5 * ((x - 0.5) * (x - 0.5) + (y - 0.5) * (y - 0.5))); // non-smooth bathymetry
    } else if (icase == 2) {
      B_n[0][n] = 0.25 - 0.25 * x * std::sin(pi * y); // non-zero bathymetry at boundary
    }
  }

  S->Get<CompositeVector>("surface-bathymetry").ScatterMasterToGhosted("node");

  for (int c = 0; c < ncells_wghost; ++c) {
    const Amanzi::AmanziGeometry::Point& xc = mesh->getCellCentroid(c);

    auto cfaces = mesh->getCellFaces(c);
    auto cnodes = mesh->getCellNodes(c);
    auto cedges = mesh->getCellEdges(c);

    int nfaces_cell = cfaces.size();

    B_c[0][c] = 0;

    // Compute cell averaged bathymetrt (Bc)
    for (int f = 0; f < nfaces_cell; ++f) {
      Amanzi::AmanziGeometry::Point x0, x1;
      int edge = cfaces[f];

      auto face_nodes = mesh->getFaceNodes(edge);
      int n0 = face_nodes[0], n1 = face_nodes[1];

      x0 = mesh->getNodeCoordinate(n0);
      x1 = mesh->getNodeCoordinate(n1);

      double area = norm((xc - x0) ^ (xc - x1)) / 2.0;

      B_c[0][c] += (area / mesh->getCellVolume(c)) * (B_n[0][n0] + B_n[0][n1]) / 2.0;
    }

    // Perturb the solution; change time period t_new to at least 10.0
    //    if ((xc[0] - 0.3)*(xc[0] - 0.3) + (xc[1] - 0.3)*(xc[1] - 0.3) < 0.1 * 0.1) {
    //      ht_c[0][c] = H_inf + 0.01;
    //     } else {
    //      ht_c[0][c] = H_inf;
    //    }

    ht_c[0][c] = H_inf;
    h_c[0][c] = ht_c[0][c] - B_c[0][c];
    vel_c[0][c] = 0.0;
    vel_c[1][c] = 0.0;
    q_c[0][c] = h_c[0][c] * vel_c[0][c];
    q_c[1][c] = h_c[0][c] * vel_c[1][c];
    p_c[0][c] = patm + rho * g * h_c[0][c];
  }

  S->Get<CompositeVector>("surface-bathymetry").ScatterMasterToGhosted("cell");
  S->Get<CompositeVector>("surface-total_depth").ScatterMasterToGhosted("cell");
  S->Get<CompositeVector>("surface-ponded_depth").ScatterMasterToGhosted("cell");
  S->Get<CompositeVector>("surface-velocity").ScatterMasterToGhosted("cell");
  S->Get<CompositeVector>("surface-discharge").ScatterMasterToGhosted("cell");
}


// ---- ---- ---- ---- ---- ---- ---- ----
// Error
// ---- ---- ---- ---- ---- ---- ---- ----
void
error(Teuchos::RCP<const Amanzi::AmanziMesh::Mesh> mesh,
      Epetra_MultiVector& ht_ex,
      Epetra_MultiVector& vel_ex,
      const Epetra_MultiVector& ht,
      const Epetra_MultiVector& vel,
      double& err_max,
      double& err_L1,
      double& hmax)
{
  int ncells_owned =
    mesh->getNumEntities(Amanzi::AmanziMesh::Entity_kind::CELL, Amanzi::AmanziMesh::Parallel_kind::OWNED);

  err_max = 0.0;
  err_L1 = 0.0;
  hmax = 0.0;

  double err_max_u = 0.0, err_max_v = 0.0, err_L1_u = 0.0, err_L1_v = 0.0, tmp_hmax = 0.0;

  for (int c = 0; c < ncells_owned; ++c) {
    double tmp = std::abs(ht_ex[0][c] - ht[0][c]);
    err_max = std::max(err_max, tmp);
    err_L1 += tmp * mesh->getCellVolume(c);
    tmp_hmax = std::sqrt(mesh->getCellVolume(c));

    tmp = std::abs(vel_ex[0][c] - vel[0][c]);
    err_max_u = std::max(err_max_u, tmp);
    err_L1_u += tmp * mesh->getCellVolume(c);

    tmp = std::abs(vel_ex[1][c] - vel[1][c]);
    err_max_v = std::max(err_max_v, tmp);
    err_L1_v += tmp * mesh->getCellVolume(c);

    hmax = std::max(tmp_hmax, hmax);
  }

  double err_max_tmp, err_L1_tmp;

  mesh->getComm()->MaxAll(&err_max, &err_max_tmp, 1);
  mesh->getComm()->SumAll(&err_L1, &err_L1_tmp, 1);

  err_max = err_max_tmp;
  err_L1 = err_L1_tmp;

  mesh->getComm()->MaxAll(&err_max_u, &err_max_tmp, 1);
  mesh->getComm()->SumAll(&err_L1_u, &err_L1_tmp, 1);

  err_max_u = err_max_tmp;
  err_L1_u = err_L1_tmp;

  mesh->getComm()->MaxAll(&err_max_v, &err_max_tmp, 1);
  mesh->getComm()->SumAll(&err_L1_v, &err_L1_tmp, 1);

  err_max_v = err_max_tmp;
  err_L1_v = err_L1_tmp;

  std::cout.precision(6);

  std::cout << std::scientific << "H err_max: " << err_max << std::endl;
  std::cout << std::scientific << "H err_L1: " << err_L1 << std::endl;

  std::cout << std::scientific << "u err_max: " << err_max_u << std::endl;
  std::cout << std::scientific << "u err_L1: " << err_L1_u << std::endl;

  std::cout << std::scientific << "v err_max: " << err_max_v << std::endl;
  std::cout << std::scientific << "v err_L1: " << err_L1_v << std::endl;
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

  if (MyPID == 0) { std::cout << "Test: 2D Shallow water: Lake at rest" << std::endl; }

  // Read parameter list
  std::string xmlFilename;
  if (icase == 1) {
    xmlFilename = "test/shallow_water_lake_at_rest_case1.xml";
  } else if (icase == 2) {
    xmlFilename = "test/shallow_water_lake_at_rest_case2.xml";
  }
  Teuchos::RCP<Teuchos::ParameterList> plist = Teuchos::getParametersFromXmlFile(xmlFilename);

  // Create a mesh framework
  ParameterList regions_list = plist->get<Teuchos::ParameterList>("regions");
  auto gm = Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(2, regions_list, *comm));

  // Creat a mesh
  auto fac_list = Teuchos::rcp(new Teuchos::ParameterList()); 
  fac_list->set<bool>("request edges", true); 
  MeshFactory meshfactory(comm, gm, fac_list);
  meshfactory.set_preference(Preference({ Framework::MSTK }));

  RCP<Mesh> mesh;
  if (icase == 1) {
    // Rectangular mesh
    mesh = meshfactory.create(0.0, 0.0, 1.0, 1.0, 25, 25);
  } else if (icase == 2) {
    // Triangular mesh
    mesh = meshfactory.create("test/triangular16.exo");
  }
  // Other polygonal meshes
  //  RCP<Mesh> mesh = meshfactory.create ("test/median15x16.exo");
  //  RCP<Mesh> mesh = meshfactory.create ("test/random40.exo");

  // Create a state

  Teuchos::ParameterList state_list = plist->sublist("state");
  RCP<State> S = rcp(new State(state_list));
  S->RegisterMesh("surface", mesh);

  Teuchos::RCP<TreeVector> soln = Teuchos::rcp(new TreeVector());

  Teuchos::ParameterList sw_list = plist->sublist("PKs").sublist("shallow water");

  // Create a shallow water PK
  ShallowWater_PK SWPK(sw_list, plist, S, soln);
  SWPK.Setup();
  S->Setup();
  S->InitializeFields();
  S->InitializeEvaluators();
  S->set_time(0.0);
  SWPK.Initialize();

  lake_at_rest_setIC(mesh, S, icase);
  S->CheckAllFieldsInitialized();

  const auto& ht = *S->Get<CompositeVector>("surface-total_depth").ViewComponent("cell");
  const auto& vel = *S->Get<CompositeVector>("surface-velocity").ViewComponent("cell");

  // Create screen io
  auto vo = Teuchos::rcp(new Amanzi::VerboseObject("ShallowWater", *plist));
  WriteStateStatistics(*S, *vo);

  // Advance in time
  double t_old(0.0), t_new(0.0), dt;

  // Initialize io
  Teuchos::ParameterList iolist;
  std::string fname;
  fname = "SW_sol";
  iolist.get<std::string>("file name base", fname);
  OutputXDMF io(iolist, mesh, true, false);

  std::string passwd("state");

  int iter = 0;
  std::vector<double> dx, Linferror, L1error, L2error;

  while (t_new < 1.0) {
    double t_out = t_new;

    if (iter % 30 == 0) {
      Epetra_MultiVector ht_ex(ht);
      Epetra_MultiVector vel_ex(vel);
      lake_at_rest_exact_field(mesh, ht_ex, vel_ex, t_out);

      IO_Fields(t_out, iter, MyPID, io, *S, nullptr, nullptr);
    }

    dt = SWPK.get_dt();

    t_new = t_old + dt;

    SWPK.AdvanceStep(t_old, t_new);
    SWPK.CommitStep(t_old, t_new, Tags::DEFAULT);

    t_old = t_new;
    iter += 1;
  } // time loop

  if (MyPID == 0) { std::cout << "Time-stepping finished. " << std::endl; }

  double t_out = t_new;

  double err_max, err_L1, hmax;
  Epetra_MultiVector ht_ex(ht);
  Epetra_MultiVector vel_ex(vel);

  lake_at_rest_exact_field(mesh, ht_ex, vel_ex, t_out);
  error(mesh, ht_ex, vel_ex, ht, vel, err_max, err_L1, hmax);

  std::cout << "hmax: " << hmax << std::endl;

  dx.push_back(hmax);
  Linferror.push_back(err_max);
  L1error.push_back(err_L1);

  CHECK_CLOSE(0.0, L1error[0], 1e-12);
  CHECK_CLOSE(0.0, Linferror[0], 1e-12);

  IO_Fields(t_out, iter, MyPID, io, *S, nullptr, nullptr);

  std::cout << "Computed error H (L_1): " << L1error[0] << std::endl;
  std::cout << "Computed error H (L_inf): " << Linferror[0] << std::endl;
}

TEST(SHALLOW_WATER_LAKE_AT_REST)
{
  RunTest(1); // rectangular mesh
  RunTest(2); // triangular mesh
}
