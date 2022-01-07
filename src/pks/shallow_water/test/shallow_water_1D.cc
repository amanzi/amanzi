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
#include "Mesh.hh"
#include "MeshFactory.hh"
#include "OutputXDMF.hh"

#include "ShallowWater_PK.hh"


using namespace Amanzi;

//--------------------------------------------------------------
// analytic solution
//--------------------------------------------------------------
void dam_break_1D_setIC(Teuchos::RCP<const Amanzi::AmanziMesh::Mesh> mesh,
                        Teuchos::RCP<Amanzi::State>& S)
{
  int ncells_owned = mesh->num_entities(Amanzi::AmanziMesh::CELL, Amanzi::AmanziMesh::Parallel_type::OWNED);

  std::string passwd = "state";

  auto& B_vec_c = *S->GetW<CompositeVector>("surface-bathymetry", Tags::DEFAULT, passwd).ViewComponent("cell");
  auto& h_vec_c = *S->GetW<CompositeVector>("surface-ponded_depth", Tags::DEFAULT, passwd).ViewComponent("cell");
  auto& ht_vec_c = *S->GetW<CompositeVector>("surface-total_depth", Tags::DEFAULT, passwd).ViewComponent("cell");

  for (int c = 0; c < ncells_owned; c++) {
    Amanzi::AmanziGeometry::Point xc = mesh->cell_centroid(c);
    if (xc[0] < 1000.) {
      h_vec_c[0][c] = 10.;
    } else {
      h_vec_c[0][c] = 0.; //0.0000000001;
    }
    ht_vec_c[0][c] = h_vec_c[0][c] + B_vec_c[0][c];
  }

  S->GetW<CompositeVector>("surface-velocity", Tags::DEFAULT, passwd).PutScalar(0.0);
  S->GetW<CompositeVector>("surface-discharge", Tags::DEFAULT, "surface-discharge").PutScalar(0.0);
}


void dam_break_1D_exact(double hL, double x0, double t, double x, double &h, double &u)
{
  double g = 9.81;
  double xA = x0 - t*std::sqrt(g*hL);
  double xB = x0 + 2.*t*std::sqrt(g*hL);

  if (0. <= x && x < xA) {
    h = hL;
    u = 0.;
  } else {
    if (xA <= x && x < xB) {
      h = 4./(9.*g)*( std::sqrt(g*hL) - (x-x0)/(2.*t) )*( std::sqrt(g*hL) - (x-x0)/(2.*t) );
      u = 2./3.*( (x-x0)/t + std::sqrt(g*hL) );
    } else {
      h = 0.;
      u = 0.;
    }
  }
}


void dam_break_1D_exact_field(Teuchos::RCP<const Amanzi::AmanziMesh::Mesh> mesh,
                              Epetra_MultiVector& hh_ex, Epetra_MultiVector& vx_ex, double t)
{
  double hL, x0, x, h, u;

  hL = 10.;
  x0 = 1000.;

  int ncells_owned = mesh->num_entities(Amanzi::AmanziMesh::CELL, Amanzi::AmanziMesh::Parallel_type::OWNED);

  for (int c = 0; c < ncells_owned; c++) {
    Amanzi::AmanziGeometry::Point xc = mesh->cell_centroid(c);

    x = xc[0];

    dam_break_1D_exact(hL, x0, t, x, h, u);
    hh_ex[0][c] = h;
    vx_ex[0][c] = u;
  }
}


void error(Teuchos::RCP<const Amanzi::AmanziMesh::Mesh> mesh,
           Epetra_MultiVector& hh_ex, Epetra_MultiVector& vx_ex,
           const Epetra_MultiVector& hh, const Epetra_MultiVector& vx,
           double& err_max, double& err_L1, double& hmax)
{
  int ncells_owned = mesh->num_entities(Amanzi::AmanziMesh::CELL, Amanzi::AmanziMesh::Parallel_type::OWNED);

  err_max = 0.;
  err_L1 = 0.;
  hmax = 0.;

  for (int c = 0; c < ncells_owned; c++) {
    double vol = mesh->cell_volume(c);
    double tmp = std::abs(hh_ex[0][c] - hh[0][c]);
    err_max = std::max(err_max, tmp);
    err_L1 += tmp * vol;
    hmax = std::sqrt(vol);
  }

  double err_max_tmp(err_max);
  double err_L1_tmp(err_L1);

  mesh->get_comm()->MaxAll(&err_max_tmp, &err_max, 1);
  mesh->get_comm()->SumAll(&err_L1_tmp, &err_L1, 1);

  if (mesh->get_comm()->MyPID() == 0) {
    std::cout << "err_max = " << err_max << std::endl;
    std::cout << "err_L1  = " << err_L1 << std::endl;
  }
}


void ConservationCheck(Teuchos::RCP<const Amanzi::AmanziMesh::Mesh> mesh, const Epetra_MultiVector& hh,
                       const Epetra_MultiVector& ht,
                       const Epetra_MultiVector& vel,
                       const Epetra_MultiVector& B, double& TE)
{
  int ncells_owned = mesh->num_entities(Amanzi::AmanziMesh::CELL, Amanzi::AmanziMesh::Parallel_type::OWNED);

  double g = 9.81;
  double KE = 0., IE = 0.;

  for (int c = 0; c < ncells_owned; c++) {
    double volume = mesh->cell_volume(c);
    KE += hh[0][c] * (vel[0][c] * vel[0][c] + vel[1][c] * vel[1][c]) * volume;
    IE += g * hh[0][c] * hh[0][c] * volume;
  }

  double TE_tmp = (KE + IE) / 2;
  mesh->get_comm()->SumAll(&TE_tmp, &TE, 1);
}


/* **************************************************************** */
TEST(SHALLOW_WATER_1D) {
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

  // create a mesh framework
  ParameterList regions_list = plist->get<Teuchos::ParameterList>("regions");
  auto gm = Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(2, regions_list, *comm));

  // create a mesh
  bool request_faces = true, request_edges = false;
  MeshFactory meshfactory(comm,gm);
  meshfactory.set_preference(Preference({Framework::MSTK}));

  RCP<Mesh> mesh = meshfactory.create(0.0, 0.0, 2000.0, 50.0, 1000, 1, request_faces, request_edges);
  // mesh = meshfactory.create("test/median63x64.exo",request_faces,request_edges); // works only with first order, no reconstruction

  // create a state
  Teuchos::ParameterList state_list = plist->sublist("state");
  RCP<State> S = rcp(new State(state_list));
  S->RegisterMesh("surface", mesh);

  Teuchos::RCP<TreeVector> soln = Teuchos::rcp(new TreeVector());

  Teuchos::ParameterList sw_list = plist->sublist("PKs").sublist("shallow water");

  // create a shallow water PK
  ShallowWater_PK SWPK(sw_list,plist,S,soln);
  SWPK.Setup(S.ptr());
  S->Setup();
  S->InitializeFields();
  S->InitializeEvaluators();
  S->set_time(0.0);
  SWPK.Initialize(S.ptr());

  dam_break_1D_setIC(mesh, S);
  S->CheckAllFieldsInitialized();

  const auto& hh = *S->Get<CompositeVector>("surface-ponded_depth").ViewComponent("cell");
  const auto& ht = *S->Get<CompositeVector>("surface-total_depth").ViewComponent("cell");
  const auto& vel = *S->Get<CompositeVector>("surface-velocity").ViewComponent("cell");
  const auto& q = *S->Get<CompositeVector>("surface-discharge").ViewComponent("cell");
  const auto& B = *S->Get<CompositeVector>("surface-bathymetry").ViewComponent("cell");

  // create pid vector
  Epetra_MultiVector pid(B);
  for (int c = 0; c < pid.MyLength(); c++) pid[0][c] = MyPID;

  // create screen io
  auto vo = Teuchos::rcp(new Amanzi::VerboseObject("ShallowWater", sw_list));
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
  double TEini, TEfin;

  while (t_new < 10.) {
    // cycle 1, time t
    double t_out = t_new;

    Epetra_MultiVector hh_ex(hh);
    Epetra_MultiVector vel_ex(vel);

    dam_break_1D_exact_field(mesh, hh_ex, vel_ex, t_out);

    if (iter % 25 == 0) {
      io.InitializeCycle(t_out, iter, "");
      io.WriteVector(*hh(0), "depth", AmanziMesh::CELL);
      io.WriteVector(*ht(0), "total_depth", AmanziMesh::CELL);
      io.WriteVector(*vel(0), "vx", AmanziMesh::CELL);
      io.WriteVector(*vel(1), "vy", AmanziMesh::CELL);
      io.WriteVector(*q(0), "qx", AmanziMesh::CELL);
      io.WriteVector(*q(1), "qy", AmanziMesh::CELL);
      io.WriteVector(*B(0), "B", AmanziMesh::CELL);
      io.WriteVector(*pid(0), "pid", AmanziMesh::CELL);

      io.WriteVector(*hh_ex(0), "hh_ex", AmanziMesh::CELL);
      io.WriteVector(*vel_ex(0), "vx_ex", AmanziMesh::CELL);
      io.FinalizeCycle();
    }

    dt = SWPK.get_dt();

    t_new = t_old + dt;

    SWPK.AdvanceStep(t_old, t_new);
    SWPK.CommitStep(t_old, t_new, S);

    t_old = t_new;
    iter++;

    if (iter == 1) {
      ConservationCheck(mesh, hh, ht, vel, B, TEini);
      if (MyPID == 0) std::cout << "cycle= " << iter << "  initial TE=" << TEini << "  dt=" << dt << std::endl;
    }
  }

  ConservationCheck(mesh, hh, ht, vel, B, TEfin);
  if (MyPID == 0) std::cout << "cycle= " << iter << "  final TE=" << TEfin << "  dt=" << dt << std::endl;
  CHECK_CLOSE(TEini, TEfin, 6e-3 * TEini);

  // error calculation at the time time
  double t_out = t_new;
  double err_max, err_L1, hmax;
  Epetra_MultiVector hh_ex(hh), vel_ex(vel);

  dam_break_1D_exact_field(mesh, hh_ex, vel_ex, t_out);

  error(mesh, hh_ex, vel_ex, hh, vel, err_max, err_L1, hmax);

  CHECK_CLOSE(1./hmax, err_max, 0.5);

  // save final state values into HDF5 file
  io.InitializeCycle(t_out, iter, "");
  io.WriteVector(*hh(0), "depth", AmanziMesh::CELL);
  io.WriteVector(*ht(0), "total_depth", AmanziMesh::CELL);
  io.WriteVector(*vel(0), "vx", AmanziMesh::CELL);
  io.WriteVector(*vel(1), "vy", AmanziMesh::CELL);
  io.WriteVector(*q(0), "qx", AmanziMesh::CELL);
  io.WriteVector(*q(1), "qy", AmanziMesh::CELL);
  io.WriteVector(*B(0), "B", AmanziMesh::CELL);
  io.WriteVector(*pid(0), "pid", AmanziMesh::CELL);

  io.WriteVector(*hh_ex(0), "hh_ex", AmanziMesh::CELL);
  io.WriteVector(*vel_ex(0), "vx_ex", AmanziMesh::CELL);
  io.FinalizeCycle();

  /*
  std::ofstream myfile;
  myfile.open("solution_" + std::to_string(MyPID) + ".txt");
  int ncells_owned = mesh->num_entities(Amanzi::AmanziMesh::CELL, Amanzi::AmanziMesh::Parallel_type::OWNED);
  for (int c = 0; c < ncells_owned; c++) {
      AmanziGeometry::Point xc = mesh->cell_centroid(c);
      myfile << xc[0] << " " << hh[0][c] << "\n";
  }
  myfile.close();
  */
  WriteStateStatistics(*S, *vo);
}
