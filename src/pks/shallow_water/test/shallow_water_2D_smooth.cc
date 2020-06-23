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
// bottom topography
//--------------------------------------------------------------
double Bathymetry(double x, double y)
{
  return 0.;
}

double hf(double x)
{
  return 2.*cos(x) + 2.*x*sin(x) + 1./8.*cos(2.*x) + x/4.*sin(2.*x) + 12./16.*x*x;
}

void vortex_2D_exact(double t, double x, double y, double &h, double &u, double &v)
{
  double g = 9.81;
  double gmm = 15.;
  double omega = M_PI;
  double u_inf = 6., v_inf = 0.;
  double H_inf = 10.;
  double xc = 5., yc = 5.;
  double H0;

  double xt = x - xc - u_inf*t;
  double yt = y - yc - v_inf*t;

  double rc = std::sqrt(xt*xt + yt*yt); // distance from the vortex core

  // if (omega*rc <= M_PI) {
  if (rc <= 1.) {
    h = H_inf + 1./g*(gmm/omega)*(gmm/omega)*(hf(omega*rc)-hf(M_PI));
    u = u_inf + gmm*(1.+cos(omega*rc))*(yc-y);
    v = v_inf + gmm*(1.+cos(omega*rc))*(x-xc);
  } else {
    h = H_inf;
    u = u_inf;
    v = v_inf;
  }
}

void vortex_2D_exact_field(Teuchos::RCP<const Amanzi::AmanziMesh::Mesh> mesh, Epetra_MultiVector& hh_ex, Epetra_MultiVector& vx_ex, Epetra_MultiVector& vy_ex, double t)
{
  double x, y, h, u, v;

  int ncells_owned = mesh->num_entities(Amanzi::AmanziMesh::CELL, Amanzi::AmanziMesh::Parallel_type::OWNED);

  for (int c = 0; c < ncells_owned; c++) {

    Amanzi::AmanziGeometry::Point xc = mesh->cell_centroid(c);

    x = xc[0];
    y = xc[1];

    vortex_2D_exact(t, x, y, h, u, v);
    hh_ex[0][c] = h;
    vx_ex[0][c] = u;
    vy_ex[0][c] = v;

  }
}

void vortex_2D_setIC(Teuchos::RCP<const Amanzi::AmanziMesh::Mesh> mesh_, Teuchos::RCP<Amanzi::State>& S_)
{
  int ncells_owned = mesh_->num_entities(Amanzi::AmanziMesh::CELL, Amanzi::AmanziMesh::Parallel_type::OWNED);

  std::string passwd_ = "state";

  Epetra_MultiVector& B_vec_c = *S_->GetFieldData("surface-bathymetry",passwd_)->ViewComponent("cell");

  for (int c = 0; c < ncells_owned; c++) {
    Amanzi::AmanziGeometry::Point xc = mesh_->cell_centroid(c);
    B_vec_c[0][c] = Bathymetry(xc[0],xc[1]);
  }

  Epetra_MultiVector& h_vec_c = *S_->GetFieldData("surface-ponded_depth",passwd_)->ViewComponent("cell");
  Epetra_MultiVector& ht_vec_c = *S_->GetFieldData("surface-total_depth",passwd_)->ViewComponent("cell");
  Epetra_MultiVector& vx_vec_c = *S_->GetFieldData("surface-velocity-x",passwd_)->ViewComponent("cell");
  Epetra_MultiVector& vy_vec_c = *S_->GetFieldData("surface-velocity-y",passwd_)->ViewComponent("cell");
  Epetra_MultiVector& qx_vec_c = *S_->GetFieldData("surface-discharge-x",passwd_)->ViewComponent("cell");
  Epetra_MultiVector& qy_vec_c = *S_->GetFieldData("surface-discharge-y",passwd_)->ViewComponent("cell");

  for (int c = 0; c < ncells_owned; c++) {
    Amanzi::AmanziGeometry::Point xc = mesh_->cell_centroid(c);
    double h, u, v;
    vortex_2D_exact(0., xc[0], xc[1], h, u, v);
    h_vec_c[0][c] = h;
    ht_vec_c[0][c] = h_vec_c[0][c] + B_vec_c[0][c];
    vx_vec_c[0][c] = u;
    vy_vec_c[0][c] = v;
    qx_vec_c[0][c] = u*h;
    qy_vec_c[0][c] = v*h;
  }
}

void error(Teuchos::RCP<const Amanzi::AmanziMesh::Mesh> mesh,
           Epetra_MultiVector& hh_ex, Epetra_MultiVector& vx_ex, Epetra_MultiVector& vy_ex, 
           const Epetra_MultiVector& hh, const Epetra_MultiVector& vx, const Epetra_MultiVector& vy,
           double& err_max, double& err_L1, double& hmax)
{
  int ncells_owned = mesh->num_entities(Amanzi::AmanziMesh::CELL, Amanzi::AmanziMesh::Parallel_type::OWNED);

  err_max = 0.;
  err_L1 = 0.;
  hmax = 0.;

  for (int c = 0; c < ncells_owned; c++) {
    err_max = std::max(err_max,std::abs(hh_ex[0][c]-hh[0][c]));
    err_L1 += std::abs(hh_ex[0][c]-hh[0][c])*mesh->cell_volume(c);
    hmax = std::sqrt(mesh->cell_volume(c));
  }

  std::cout << "err_max = " << err_max << std::endl;
  std::cout << "err_L1  = " << err_L1 << std::endl;
}


/* **************************************************************** */
TEST(SHALLOW_WATER_2D_SMOOTH) {
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
  MeshFactory meshfactory(comm,gm);
  meshfactory.set_preference(Preference({Framework::MSTK, Framework::STK}));
  if (MyPID == 0) std::cout << "Mesh factory created." << std::endl;

  double order = 1.;

  std::vector<double> dx, Linferror, L1error, L2error;

  for (int NN = 20; NN <= 160; NN *= 2) {

    RCP<const Mesh> mesh;
    mesh = meshfactory.create(0.0, 0.0, 10.0, 10.0, NN, NN, request_faces, request_edges);
    // mesh = meshfactory.create("test/median63x64.exo",request_faces,request_edges); 
    // works only with first order, no reconstruction
    if (MyPID == 0) std::cout << "Mesh created." << std::endl;

    // create a state
    RCP<State> S = rcp(new State());
    S->RegisterMesh("surface",rcp_const_cast<Mesh>(mesh));
    S->set_time(0.0);
    if (MyPID == 0) std::cout << "State created." << std::endl;

    Teuchos::RCP<TreeVector> soln = Teuchos::rcp(new TreeVector());

    Teuchos::ParameterList pk_tree = plist->sublist("PK tree").sublist("shallow water");

    // create a shallow water PK
    ShallowWater_PK SWPK(pk_tree,plist,S,soln);
    SWPK.Setup(S.ptr());
    S->Setup();
    S->InitializeFields();
    S->InitializeEvaluators();
    SWPK.Initialize(S.ptr());
    vortex_2D_setIC(mesh,S);
    if (MyPID == 0) std::cout << "Shallow water PK created." << std::endl;

    // create screen io
    auto vo = Teuchos::rcp(new Amanzi::VerboseObject("ShallowWater", *plist));
    WriteStateStatistics(S.ptr(), vo);

    // advance in time
    double t_old(0.0), t_new(0.0), dt;

    // initialize io
    Teuchos::ParameterList iolist;
    std::string fname;
    fname = "SW_sol";
    iolist.get<std::string>("file name base", fname);
    OutputXDMF io(iolist, mesh, true, false);

    std::string passwd("state");

    auto& h_vec = *S->GetFieldData("surface-ponded_depth",passwd)->ViewComponent("cell");
    auto& u_vec = *S->GetFieldData("surface-velocity-x",passwd)->ViewComponent("cell");
    auto& v_vec = *S->GetFieldData("surface-velocity-y",passwd)->ViewComponent("cell");

    int iter = 0;
    bool flag = true;

    while (t_new < 0.001) {
      // cycle 1, time t
      double t_out = t_new;

      const Epetra_MultiVector& hh = *S->GetFieldData("surface-ponded_depth",passwd)->ViewComponent("cell");
      const Epetra_MultiVector& ht = *S->GetFieldData("surface-total_depth",passwd)->ViewComponent("cell");
      const Epetra_MultiVector& vx = *S->GetFieldData("surface-velocity-x",passwd)->ViewComponent("cell");
      const Epetra_MultiVector& vy = *S->GetFieldData("surface-velocity-y",passwd)->ViewComponent("cell");
      const Epetra_MultiVector& qx = *S->GetFieldData("surface-discharge-x",passwd)->ViewComponent("cell");
      const Epetra_MultiVector& qy = *S->GetFieldData("surface-discharge-y",passwd)->ViewComponent("cell");
      const Epetra_MultiVector& B  = *S->GetFieldData("surface-bathymetry",passwd)->ViewComponent("cell");
      const Epetra_MultiVector& pid = *S->GetFieldData("surface-PID",passwd)->ViewComponent("cell");

      Epetra_MultiVector hh_ex(hh);
      Epetra_MultiVector vx_ex(vx);
      Epetra_MultiVector vy_ex(vy);

      vortex_2D_exact_field(mesh, hh_ex, vx_ex, vy_ex, t_out);

      io.InitializeCycle(t_out, iter);
      io.WriteVector(*hh(0), "depth", AmanziMesh::CELL);
      io.WriteVector(*ht(0), "total_depth", AmanziMesh::CELL);
      io.WriteVector(*vx(0), "vx", AmanziMesh::CELL);
      io.WriteVector(*vy(0), "vy", AmanziMesh::CELL);
      io.WriteVector(*qx(0), "qx", AmanziMesh::CELL);
      io.WriteVector(*qy(0), "qy", AmanziMesh::CELL);
      io.WriteVector(*B(0), "B", AmanziMesh::CELL);
      io.WriteVector(*pid(0), "pid", AmanziMesh::CELL);
      io.WriteVector(*hh_ex(0), "hh_ex", AmanziMesh::CELL);
      io.WriteVector(*vx_ex(0), "vx_ex", AmanziMesh::CELL);
      io.WriteVector(*vy_ex(0), "vy_ex", AmanziMesh::CELL);
      io.FinalizeCycle();

      dt = SWPK.get_dt();

      if (iter < 10) dt = 0.01*dt;

      t_new = t_old + dt;

      SWPK.AdvanceStep(t_old, t_new);
      SWPK.CommitStep(t_old, t_new, S);

      t_old = t_new;
      iter++;
    }

    if (MyPID == 0) std::cout << "Time-stepping finished." << std::endl;

    // cycle 1, time t
    double t_out = t_new;

    const Epetra_MultiVector& hh = *S->GetFieldData("surface-ponded_depth",passwd)->ViewComponent("cell");
    const Epetra_MultiVector& ht = *S->GetFieldData("surface-total_depth",passwd)->ViewComponent("cell");
    const Epetra_MultiVector& vx = *S->GetFieldData("surface-velocity-x",passwd)->ViewComponent("cell");
    const Epetra_MultiVector& vy = *S->GetFieldData("surface-velocity-y",passwd)->ViewComponent("cell");
    const Epetra_MultiVector& qx = *S->GetFieldData("surface-discharge-x",passwd)->ViewComponent("cell");
    const Epetra_MultiVector& qy = *S->GetFieldData("surface-discharge-y",passwd)->ViewComponent("cell");
    const Epetra_MultiVector& B  = *S->GetFieldData("surface-bathymetry",passwd)->ViewComponent("cell");
    const Epetra_MultiVector& pid = *S->GetFieldData("surface-PID",passwd)->ViewComponent("cell");

    Epetra_MultiVector hh_ex(hh);
    Epetra_MultiVector vx_ex(vx);
    Epetra_MultiVector vy_ex(vy);

    vortex_2D_exact_field(mesh, hh_ex, vx_ex, vy_ex, t_out);

    double err_max, err_L1, hmax;

    error(mesh, hh_ex, vx_ex, vy_ex, hh, vx, vy, err_max, err_L1, hmax);

    dx.push_back(hmax);
    Linferror.push_back(err_max);
    L1error.push_back(err_L1);

    io.InitializeCycle(t_out, iter);
    io.WriteVector(*hh(0), "depth", AmanziMesh::CELL);
    io.WriteVector(*ht(0), "total_depth", AmanziMesh::CELL);
    io.WriteVector(*vx(0), "vx", AmanziMesh::CELL);
    io.WriteVector(*vy(0), "vy", AmanziMesh::CELL);
    io.WriteVector(*qx(0), "qx", AmanziMesh::CELL);
    io.WriteVector(*qy(0), "qy", AmanziMesh::CELL);
    io.WriteVector(*B(0), "B", AmanziMesh::CELL);
    io.WriteVector(*pid(0), "pid", AmanziMesh::CELL);
    io.WriteVector(*hh_ex(0), "hh_ex", AmanziMesh::CELL);
    io.WriteVector(*vx_ex(0), "vx_ex", AmanziMesh::CELL);
    io.WriteVector(*vy_ex(0), "vy_ex", AmanziMesh::CELL);
    io.FinalizeCycle();
  } // NN

  order = Amanzi::Utils::bestLSfit(dx, L1error);

  std::cout << "computed order = " << order << std::endl;

  CHECK_CLOSE(1.5, order, 0.2);
}
