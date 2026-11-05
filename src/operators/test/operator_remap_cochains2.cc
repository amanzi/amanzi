/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Operators
*/

#include <array>
#include <cstdlib>
#include <cmath>
#include <functional>
#include <iostream>
#include <vector>

// TPLs
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "UnitTest++.h"

// Amanzi
#include "Mesh.hh"
#include "MeshFactory.hh"
#include "RemapCochains.hh"
#include "OutputXDMF.hh"

// Amanzi::Operators
#include "MeshDeformation.hh"

using namespace Amanzi;

static AmanziGeometry::Point fun0_v(const AmanziGeometry::Point& p) {
 return Amanzi::AmanziGeometry::Point(1.0, 2.0, 3.0);
}
static AmanziGeometry::Point fun0_u(const AmanziGeometry::Point& p) {
 double x = p[0];
 double y = p[1];
 double z = p[2];
 return 0.5 * AmanziGeometry::Point(2 * z - 3 * y, 3 * x - z, y - 2 * x);
}

static AmanziGeometry::Point fun1_v(const AmanziGeometry::Point& p) {
 return AmanziGeometry::Point(p[0], 2 * p[1], -3 * p[2]);
}
static AmanziGeometry::Point fun1_u(const AmanziGeometry::Point& p) {
 double x = p[0];
 double y = p[1];
 double z = p[2];
 return AmanziGeometry::Point(0.0, -3 * x * z, -2 * x * y);
}

static AmanziGeometry::Point fun2_v(const AmanziGeometry::Point& p) {
 double x = p[0];
 double y = p[1];
 double z = p[2];
 return AmanziGeometry::Point(-y * std::cos(y * z), -z * std::cos(z * x), -x * std::cos(x * y));
}
static AmanziGeometry::Point fun2_u(const AmanziGeometry::Point& p) {
 double x = p[0];
 double y = p[1];
 double z = p[2];
 return AmanziGeometry::Point(std::sin(x * y), std::sin(y * z), std::sin(z * x));
}

static AmanziGeometry::Point tornado_v(const AmanziGeometry::Point& p) {
  constexpr double W = 1.0;
  constexpr double C = 0.06;
  constexpr double a = 0.04;
  constexpr double b = 0.08;

  const double X = p[0] - 0.5;
  const double Y = p[1] - 0.5;
  const double z = p[2];
  const double s = X * X + Y * Y;

  const double g = std::sin(M_PI * z);
  const double gp = M_PI * std::cos(M_PI * z);
  const double ea = std::exp(-s / a);
  const double eb = std::exp(-s / b);

  // v = curl(u)
  const double vx = -2.0 * C * g * Y * eb / b - 0.5 * W * gp * X * ea;
  const double vy =  2.0 * C * g * X * eb / b - 0.5 * W * gp * Y * ea;
  const double vz = W * g * ea * (1.0 - s / a);

  return AmanziGeometry::Point(vx, vy, vz);
}
static AmanziGeometry::Point tornado_u(const AmanziGeometry::Point& p) {
  constexpr double W = 1.0;    // axial/updraft strength
  constexpr double C = 0.06;   // swirl-potential strength
  constexpr double a = 0.04;   // updraft-core radius squared
  constexpr double b = 0.08;   // swirl-core radius squared

  const double X = p[0] - 0.5;
  const double Y = p[1] - 0.5;
  const double z = p[2];
  const double s = X * X + Y * Y;
  const double g = std::sin(M_PI * z);

  const double F = W * g * std::exp(-s / a);
  const double psi = C * g * std::exp(-s / b);

  return AmanziGeometry::Point(-0.5 * Y * F, 0.5 * X * F, psi);
}


/* *****************************************************************
* Simplified mesh deformation
***************************************************************** */
static void 
Deform3D(const Teuchos::RCP<AmanziMesh::Mesh>& mesh1,
         const Teuchos::RCP<AmanziMesh::Mesh>& mesh0,
         double t,
         std::function<AmanziGeometry::Point(AmanziGeometry::Point, double)> fun)
{
  CompositeVectorSpace cvs;
  cvs.SetMesh(mesh1)->SetGhosted(true)->AddComponent("node", AmanziMesh::Entity_kind::NODE, 3);

  int nnodes = mesh1->getNumEntities(AmanziMesh::Entity_kind::NODE, AmanziMesh::Parallel_kind::ALL);
  AmanziMesh::Entity_ID_View nodeids("nodeids", nnodes);
  AmanziMesh::Point_View new_positions("newpos", nnodes);

  for (int v = 0; v < nnodes; ++v) {
    const auto& xv = mesh0->getNodeCoordinate(v);
    nodeids[v] = v;
    new_positions[v] = fun(xv, t);
  }
  AmanziMesh::deform(*mesh1, nodeids, new_positions);
}


static AmanziGeometry::Point deform1(const AmanziGeometry::Point& xv, double t)
{
  const double a = (t < 0.5) ? t / 5 : (1 - t) / 5;

  AmanziGeometry::Point yv(xv);
  yv[0] += a * std::sin(2 * M_PI * xv[0]) * std::sin(2 * M_PI * xv[1]);
  yv[1] += a * std::sin(2 * M_PI * xv[0]) * std::sin(2 * M_PI * xv[1]);
  return yv;
}


/* *****************************************************************
* v = curl u
***************************************************************** */
double
TestRemap(std::function<AmanziGeometry::Point(AmanziGeometry::Point)> fun_u,
          std::function<AmanziGeometry::Point(AmanziGeometry::Point)> fun_v,
          int remap_order = 1,
          int n = 8,
          int mt = 1)
{
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;

  auto comm = Amanzi::getDefaultComm();
  int MyPID = comm->MyPID();
  if (MyPID == 0) std::cout << "\nTest: remap 2-cochains\n";

  // create two meshes
  auto mlist = Teuchos::rcp(new Teuchos::ParameterList());
  mlist->set<bool>("request faces", true);
  mlist->set<bool>("request edges", true);
  mlist->sublist("verbose object").set<std::string>("verbosity level", "none");

  MeshFactory factory(comm, Teuchos::null, mlist);
  factory.set_preference(Preference({Framework::MSTK}));
  Teuchos::RCP<Mesh> mesh_base = factory.create(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, n, n + 1, n + 2);
  Teuchos::RCP<Mesh> mesh0 = factory.create(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, n, n + 1, n + 2);
  Teuchos::RCP<Mesh> mesh1 = factory.create(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, n, n + 1, n + 2);

  int ncells_owned = mesh0->getNumEntities(AmanziMesh::CELL, AmanziMesh::Parallel_kind::OWNED);
  int nfaces_owned = mesh0->getNumEntities(AmanziMesh::FACE, AmanziMesh::Parallel_kind::OWNED);

  int nfaces_wghost = mesh0->getNumEntities(AmanziMesh::FACE, AmanziMesh::Parallel_kind::ALL);
  int nedges_wghost = mesh0->getNumEntities(AmanziMesh::EDGE, AmanziMesh::Parallel_kind::ALL);

  // initialize I/O
  // Teuchos::ParameterList iolist;
  // iolist.get<std::string>("file name base", "plot2");
  // OutputXDMF io(iolist, mesh0, true, true);

  // create and initialize cell-based fields
  CompositeVectorSpace cvs0, cvs1;
  cvs0.SetMesh(mesh0)->SetGhosted(true)->AddComponent("face", AmanziMesh::Entity_kind::FACE, 1);
  cvs1.SetMesh(mesh1)->SetGhosted(true)->AddComponent("face", AmanziMesh::Entity_kind::FACE, 1);

  Teuchos::RCP<CompositeVector> u0 = Teuchos::rcp(new CompositeVector(cvs0));
  Teuchos::RCP<CompositeVector> u1 = Teuchos::rcp(new CompositeVector(cvs1));

  auto& u0_f = *u0->ViewComponent("face", true);
  auto& u1_f = *u1->ViewComponent("face", true);

  // -- edge data
  std::vector<double> u0_e(nedges_wghost);
  for (int e = 0; e < nedges_wghost; ++e) {
    const auto& xe = mesh0->getEdgeCentroid(e);
    const auto& tau = mesh0->getEdgeVector(e);
    u0_e[e] = fun_u(xe) * tau;
  } 

  // -- face data 
  for (int f = 0; f < nfaces_wghost; ++f) {
    const auto& [edges, dirs] = mesh0->getFaceEdgesAndDirections(f);
    int nedges = edges.size();
    for (int n = 0; n < nedges; ++n) {
      int e = edges[n];
      u0_f[0][f] += u0_e[e] * dirs[n];
    }
  }

  // ------- cyclic remap -------
  double l2_err, inf_err;
  double dtau = 1.0 / mt;

  for (int m = 0; m < mt; ++m) {
    Deform3D(mesh0, mesh_base, m * dtau, deform1);
    Deform3D(mesh1, mesh_base, (m + 1) * dtau, deform1);

    RemapCochains remap(mesh0, mesh1);
    if (remap_order == 1)
      remap.LowOrder<2>(*u0, *u1);
    else 
      remap.HighOrder<2>(*u0, *u1);

    // errors
    l2_err = 0.0;
    inf_err = 0.0;
    for (int f = 0; f < nfaces_owned; ++f) {
      const auto& xf = mesh1->getFaceCentroid(f);
      const auto& normal = mesh1->getFaceNormal(f);
      double area = mesh1->getFaceArea(f);

      double err = std::fabs(fun_v(xf) * normal - u1_f[0][f]);
      inf_err = std::max(inf_err, err / area);
      l2_err += err * err / std::sqrt(area);
    }
    l2_err = std::sqrt(l2_err);

    // divergence errors
    double div_err(0.0);
    for (int c = 0; c < ncells_owned; ++c) {
      const auto& [faces, dirs] = mesh0->getCellFacesAndDirections(c);
      int nfaces = faces.size();

      double sum0(0.0), sum1(0.0);
      for (int n = 0; n < nfaces; ++n) {
        int f = faces[n];
        sum0 += u0_f[0][f] * dirs[n];
        sum1 += u1_f[0][f] * dirs[n];
      }
      double vol = mesh1->getCellVolume(c);
      div_err += sum1 * sum1 / (vol * vol);
    }
    div_err = std::sqrt(div_err / ncells_owned);
    CHECK_CLOSE(0.0, div_err, 1e-12);

    // energy errors
    const auto& u0r = remap.getCellFieldFromFaces(mesh0, *u0);
    const auto& u1r = remap.getCellFieldFromFaces(mesh1, *u1);
    auto& u0r_c = *u0r->ViewComponent("cell");
    auto& u1r_c = *u1r->ViewComponent("cell");

    double energy0(0.0), energy1(0.0), energy_base;
    for (int c = 0; c < ncells_owned; ++c) {
      double tmp0(0.0), tmp1(0.0);
      for (int i = 0; i < 3; ++i) {
        tmp0 += u0r_c[i][c] * u0r_c[i][c];
        tmp1 += u1r_c[i][c] * u1r_c[i][c];
      }
      energy0 += tmp0 * mesh0->getCellVolume(c);
      energy1 += tmp1 * mesh1->getCellVolume(c);
    }
    energy0 = std::sqrt(energy0);
    energy1 = std::sqrt(energy1);
    if (m == 0) energy_base = energy0;

    if (MyPID == 0) {
      printf("L2=%12.6e  Inf=%12.6e  Div=%12.6e  E0=%12.6e dE=%12.6e\n",
        l2_err, inf_err, div_err, energy_base, energy1 - energy_base);
    }

    // swap data
    u0_f = u1_f;

    // io
    // io.InitializeCycle(m * dtau, m, "");
    // io.WriteVector(*u0r_c(0), "solution", AmanziMesh::Entity_kind::CELL);
    // io.FinalizeCycle();

    /*
    if (m == 0) {
      std::cout << "# mesh_size: 16 16 16" << std::endl;
      for (int c = 0; c < ncells_owned; ++c) {
        std::cout << mesh0->getCellCentroid(c) << " " << u0r_c[0][c] << " " << u0r_c[2][c] << " " << u0r_c[2][c] << std::endl;
      }
    }
    if (m == mt - 1) {
      std::cout << "# mesh_size: 16 16 16" << std::endl;
      for (int c = 0; c < ncells_owned; ++c) {
        std::cout << mesh1->getCellCentroid(c) << " " << u1r_c[0][c] << " " << u1r_c[2][c] << " " << u1r_c[2][c] << std::endl;
      }
    }
    */
  }

  return l2_err;
}

TEST(REMAP_2_COCHAINS) {
  int low(1), high(2);
  double l2_err = TestRemap(fun0_u, fun0_v, low, 2, 2);
  CHECK_CLOSE(0.0, l2_err, 1e-12);
  l2_err = TestRemap(fun1_u, fun1_v, high, 8, 2);
  CHECK_CLOSE(0.0, l2_err, 2e-3);

  // TestRemap(tornado_u, tornado_v, high, 12, 4);
}
