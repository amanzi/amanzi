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
#include "Tensor.hh"

// Amanzi::Operators
#include "MeshDeformation.hh"

using namespace Amanzi;

static AmanziGeometry::Point fun0_v(const AmanziGeometry::Point& p) {
 return AmanziGeometry::Point(1.0, 2.0, 3.0);
}
static double fun0_u(const AmanziGeometry::Point& p) {
 double x = p[0];
 double y = p[1];
 double z = p[2];
 return x + 2 * y + 3 * z;
}

static AmanziGeometry::Point fun1_v(const AmanziGeometry::Point& p) {
 return AmanziGeometry::Point(p[0], 2 * p[1], -3 * p[2]);
}
static double fun1_u(const AmanziGeometry::Point& p) {
 double x = p[0];
 double y = p[1];
 double z = p[2];
 return 0.5 * x * x + y * y - 1.5 * z * z;
}

static AmanziGeometry::Point fun2_v(const AmanziGeometry::Point& p) {
 double x = p[0];
 double y = p[1];
 double z = p[2];
 double val = std::cos(x * y * z);
 return val * AmanziGeometry::Point(y * z, x * z, x * y);
}
static double fun2_u(const AmanziGeometry::Point& p) {
 double x = p[0];
 double y = p[1];
 double z = p[2];
 return std::sin(x * y * z);
}


/* *****************************************************************
* v = grad u
***************************************************************** */
double 
TestRemap(std::function<double(AmanziGeometry::Point)> fun_u,
          std::function<AmanziGeometry::Point(AmanziGeometry::Point)> fun_v,
          int n = 8)
{
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;

  auto comm = Amanzi::getDefaultComm();
  int MyPID = comm->MyPID();
  if (MyPID == 0) std::cout << "\nTest: remap 1-cochains\n";

  // create two meshes
  auto mlist = Teuchos::rcp(new Teuchos::ParameterList());
  mlist->set<bool>("request faces", true);
  mlist->set<bool>("request edges", true);
  mlist->sublist("verbose object").set<std::string>("verbosity level", "none");

  MeshFactory factory(comm, Teuchos::null, mlist);
  factory.set_preference(Preference({Framework::MSTK}));
  Teuchos::RCP<Mesh> mesh0 = factory.create(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, n, n + 1, n + 2);
  Teuchos::RCP<Mesh> mesh1 = factory.create(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, n, n + 1, n + 2);
  DeformMesh(mesh1, 6, 0.01);

  int ncells_owned = mesh0->getNumEntities(AmanziMesh::CELL, AmanziMesh::Parallel_kind::OWNED);
  int nfaces_owned = mesh0->getNumEntities(AmanziMesh::FACE, AmanziMesh::Parallel_kind::OWNED);
  int nedges_owned = mesh0->getNumEntities(AmanziMesh::EDGE, AmanziMesh::Parallel_kind::OWNED);

  int nedges_wghost = mesh0->getNumEntities(AmanziMesh::EDGE, AmanziMesh::Parallel_kind::ALL);
  int nnodes_wghost = mesh0->getNumEntities(AmanziMesh::NODE, AmanziMesh::Parallel_kind::ALL);

  // check planarity
  double err, err_avg(0.0), err_max(0.0);
  WhetStone::Tensor T(3, 2);
  for (int f = 0; f < nfaces_owned; ++f) {
    const auto& nodes = mesh1->getFaceNodes(f);
    const auto& x0 = mesh1->getNodeCoordinate(nodes[0]);
    const auto& x1 = mesh1->getNodeCoordinate(nodes[1]);
    const auto& x2 = mesh1->getNodeCoordinate(nodes[2]);
    const auto& x3 = mesh1->getNodeCoordinate(nodes[3]);

    T.SetColumn(0, x0 - x1);
    T.SetColumn(1, x0 - x2);
    T.SetColumn(2, x0 - x3);

    err = std::fabs(T.Det() / (norm(x0 - x1) * norm(x0 - x2) * norm(x0 - x3)));
    err_max = std::max(err_max, err); 
    err_avg += err;
  }
  if (err_max > 1e-12) {
    printf("Mesh face planarity error: mean=%12.8g  max=%12.8g\n", err_avg / nfaces_owned, err_max);
  }

  // create and initialize cell-based fields
  CompositeVectorSpace cvs0, cvs1;
  cvs0.SetMesh(mesh0)->SetGhosted(true)->AddComponent("edge", AmanziMesh::Entity_kind::EDGE, 1);
  cvs1.SetMesh(mesh1)->SetGhosted(true)->AddComponent("edge", AmanziMesh::Entity_kind::EDGE, 1);

  Teuchos::RCP<CompositeVector> u0 = Teuchos::rcp(new CompositeVector(cvs0));
  Teuchos::RCP<CompositeVector> u1 = Teuchos::rcp(new CompositeVector(cvs1));

  auto& u0_e = *u0->ViewComponent("edge", true);
  auto& u1_e = *u1->ViewComponent("edge", true);

  // -- node data
  std::vector<double> u0_n(nnodes_wghost);
  for (int n = 0; n < nnodes_wghost; ++n) {
    const auto& xn = mesh0->getNodeCoordinate(n);
    u0_n[n] = fun_u(xn);
  } 

  // -- edge data 
  for (int e = 0; e < nedges_wghost; ++e) {
    const auto& nodes = mesh0->getEdgeNodes(e);
    u0_e[0][e] = u0_n[nodes[1]] - u0_n[nodes[0]];
  }

  // create remap object
  RemapCochains remap(mesh0, mesh1);
  remap.LowOrder<1>(*u0, *u1);

  // errors
  double l2_err(0.0), inf_err(0.0);
  for (int e = 0; e < nedges_owned; ++e) {
    const auto& xe = mesh1->getEdgeCentroid(e);
    const auto& tau = mesh1->getEdgeVector(e);
    double len = mesh1->getEdgeLength(e);

    double err = std::fabs(fun_v(xe) * tau - u1_e[0][e]);
    inf_err = std::max(inf_err, err);
    l2_err = len * err * err;
  }
  l2_err = std::sqrt(l2_err);

  double curl_err(0.0);
  for (int f = 0; f < nfaces_owned; ++f) {
    const auto& [edges, dirs] = mesh0->getFaceEdgesAndDirections(f);
    int nedges = edges.size();

    double sum0(0.0), sum1(0.0);
    for (int n = 0; n < nedges; ++n) {
      int e = edges[n];
      sum0 += u0_e[0][e] * dirs[n];
      sum1 += u1_e[0][e] * dirs[n];
    }
    double area = mesh1->getFaceArea(f);
    curl_err += sum1 * sum1 / (area * area);
  }
  curl_err = std::sqrt(curl_err / nfaces_owned);
  CHECK_CLOSE(0.0, curl_err, 1e-12);

  // energies
  const auto& u0r = remap.getCellFieldFromEdges(mesh0, *u0);
  const auto& u1r = remap.getCellFieldFromEdges(mesh1, *u1);
  auto u0r_c = *u0r->ViewComponent("cell");
  auto u1r_c = *u1r->ViewComponent("cell");

  double energy0(0.0), energy1(0.0);
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

  if (MyPID == 0) {
    printf("L2=%12.8g  Inf=%12.8g  Curl=%12.8g  E0=%12.8g dE=%12.8g\n",
      l2_err, inf_err, curl_err, energy0, energy1 - energy0);
  }

  return l2_err;
}

TEST(REMAP_1_COCHAINS) {
  double l2_err = TestRemap(fun0_u, fun0_v);
  CHECK_CLOSE(0.0, l2_err, 1e-12);
  TestRemap(fun1_u, fun1_v);

  TestRemap(fun2_u, fun2_v, 4);
  TestRemap(fun2_u, fun2_v, 8);
}
