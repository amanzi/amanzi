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


/* *****************************************************************
* v = curl u
***************************************************************** */
double 
TestRemap(std::function<AmanziGeometry::Point(AmanziGeometry::Point)> fun_u,
          std::function<AmanziGeometry::Point(AmanziGeometry::Point)> fun_v,
          int n = 8,
          double t = 0.02)
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
  Teuchos::RCP<Mesh> mesh0 = factory.create(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, n, n + 1, n + 2);
  Teuchos::RCP<Mesh> mesh1 = factory.create(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, n, n + 1, n + 2);
  DeformMesh(mesh1, 6, t);

  int ncells_owned = mesh0->getNumEntities(AmanziMesh::CELL, AmanziMesh::Parallel_kind::OWNED);
  int nfaces_owned = mesh0->getNumEntities(AmanziMesh::FACE, AmanziMesh::Parallel_kind::OWNED);

  int nfaces_wghost = mesh0->getNumEntities(AmanziMesh::FACE, AmanziMesh::Parallel_kind::ALL);
  int nedges_wghost = mesh0->getNumEntities(AmanziMesh::EDGE, AmanziMesh::Parallel_kind::ALL);

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

  // create remap object
  RemapCochains remap(mesh0, mesh1);
  remap.LowOrder<2>(*u0, *u1);

  // errors
  double l2_err(0.0), inf_err(0.0);
  for (int f = 0; f < nfaces_owned; ++f) {
    const auto& xf = mesh1->getFaceCentroid(f);
    const auto& normal = mesh1->getFaceNormal(f);
    double area = mesh1->getFaceArea(f);

    double err = std::fabs(fun_v(xf) * normal - u1_f[0][f]);
    inf_err = std::max(inf_err, err);
    l2_err = err * err / std::sqrt(area);
  }
  l2_err = std::sqrt(l2_err);

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

  // energies
  const auto& u0r = remap.getCellFieldFromFaces(mesh0, *u0);
  const auto& u1r = remap.getCellFieldFromFaces(mesh1, *u1);

  double energy0(0.0), energy1(0.0);
  for (int c = 0; c < ncells_owned; ++c) {
    energy0 += L22(u0r[c]) * mesh0->getCellVolume(c);
    energy1 += L22(u1r[c]) * mesh1->getCellVolume(c);
  }
  energy0 = std::sqrt(energy0);
  energy1 = std::sqrt(energy1);

  if (MyPID == 0) {
    printf("L2=%12.8g  Inf=%12.8g  Div=%12.8g  E0=%12.8g dE=%12.8g\n",
      l2_err, inf_err, div_err, energy0, energy1 - energy0);
  }

  return l2_err;
}

TEST(REMAP_2_COCHAINS) {
  double l2_err = TestRemap(fun0_u, fun0_v, 2);
  CHECK_CLOSE(0.0, l2_err, 1e-12);
  TestRemap(fun1_u, fun1_v);

  TestRemap(fun2_u, fun2_v, 4, 0.02);
  TestRemap(fun2_u, fun2_v, 8, 0.01);
}
