/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  WhetStone, version 2.0
  Release name: naka-to.

  The nonliner finite volume method.
*/

#include <cstdlib>
#include <cmath>
#include <iostream>

// TPLs
#include "UnitTest++.h"

// Amanzi
#include "MeshFactory.hh"
#include "Mesh.hh"

// WhetStone
#include "nlfv.hh"
#include "Tensor.hh"


/* ****************************************************************
* Test positive decomposition of a given vector in 2D.
**************************************************************** */
TEST(NLFV_POSITIVE_DECOMPOSITION_2D)
{
  using namespace Amanzi;

  std::cout << "Test: positive decomposition of a 2D vector" << std::endl;

  // create basis vectors
  std::vector<AmanziGeometry::Point> tau;
  for (int i = 0; i < 6; i++) {
    double theta = (i + 0.2) * M_PI / 3;
    double a = double(i) / (i + 1.0);
    AmanziGeometry::Point p(cos(theta), a * sin(theta));
    tau.push_back(p);
    std::cout << "tau[" << i << "] = " << p << std::endl;
  }

  int ids[2];
  double ws[2];
  WhetStone::NLFV nlfv;
  AmanziGeometry::Point conormal(2), v(2);

  for (int i = 0; i < 6; i++) {
    conormal = tau[i];
    conormal[1] += 0.1;

    int ierr = nlfv.PositiveDecomposition(i, tau, conormal, ws, ids);

    std::cout << "cornormal = " << conormal << "\nws: " << ws[0] << " " << ws[1]
              << "\nids: " << ids[0] << " " << ids[1] << std::endl;
    v = ws[0] * tau[ids[0]] + ws[1] * tau[ids[1]];

    CHECK_CLOSE(conormal[0], v[0], 1e-12);
    CHECK_CLOSE(conormal[1], v[1], 1e-12);
    CHECK(ierr == 0);
  }
}


/* ****************************************************************
* Test positive decomposition of a given vector in 3D.
**************************************************************** */
TEST(NLFV_POSITIVE_DECOMPOSITION_3D)
{
  using namespace Amanzi;

  std::cout << "\nTest: positive decomposition of a 3D vector" << std::endl;

  // create basis vectors
  int n(0);
  double h[3] = { 0.9, 1.0, 1.2 };
  std::vector<AmanziGeometry::Point> tau;
  for (int i = -1; i < 2; i += 2) {
    for (int j = -1; j < 2; j += 2) {
      for (int k = -1; k < 2; k += 2) {
        AmanziGeometry::Point p(i * h[0], j * h[1], k * h[2]);
        tau.push_back(p);
        std::cout << "tau[" << n++ << "] = " << p << std::endl;
      }
    }
  }

  int ids[3];
  double ws[3];
  WhetStone::NLFV nlfv;
  AmanziGeometry::Point conormal(3), v(3);

  for (int i = 0; i < 6; i++) {
    conormal = tau[i];
    conormal[0] += 0.1 / (i + 1);
    conormal[1] += 0.2 / (i + 1);
    conormal[2] += 0.3 / (i + 1);

    int ierr = nlfv.PositiveDecomposition(i, tau, conormal, ws, ids);

    std::cout << "cornormal = " << conormal << "\nws: " << ws[0] << " " << ws[1] << " " << ws[2]
              << "\nids: " << ids[0] << " " << ids[1] << " " << ids[2] << std::endl;
    v = ws[0] * tau[ids[0]] + ws[1] * tau[ids[1]] + ws[2] * tau[ids[2]];

    CHECK_CLOSE(conormal[0], v[0], 1e-12);
    CHECK_CLOSE(conormal[1], v[1], 1e-12);
    CHECK_CLOSE(conormal[2], v[2], 1e-12);
    CHECK(ierr == 0);
  }
}


/* ****************************************************************
* Test harmonic averaging points in 2D.
**************************************************************** */
TEST(HARMONIC_AVERAGING_POINT_2D)
{
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;

  std::cout << "\nTest: Harmonic averagin point in 2D" << std::endl;
#ifdef HAVE_MPI
  auto comm = Amanzi::getDefaultComm();
#else
  auto comm = Amanzi::getCommSelf();
#endif

  // initialize a two-cell mesh (quad and triangle)
  Preference pref;
  pref.clear();
  pref.push_back(Framework::MSTK);

  MeshFactory meshfactory(comm);
  meshfactory.set_preference(pref);
  Teuchos::RCP<Mesh> mesh = meshfactory.create("test/two_cell2_dist.exo");

  // instantiate the toolset and populate data
  WhetStone::NLFV nlfv(mesh);

  int f(1), c1(0), c2(1);
  double w;
  AmanziGeometry::Point xa(0.8, 0.0), xb(0.7, 1.0); // end-points of face f
  AmanziGeometry::Point p(2), v(2), u(2), xab(2), xcc(2);
  const AmanziGeometry::Point& xc1 = mesh->getCellCentroid(c1);
  const AmanziGeometry::Point& xc2 = mesh->getCellCentroid(c2);

  // identity tensor: conormal = normal
  {
    double tmp1, tmp2;
    AmanziGeometry::Point conormal1(2), conormal2(2);
    conormal1 = mesh->getFaceNormal(f);
    conormal2 = mesh->getFaceNormal(f);

    nlfv.HarmonicAveragingPoint(f, c1, c2, conormal1, conormal2, p, w);
    std::cout << "hap: " << p << " weight=" << w << std::endl;
    v = w * xc1 + (1.0 - w) * xc2;

    // hap is intersection of two lines xa-xb and xc1-xc2
    xab = xb - xa;
    xcc = xc2 - xc1;
    tmp1 = ((xa - xc1) ^ xcc)[0];
    tmp2 = (xab ^ xcc)[0];
    u = xa - xab * (tmp1 / tmp2);

    CHECK(norm(v - p) < 1e-12);
    CHECK(norm(u - p) < 1e-12);
  }

  // rotated normal
  {
    AmanziGeometry::Point conormal1(1.0, 0.2), conormal2(1.0, 0.2);
    nlfv.HarmonicAveragingPoint(f, c1, c2, conormal1, conormal2, p, w);
    std::cout << "hap: " << p << " weight=" << w << "\n\n";

    v = w * xc1 + (1.0 - w) * xc2;
    CHECK(norm(v - p) < 1e-12);

    double tmp = ((p - xa) ^ (xb - xa))[0];
    CHECK(fabs(tmp) < 1e-12);
  }

  // full tensors
  {
    WhetStone::Tensor K1(2, 2), K2(2, 2);
    AmanziGeometry::Point conormal1(2), conormal2(2);
    const AmanziGeometry::Point& normal = mesh->getFaceNormal(f);

    K1(0, 0) = 2.0;
    K1(0, 1) = K1(1, 0) = 1.0;
    K1(1, 1) = 1.0;

    K2(0, 0) = 1.0;
    K2(0, 1) = K2(1, 0) = -1.0;
    K2(1, 1) = 2.0;

    conormal1 = K1 * normal;
    conormal2 = K2 * normal;

    nlfv.HarmonicAveragingPoint(f, c1, c2, conormal1, conormal2, p, w);
    std::cout << "hap: " << p << " weight=" << w << "\n\n";

    double tmp = ((p - xa) ^ (xb - xa))[0];

    CHECK_CLOSE(0.35985587749587, w, 1e-12);
    CHECK(fabs(tmp) < 1e-12);
  }

  // symmetry
  {
    double w1, w2;
    AmanziGeometry::Point p1(2), p2(2);
    AmanziGeometry::Point conormal1(1.0, 0.2), conormal2(1.0, 0.0);

    nlfv.HarmonicAveragingPoint(f, c1, c2, conormal1, conormal2, p1, w1);
    std::cout << "hap: " << p1 << " weight1=" << w1 << std::endl;

    nlfv.HarmonicAveragingPoint(f, c2, c1, conormal2, conormal1, p2, w2);
    std::cout << "hap: " << p2 << " weight2=" << w2 << std::endl;

    CHECK_CLOSE(0.0, norm(p1 - p2), 1e-12);
    CHECK_CLOSE(1.0, w1 + w2, 1e-12);
  }
}
