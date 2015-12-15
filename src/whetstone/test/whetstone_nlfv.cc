/*
  WhetStone, version 2.0
  Release name: naka-to.

  Copyright 2010-2012 held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

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
TEST(NLFV_POSITIVE_DECOMPOSITION_2D) {
  using namespace Amanzi;

  std::cout << "Test: positive decomposition of a vector" << std::endl;
 
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

    std::cout << "cornormal = " << conormal
              << "\nws: " << ws[0] << " " << ws[1]
              << "\nids: " << ids[0] << " " << ids[1] << std::endl;
    v = ws[0] * tau[ids[0]] + ws[1] * tau[ids[1]];

    CHECK_CLOSE(conormal[0], v[0], 1e-12);
    CHECK_CLOSE(conormal[1], v[1], 1e-12);
    CHECK(ierr == 0);
  }
}


/* ****************************************************************
* Test harmonic averaging points in 2D.
**************************************************************** */
TEST(HARMONIC_AVERAGING_POINT_2D) {
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;

  std::cout << "\nTest: Harmonic averagin point in 2D" << std::endl;
#ifdef HAVE_MPI
  Epetra_MpiComm *comm = new Epetra_MpiComm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm *comm = new Epetra_SerialComm();
#endif

  // initialize a two-cell mesh (quad and triangle)
  FrameworkPreference pref;
  pref.clear();
  pref.push_back(MSTK);

  MeshFactory meshfactory(comm);
  meshfactory.preference(pref);
  Teuchos::RCP<Mesh> mesh = meshfactory("test/two_cell2.exo"); 
 
  // instantiate the toolset and populate data
  WhetStone::NLFV nlfv(mesh);

  int f(1), c1(0), c2(1);
  double w;
  AmanziGeometry::Point p(2), v(2);
  const AmanziGeometry::Point& xf = mesh->face_centroid(f);
  const AmanziGeometry::Point& xc1 = mesh->cell_centroid(c1);
  const AmanziGeometry::Point& xc2 = mesh->cell_centroid(c2);

  // identity tensor: conormal = normal
  {
    AmanziGeometry::Point conormal1(2), conormal2(2);
    conormal1 = mesh->face_normal(f);
    conormal2 = mesh->face_normal(f);

    nlfv.HarmonicAveragingPoint(f, c1, c2, conormal1, conormal2, p, w);
    std::cout << "hap: " << p << " weight=" << w << std::endl;
    v = w * xc1 + (1.0 - w) * xc2;

    CHECK_CLOSE(xf[0], p[0], 1e-12);
    CHECK(norm(v - p) < 1e-12);
  }

  // rotated normal
  {
    AmanziGeometry::Point conormal1(1.0, 0.2), conormal2(1.0, 0.2);
    nlfv.HarmonicAveragingPoint(f, c1, c2, conormal1, conormal2, p, w);
    std::cout << "hap: " << p << " weight=" << w << std::endl;
    v = w * xc1 + (1.0 - w) * xc2;

    CHECK_CLOSE(xf[0], p[0], 1e-12);
    CHECK(norm(v - p) < 1e-12);
  }

  // full tensors
  {
    WhetStone::Tensor K1(2, 2), K2(2, 2); 
    AmanziGeometry::Point conormal1(2), conormal2(2);
    const AmanziGeometry::Point& normal = mesh->face_normal(f);

    K1(0, 0) = 2.0;
    K1(0, 1) = K1(1, 0) = 1.0;
    K1(1, 1) = 1.0;

    K2(0, 0) = 1.0;
    K2(0, 1) = K2(1, 0) = -1.0;
    K2(1, 1) = 2.0;

    conormal1 = K1 * normal;
    conormal2 = K2 * normal;

    nlfv.HarmonicAveragingPoint(f, c1, c2, conormal1, conormal2, p, w);
    std::cout << "hap: " << p << " weight=" << w << std::endl;

    CHECK_CLOSE(0.363636363636364, w, 1e-12);
    CHECK_CLOSE(xf[0], p[0], 1e-12);
  }
}

