/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <vector>

#include "UnitTest++.h"

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_LAPACK.hpp"

#include "MeshFactory.hh"
#include "MeshAudit.hh"

#include "Mesh.hh"
#include "Point.hh"

#include "MFD3D_Elasticity.hh"
#include "MFD3D_ElasticityWeakSymmetry.hh"
#include "Tensor.hh"


/* **************************************************************** */
TEST(ELASTICITY_STIFFNESS_2D)
{
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::WhetStone;

  std::cout << "\nTest: Stiffness matrix for elasticity in 2D" << std::endl;
#ifdef HAVE_MPI
  auto comm = Amanzi::getDefaultComm();
#else
  auto comm = Amanzi::getCommSelf();
#endif

  MeshFactory meshfactory(comm);
  meshfactory.set_preference(Preference({ Framework::MSTK }));
  // RCP<Mesh> mesh = meshfactory.create(0.0, 0.0, 1.0, 1.0, 1, 1);
  RCP<Mesh> mesh = meshfactory.create("test/one_pentagon.exo");

  Teuchos::ParameterList plist;
  plist.set<std::string>("base", "cell");
  MFD3D_Elasticity mfd(plist, mesh);

  int nnodes = 5, nrows = nnodes * 2, cell = 0;
  DenseMatrix A(nrows, nrows);

  for (int method = 0; method < 2; method++) {
    Tensor T;
    double lambda(1.0), mu(0.0);
    if (method == 0) {
      T.Init(2, 1);
      T(0, 0) = 1.0;
    } else if (method == 1) {
      mu = 3.0;
      lambda = 0.0;
      T.Init(2, 4);
      T(0, 0) = T(1, 1) = lambda + 2 * mu;
      T(0, 1) = T(1, 0) = lambda;
      T(2, 2) = T(3, 3) = 2 * mu;
    }

    // mfd.StiffnessMatrix(cell, T, A);
    mfd.StiffnessMatrixOptimized(cell, T, A);
    // mfd.StiffnessMatrixMMatrix(cell, T, A);
    // std::cout << "Number of simplex itrs=" << mfd.simplex_num_itrs() << std::endl;
    // std::cout << "Functional value=" << mfd.simplex_functional() << std::endl;

    printf("Stiffness matrix for cell %3d\n", cell);
    PrintMatrix(A, "%8.4f ");

    // verify SPD propery
    for (int i = 0; i < nrows; i++) CHECK(A(i, i) > 0.0);

    // verify exact integration property
    auto nodes = mesh->getCellNodes(cell);

    int d = mesh->getSpaceDimension();
    Point p(d);

    double xx[nrows], yy[nrows];
    for (int i = 0; i < nnodes; i++) {
      int v = nodes[i];
      p = mesh->getNodeCoordinate(v);
      xx[2 * i] = p[0];
      xx[2 * i + 1] = 0.0;

      yy[2 * i] = 0.0;
      yy[2 * i + 1] = p[1];
    }

    double vxx = 0.0, vxy = 0.0, volume = mesh->getCellVolume(cell);
    for (int i = 0; i < nrows; i++) {
      for (int j = 0; j < nrows; j++) {
        vxx += A(i, j) * xx[i] * xx[j];
        vxy += A(i, j) * yy[i] * xx[j];
      }
    }
    CHECK_CLOSE((2 * mu + lambda) * volume, vxx, 1e-10);
    CHECK_CLOSE(0.0, vxy, 1e-10);
  }
}


/* **************************************************************** */
TEST(ELASTICITY_STIFFNESS_3D)
{
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::WhetStone;

  std::cout << "\nTest: Stiffness matrix for Elasticity in 3D" << std::endl;
  auto comm = Amanzi::getDefaultComm();

  MeshFactory meshfactory(comm);
  meshfactory.set_preference(Preference({ Framework::MSTK }));
  // RCP<Mesh> mesh = meshfactory.create(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1, 1, 1);
  RCP<Mesh> mesh = meshfactory.create("test/one_trapezoid.exo");

  Teuchos::ParameterList plist;
  plist.set<std::string>("base", "cell");
  MFD3D_Elasticity mfd(plist, mesh);

  int nrows = 24, nnodes = 8, cell = 0;
  Tensor T(3, 1);
  T(0, 0) = 1;

  DenseMatrix A(nrows, nrows);
  mfd.StiffnessMatrix(cell, T, A);

  printf("Stiffness matrix for cell %3d\n", cell);
  PrintMatrix(A, "%8.4f ");

  // verify SPD propery
  for (int i = 0; i < nrows; i++) CHECK(A(i, i) > 0.0);

  // verify exact integration property
  auto nodes = mesh->getCellNodes(cell);

  int d = mesh->getSpaceDimension();
  Point p(d);

  double xx[nrows], yy[nrows];
  for (int i = 0; i < nnodes; i++) {
    int v = nodes[i];
    p = mesh->getNodeCoordinate(v);
    xx[3 * i] = p[0];
    xx[3 * i + 1] = 0.0;
    xx[3 * i + 2] = 0.0;

    yy[3 * i] = 0.0;
    yy[3 * i + 1] = p[1];
    yy[3 * i + 2] = 0.0;
  }

  double vxx = 0.0, vxy = 0.0, volume = mesh->getCellVolume(cell);
  for (int i = 0; i < nrows; i++) {
    for (int j = 0; j < nrows; j++) {
      vxx += A(i, j) * xx[i] * xx[j];
      vxy += A(i, j) * xx[i] * yy[j];
    }
  }
  CHECK_CLOSE(vxx, volume, 1e-10);
  CHECK_CLOSE(vxy, 0.0, 1e-10);
}


/* **************************************************************** */
TEST(ELASTICITY_WEAK_SYMMETRY_2D)
{
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::WhetStone;

  std::cout << "\nTest: Mass matrix for elasticity with weak symmetry in 2D" << std::endl;
#ifdef HAVE_MPI
  auto comm = Amanzi::getDefaultComm();
#else
  auto comm = Amanzi::getCommSelf();
#endif

  MeshFactory meshfactory(comm);
  meshfactory.set_preference(Preference({ Framework::MSTK }));
  // RCP<Mesh> mesh = meshfactory.create(0.0, 0.0, 1.0, 1.0, 1, 1);
  RCP<Mesh> mesh = meshfactory.create("test/one_pentagon.exo");

  Teuchos::ParameterList plist;
  plist.set<std::string>("base", "cell");
  MFD3D_ElasticityWeakSymmetry mfd(plist, mesh);

  DenseMatrix M, G, A;
  for (int method = 0; method < 2; method++) {
    Tensor T;
    double lambda(0.0), mu(1.0);
    if (method == 0) {
      mu = 2.0;
      lambda = 0.0;
      T.Init(2, 1);
      T(0, 0) = mu;
    } else if (method == 1) {
      mu = 3.0;
      lambda = 1.0;
      T.Init(2, 4);
      T(0, 0) = T(1, 1) = lambda + mu;
      T(0, 1) = T(1, 0) = lambda;
      T(2, 2) = T(3, 3) = mu;
    }
    Tensor Tinv(T);
    Tinv.Inverse();

    mfd.MassMatrix(0, Tinv, M);
    mfd.RotationMatrix(0, G);
    mfd.StiffnessMatrix(0, Tinv, A);

    printf("Mass matrix for cell 0\n");
    PrintMatrix(M, "%8.4f ", 6);

    // verify SPD propery
    int nrows = M.NumRows();
    for (int i = 0; i < nrows; i++) CHECK(M(i, i) > 0.0);

    int mrows = A.NumRows();
    for (int i = 0; i < mrows; i++) CHECK(A(i, i) > 0.0);

    // verify exact integration property for mass matrix
    int orientation;
    int nfaces = mesh->getNumEntities(AmanziMesh::Entity_kind::FACE, Parallel_kind::ALL);
    std::vector<double> xx(nrows), yy(nrows), rr(nrows);

    Tensor Tx(2, 2), Ty(2, 2);
    Tx(0, 0) = 1.0;
    Ty(1, 0) = 1.0;
    Tensor T1 = T * Tx;
    Tensor T2 = T * Ty;

    const auto& faces = mesh->getCellFaces(0);
    for (int i = 0; i < nfaces; i++) {
      int f = faces[i];
      const auto& normal = mesh->getFaceNormal(f, 0, &orientation);
      double area = mesh->getFaceArea(f);

      const auto p1 = T1 * normal;
      xx[2 * i] = p1[0] / area;
      xx[2 * i + 1] = p1[1] / area;

      const auto p2 = T2 * normal;
      yy[2 * i] = p2[0] / area;
      yy[2 * i + 1] = p2[1] / area;

      const auto p3 = Ty * normal;
      rr[2 * i] = p3[0] / area;
      rr[2 * i + 1] = p3[1] / area;
    }

    double vxx(0.0), vxy(0.0), rot(0.0);
    double volume = mesh->getCellVolume(0);
    for (int i = 0; i < nrows; i++) {
      rot += G(i, 0) * rr[i];
      for (int j = 0; j < nrows; j++) {
        vxx += M(i, j) * xx[i] * xx[j];
        vxy += M(i, j) * yy[i] * xx[j];
      }
    }
    CHECK_CLOSE(-volume, rot, 1e-10);
    CHECK_CLOSE((mu + lambda) * volume, vxx, 1e-10);
    CHECK_CLOSE(0.0, vxy, 1e-10);

    // verify exact integration property for mass matrix
    std::vector<double> vx(mrows);
    for (int i = 0; i < nfaces; i++) {
      int f = faces[i];
      const auto& xf = mesh->getFaceCentroid(f);
      vx[2 * i] = 1.0;
      vx[2 * i + 1] = 1.0;
    }

    const auto& xc = mesh->getCellCentroid(0);
    vx[4 * nfaces] = 1.0;
    vx[4 * nfaces + 1] = 1.0;

    double axx(0.0);
    for (int i = 0; i < mrows; i++) {
      for (int j = 0; j < mrows; j++) {
        axx += A(i, j) * vx[i] * vx[j];
      }
    }
    CHECK_CLOSE(0.0, axx, 1e-10);
  }
}


/* **************************************************************** */
void RunWeakSymmetry3D(const std::string& filename)
{
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::WhetStone;

  std::cout << "\nTest: 3D mass matrix for elasticity with weak symmetry: " 
            << filename << std::endl;
#ifdef HAVE_MPI
  auto comm = Amanzi::getDefaultComm();
#else
  auto comm = Amanzi::getCommSelf();
#endif

  auto fac_list = Teuchos::rcp(new ParameterList());
  fac_list->set<bool>("request edges", true);
  MeshFactory meshfactory(comm, Teuchos::null, fac_list);
  meshfactory.set_preference(Preference({ Framework::MSTK }));
  // RCP<Mesh> mesh = meshfactory.create(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1, 1, 1);
  RCP<Mesh> mesh = meshfactory.create(filename);

  Teuchos::ParameterList plist;
  plist.set<std::string>("base", "cell");
  MFD3D_ElasticityWeakSymmetry mfd(plist, mesh);

  DenseMatrix M, G, A;
  for (int method = 0; method < 2; method++) {
    Tensor T;
    double lambda, mu;
    if (method == 0) {
      mu = 2.0;
      lambda = 0.0;
      T.Init(3, 1);
      T(0, 0) = mu;
    } else if (method == 1) {
      mu = 3.0;
      lambda = 1.0;
      T.Init(3, 4);
      for (int i = 0; i < 9; ++i) T(i, i) = mu;
      T(0, 0) = T(1, 1) = T(2, 2) += lambda;
      T(0, 1) = T(1, 0) = lambda;
    }
    Tensor Tinv(T);
    Tinv.Inverse();

    mfd.MassMatrix(0, Tinv, M);
    mfd.RotationMatrix(0, G);
    mfd.StiffnessMatrix(0, Tinv, A);

    printf("Mass matrix for cell 0\n");
    PrintMatrix(M, "%8.4f ", 12);

    // verify SPD propery
    int nrows = M.NumRows();
    for (int i = 0; i < nrows; i++) CHECK(M(i, i) > 0.0);

    // verify exact integration property
    int orientation;
    int nfaces = mesh->getNumEntities(AmanziMesh::Entity_kind::FACE, Parallel_kind::ALL);
    std::vector<double> xx(nrows), yy(nrows), zz(nrows), rr(nrows);


    Tensor Tx(3, 2), Ty(3, 2), Tz(3, 2);
    Tx(0, 0) = 1.0;
    Ty(1, 0) = 1.0;
    Tz(2, 2) = 1.0;
    Tensor T1 = T * Tx;
    Tensor T2 = T * Ty;
    Tensor T3 = T * Tz;

    const auto& faces = mesh->getCellFaces(0);
    for (int i = 0; i < nfaces; i++) {
      int f = faces[i];
      const auto& normal = mesh->getFaceNormal(f, 0, &orientation);
      double area = mesh->getFaceArea(f);

      auto p = T1 * normal;
      for (int k = 0; k < 3; ++k) xx[3 * i + k] = p[k] / area;

      p = T2 * normal;
      for (int k = 0; k < 3; ++k) yy[3 * i + k] = p[k] / area;

      p = T3 * normal;
      for (int k = 0; k < 3; ++k) zz[3 * i + k] = p[k] / area;

      p = Ty * normal;
      for (int k = 0; k < 3; ++k) rr[3 * i + k] = p[k] / area;
    }

    double vxx(0.0), vxy(0.0), vzz(0.0), rot(0.0);
    double volume = mesh->getCellVolume(0);

    for (int i = 0; i < nrows; i++) {
      rot += G(i, 0) * rr[i];
      for (int j = 0; j < nrows; j++) {
        vxx += M(i, j) * xx[i] * xx[j];
        vxy += M(i, j) * yy[i] * xx[j];
        vzz += M(i, j) * zz[i] * zz[j];
      }
    }
    CHECK_CLOSE(-volume, rot, 1e-10);
    CHECK_CLOSE((mu + lambda) * volume, vxx, 1e-10);
    CHECK_CLOSE(0.0, vxy, 1e-10);
    CHECK_CLOSE((mu + lambda) * volume, vzz, 1e-10);
  }
}

TEST(ELASTICITY_WEAK_SYMMETRY_3D_UNIT_CUBE) {
  RunWeakSymmetry3D("test/cube_unit.exo");
}

TEST(ELASTICITY_WEAK_SYMMETRY_3D_UNIT_CUBE_ROTATED) {
  RunWeakSymmetry3D("test/cube_unit_rotated.exo");
}

TEST(ELASTICITY_WEAK_SYMMETRY_3D_PARALLEPIPED_ROTATED) {
  RunWeakSymmetry3D("test/parallepiped_rotated.exo");
}
