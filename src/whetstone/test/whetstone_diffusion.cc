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

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "UnitTest++.h"

#include "Mesh.hh"
#include "MeshFactory.hh"
#include "MeshAudit.hh"
#include "Point.hh"

#include "MFD3D_CrouzeixRaviart.hh"
#include "MFD3D_Diffusion.hh"
#include "MFD3D_GeneralizedDiffusion.hh"
#include "MFD3D_Lagrange.hh"
#include "Tensor.hh"


/* **************************************************************** */
TEST(DARCY_MASS_2D)
{
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::WhetStone;

  std::cout << "Test: Mass matrix for Darcy in 2D" << std::endl;
  auto comm = Amanzi::getDefaultComm();

  MeshFactory meshfactory(comm);
  meshfactory.set_preference(Preference({ Framework::MSTK }));
  Teuchos::RCP<Mesh> mesh = meshfactory.create(0.0, 0.0, 0.5, 1.0, 1, 1);

  MFD3D_Diffusion mfd(mesh);
  DeRham_Face drc(mfd);

  int nfaces = 4, cell = 0;
  DenseMatrix M;

  for (int method = 0; method < 3; method++) {
    Tensor T(2, 2);
    T(0, 0) = 1.0;
    T(1, 1) = 1.0;
    T(0, 1) = 0.1;
    T(1, 0) = 0.1;

    M.Reshape(nfaces, nfaces);
    if (method == 0) {
      mfd.MassMatrix(cell, T, M);
    } else if (method == 1) {
      T(0, 1) += 0.1;
      mfd.MassMatrixNonSymmetric(cell, T, M);
    } else if (method == 2) {
      drc.MassMatrix(cell, T, M);
    }

    printf("Mass matrix for cell %3d\n", cell);
    PrintMatrix(M, "%8.4f ");

    // verify SPD propery
    for (int i = 0; i < nfaces; i++) CHECK(M(i, i) > 0.0);

    // verify exact integration property
    auto [faces, dirs] = mesh->getCellFacesAndDirections(cell);

    double xi, yi, xj;
    double vxx = 0.0, vxy = 0.0, volume = mesh->getCellVolume(cell);
    for (int i = 0; i < nfaces; i++) {
      int f1 = faces[i];
      for (int j = 0; j < nfaces; j++) {
        int f2 = faces[j];

        xi = mesh->getFaceNormal(f1)[0] * dirs[i];
        yi = mesh->getFaceNormal(f1)[1] * dirs[i];
        xj = mesh->getFaceNormal(f2)[0] * dirs[j];

        if (method == 0 || method == 2) {
          xi /= mesh->getFaceArea(f1);
          yi /= mesh->getFaceArea(f1);
          xj /= mesh->getFaceArea(f2);
        }

        vxx += M(i, j) * xi * xj;
        vxy += M(i, j) * yi * xj;
      }
    }

    CHECK_CLOSE(T(0, 0) * volume, vxx, 1e-10);
    CHECK_CLOSE(T(1, 0) * volume, vxy, 1e-10);
  }
}


/* **************************************************************** */
TEST(DARCY_MASS_3D)
{
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::WhetStone;

  std::cout << "\n\nTest: Mass matrix for Darcy in 3D" << std::endl;
#ifdef HAVE_MPI
  auto comm = Amanzi::getDefaultComm();
#else
  auto comm = Amanzi::getCommSelf();
#endif

  Preference pref;
  pref.clear();
  pref.push_back(Framework::SIMPLE);

  MeshFactory meshfactory(comm);
  meshfactory.set_preference(pref);
  Teuchos::RCP<Mesh> mesh = meshfactory.create(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1, 1, 1);

  MFD3D_Diffusion mfd(mesh);
  DeRham_Face drc(mfd);

  int nfaces = 6, cell = 0;
  DenseMatrix M;

  for (int method = 0; method < 2; method++) {
    Tensor T(3, 1);
    T(0, 0) = 1.0;

    if (method == 0) {
      mfd.MassMatrix(cell, T, M);
    } else if (method == 1) {
      drc.MassMatrix(cell, T, M);
    }

    printf("Mass matrix for cell %3d\n", cell);
    PrintMatrix(M, "%8.4f ");

    // verify SPD propery
    for (int i = 0; i < nfaces; ++i) CHECK(M(i, i) > 0.0);

    // verify exact integration property
    auto [faces, dirs] = mesh->getCellFacesAndDirections(cell);

    double xi, xj, yj;
    double vxx = 0.0, vxy = 0.0, volume = mesh->getCellVolume(cell);
    for (int i = 0; i < nfaces; i++) {
      int f1 = faces[i];
      for (int j = 0; j < nfaces; j++) {
        int f2 = faces[j];

        xi = mesh->getFaceNormal(f1)[0] * dirs[i];
        xj = mesh->getFaceNormal(f2)[0] * dirs[j];
        yj = mesh->getFaceNormal(f2)[1] * dirs[j];

        if (method == 1) {
          xi /= mesh->getFaceArea(f1);
          xj /= mesh->getFaceArea(f2);
          yj /= mesh->getFaceArea(f2);
        }

        vxx += M(i, j) * xi * xj;
        vxy += M(i, j) * xi * yj;
      }
    }

    CHECK_CLOSE(volume, vxx, 1e-10);
    CHECK_CLOSE(0.0, vxy, 1e-10);
  }
}


/* **************************************************************** */
TEST(DARCY_MASS_3D_GENERALIZED_POLYHEDRON)
{
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::WhetStone;

  std::cout << "\n\nTest: Mass matrix for generalized polyhedra" << std::endl;
  auto comm = Amanzi::getDefaultComm();

  MeshFactory meshfactory(comm);
  meshfactory.set_preference(Preference({ Framework::MSTK }));
  Teuchos::RCP<Mesh> mesh = meshfactory.create("test/hex_random.exo");

  Teuchos::ParameterList plist;
  MFD3D_GeneralizedDiffusion mfd(plist, mesh);

  int nfaces = 6, cell = 0;
  double volume = mesh->getCellVolume(cell);
  DenseMatrix N, R, M;
  DenseMatrix B(3, 3);

  Tensor T(3, 2);
  T(0, 0) = 0.5;
  T(1, 1) = 0.344827586206897;
  T(2, 2) = 0.103448275862069;
  T(1, 2) = T(2, 1) = -0.03448275862069;

  // consistency condition
  mfd.L2consistency(cell, T, N, M, true);
  mfd.L2consistencyInverse(cell, T, R, M, true);

  B.Multiply(N, R, true);
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      double tmp = (i == j) ? volume : 0.0;
      CHECK_CLOSE(B(i, j), tmp, 1e-6);
    }
  }

  // mass matrix
  mfd.MassMatrix(cell, T, M);

  printf("Mass matrix for cell %3d  volume=%12.4f\n", cell, volume);
  PrintMatrix(M, "%8.4f ");

  // verify SPD propery
  for (int i = 0; i < nfaces; ++i) CHECK(M(i, i) > 0.0);
}


/* **************************************************************** */
TEST(DARCY_INVERSE_MASS_3D)
{
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::WhetStone;

  std::cout << "\n\nTest: Inverse mass matrix for Darcy" << std::endl;
#ifdef HAVE_MPI
  auto comm = Amanzi::getDefaultComm();
#else
  auto comm = Amanzi::getCommSelf();
#endif

  Preference pref;
  pref.clear();
  pref.push_back(Framework::SIMPLE);

  MeshFactory meshfactory(comm);
  meshfactory.set_preference(pref);
  Teuchos::RCP<Mesh> mesh = meshfactory.create(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1, 2, 3);

  MFD3D_Diffusion mfd(mesh);
  DeRham_Face drc(mfd);

  int nfaces = 6, cell = 0;
  Tensor T(3, 1); // tensor of rank 1
  T(0, 0) = 1.0;

  DenseMatrix W;
  for (int method = 0; method < 5; method++) {
    if (method == 0) {
      mfd.MassMatrixInverse(cell, T, W);
    } else if (method == 1) {
      mfd.MassMatrixInverseOptimized(cell, T, W);
    } else if (method == 2) {
      mfd.MassMatrixInverseSO(cell, T, W);
    } else if (method == 3) {
      mfd.MassMatrixInverseMMatrixHex(cell, T, W);
    } else if (method == 4) {
      drc.MassMatrixInverse(cell, T, W);
    }

    printf("Inverse of mass matrix for method=%d\n", method);
    PrintMatrix(W, "%8.4f ");

    // verify SPD some propeties
    for (int i = 0; i < nfaces; i++) CHECK(W(i, i) > 0.0);

    // verify exact integration property
    W.Inverse();

    auto [faces, dirs] = mesh->getCellFacesAndDirections(cell);

    double xi, yi, xj;
    double vxx = 0.0, vxy = 0.0, volume = mesh->getCellVolume(cell);
    for (int i = 0; i < nfaces; i++) {
      int f1 = faces[i];
      for (int j = 0; j < nfaces; j++) {
        int f2 = faces[j];

        xi = mesh->getFaceNormal(f1)[0] * dirs[i];
        yi = mesh->getFaceNormal(f1)[1] * dirs[i];
        xj = mesh->getFaceNormal(f2)[0] * dirs[j];

        if (method == 4) {
          xi /= mesh->getFaceArea(f1);
          yi /= mesh->getFaceArea(f1);
          xj /= mesh->getFaceArea(f2);
        }

        vxx += W(i, j) * xi * xj;
        vxy += W(i, j) * yi * xj;
      }
    }
    CHECK_CLOSE(vxx, volume, 1e-10);
    CHECK_CLOSE(vxy, 0.0, 1e-10);
  }
}


/* **************************************************************** */
TEST(DARCY_FULL_TENSOR_2D)
{
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::WhetStone;

  std::cout << "\n\nTest: Inverse mass matrix and full tensor in 2D" << std::endl;
#ifdef HAVE_MPI
  auto comm = Amanzi::getDefaultComm();
#else
  auto comm = Amanzi::getCommSelf();
#endif

  MeshFactory meshfactory(comm);
  meshfactory.set_preference(Preference({ Framework::MSTK }));
  Teuchos::RCP<Mesh> mesh = meshfactory.create("test/two_cell2.exo");

  DenseMatrix W;
  MFD3D_Diffusion mfd(mesh);

  for (int cell = 0; cell < 2; cell++) {
    for (int method = 0; method < 7; method++) {
      Tensor T(2, 2); // tensor of rank 2
      T(0, 0) = 1.0;
      T(1, 1) = 2.0;
      T(0, 1) = T(1, 0) = 1.0;

      int ok, nfaces = mesh->getCellNumFaces(cell);

      if (method == 0) {
        mfd.MassMatrixInverse(cell, T, W);
      } else if (method == 1) {
        mfd.MassMatrixInverseOptimized(cell, T, W);
      } else if (method == 2) {
        mfd.MassMatrixInverseSO(cell, T, W);
      } else if (method == 3) {
        if (nfaces != 4) continue;
        mfd.MassMatrixInverseMMatrixHex(cell, T, W);
      } else if (method == 4) {
        ok = mfd.MassMatrixInverseMMatrix(cell, T, W);
        std::cout << "Number of simplex itrs=" << mfd.simplex_num_itrs() << " code=" << ok
                  << std::endl;
        std::cout << "Functional value=" << mfd.simplex_functional() << std::endl;
      } else if (method == 5) {
        double kmean = 1.0;
        AmanziGeometry::Point kgrad(0.1, 0.2);
        mfd.MassMatrixInverseDivKScaled(cell, T, kmean, kgrad, W);
      } else if (method == 6) {
        T(0, 1) += 0.1;
        mfd.MassMatrixInverseNonSymmetric(cell, T, W);
      }

      printf("Inverse of mass matrix for method=%d\n", method);
      PrintMatrix(W, "%8.4f ");

      // verify PD propery
      for (int i = 0; i < nfaces; i++) CHECK(W(i, i) > 0.0);

      // verify exact integration property
      W.Inverse();

      auto [faces, dirs] = mesh->getCellFacesAndDirections(cell);

      AmanziGeometry::Point v(1.0, 2.0);
      double xi, xj;
      double vxx = 0.0, volume = mesh->getCellVolume(cell);
      for (int i = 0; i < nfaces; i++) {
        xi = (v * mesh->getFaceNormal(faces[i])) * dirs[i];
        for (int j = 0; j < nfaces; j++) {
          xj = (v * mesh->getFaceNormal(faces[j])) * dirs[j];
          vxx += W(i, j) * xi * xj;
        }
      }
      CHECK_CLOSE(2 * volume, vxx, 1e-10);

      // additional tests for triangle: integral with v2 = RT basis function
      if (cell == 1) {
        if (method == 2 || method == 5) continue;
        vxx = 0.0;
        for (int i = 0; i < nfaces; i++) {
          xi = (v * mesh->getFaceNormal(faces[i])) * dirs[i];
          for (int j = 0; j < nfaces; j++) {
            xj = ((j == 2) ? mesh->getFaceArea(faces[j]) : 0.0) * dirs[j];
            vxx += W(i, j) * xi * xj;
          }
        }
        double vxx_exact = (method == 5) ? -0.054167 : -0.05;
        CHECK_CLOSE(vxx_exact, vxx, 1e-10);
      }
    }
  }
}


/* **************************************************************** */
TEST(DARCY_FULL_TENSOR_3D)
{
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::WhetStone;

  std::cout << "\nTest: Inverse mass matrix and full tensor in 3D" << std::endl;
#ifdef HAVE_MPI
  auto comm = Amanzi::getDefaultComm();
#else
  auto comm = Amanzi::getCommSelf();
#endif

  Preference pref;
  pref.clear();
  pref.push_back(Framework::SIMPLE);

  MeshFactory meshfactory(comm);
  meshfactory.set_preference(pref);
  Teuchos::RCP<Mesh> mesh = meshfactory.create(0.0, 0.0, 0.0, 1.1, 1.0, 1.0, 3, 2, 1);

  MFD3D_Diffusion mfd(mesh);
  DeRham_Face drc(mfd);

  int nfaces = 6, cell = 0;
  Tensor T(3, 2); // tensor of rank 2
  T(0, 0) = 1.0;
  T(1, 1) = 2.0;
  T(2, 2) = 3.0;
  T(0, 1) = T(1, 0) = 1.0;
  T(1, 2) = T(2, 1) = 1.0;

  DenseMatrix W;
  for (int method = 0; method < 6; method++) {
    if (method == 0) {
      mfd.MassMatrixInverse(cell, T, W);
    } else if (method == 1) {
      mfd.MassMatrixInverseOptimized(cell, T, W);
    } else if (method == 2) {
      mfd.MassMatrixInverseSO(cell, T, W);
    } else if (method == 3) {
      mfd.MassMatrixInverseMMatrixHex(cell, T, W);
    } else if (method == 4) {
      mfd.MassMatrixInverseMMatrix(cell, T, W);
      std::cout << "Number of simplex itrs=" << mfd.simplex_num_itrs() << std::endl;
      std::cout << "Functional value=" << mfd.simplex_functional() << std::endl;
    } else if (method == 5) {
      drc.MassMatrixInverse(cell, T, W);
    }

    printf("Inverse of mass matrix for method=%d\n", method);
    PrintMatrix(W, "%8.4f ");

    // verify SPD propery
    for (int i = 0; i < nfaces; i++) CHECK(W(i, i) > 0.0);

    // verify exact integration property
    W.Inverse();

    auto [faces, dirs] = mesh->getCellFacesAndDirections(cell);

    AmanziGeometry::Point v(1.0, 2.0, 3.0);
    double xi, xj;
    double vxx = 0.0, volume = mesh->getCellVolume(cell);
    for (int i = 0; i < nfaces; i++) {
      int f1 = faces[i];
      for (int j = 0; j < nfaces; j++) {
        int f2 = faces[j];

        xi = (v * mesh->getFaceNormal(f1)) * dirs[i];
        xj = (v * mesh->getFaceNormal(f2)) * dirs[j];

        if (method == 5) {
          xi /= mesh->getFaceArea(f1);
          xj /= mesh->getFaceArea(f2);
        }

        vxx += W(i, j) * xi * xj;
      }
    }
    CHECK_CLOSE(vxx, 4 * volume, 1e-10);
  }
}


/* **************************************************************** */
TEST(DARCY_STIFFNESS_2D_NODE)
{
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::WhetStone;

  std::cout << "\nTest: Stiffness matrix for Darcy in 2D" << std::endl;
#ifdef HAVE_MPI
  auto comm = Amanzi::getDefaultComm();
#else
  auto comm = Amanzi::getCommSelf();
#endif

  MeshFactory meshfactory(comm);
  meshfactory.set_preference(Preference({ Framework::MSTK }));
  RCP<Mesh> mesh = meshfactory.create("test/one_pentagon.exo");

  MFD3D_Diffusion mfd(mesh);

  int nnodes = 5, cell = 0;
  Tensor T(2, 1);
  T(0, 0) = 1;

  DenseMatrix A;
  for (int method = 0; method < 3; method++) {
    if (method == 0) {
      mfd.StiffnessMatrix(cell, T, A);
    } else if (method == 1) {
      mfd.StiffnessMatrixMMatrix(cell, T, A);
      std::cout << "Number of simplex itrs=" << mfd.simplex_num_itrs() << std::endl;
      std::cout << "Functional value=" << mfd.simplex_functional() << std::endl;
    } else if (method == 2) {
      mfd.StiffnessMatrixOptimized(cell, T, A);
    }

    printf("Stiffness matrix for cell %3d\n", cell);
    PrintMatrix(A, "%8.4f ");

    // verify SPD propery
    for (int i = 0; i < nnodes; i++) CHECK(A(i, i) > 0.0);

    // verify exact integration property
    AmanziMesh::Entity_ID_List nodes;
    nodes = mesh->getCellNodes(cell);

    int d = mesh->getSpaceDimension();
    Point p(d);

    double xi, yi, xj;
    double vxx = 0.0, vxy = 0.0, volume = mesh->getCellVolume(cell);
    for (int i = 0; i < nnodes; i++) {
      int v = nodes[i];
      p = mesh->getNodeCoordinate(v);
      xi = p[0];
      yi = p[1];
      for (int j = 0; j < nnodes; j++) {
        v = nodes[j];
        p = mesh->getNodeCoordinate(v);
        xj = p[0];
        vxx += A(i, j) * xi * xj;
        vxy += A(i, j) * yi * xj;
      }
    }
    CHECK_CLOSE(vxx, volume, 1e-10);
    CHECK_CLOSE(vxy, 0.0, 1e-10);
  }
}


/* **************************************************************** */
TEST(DARCY_STIFFNESS_2D_EDGE)
{
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::WhetStone;

  std::cout << "\nTest: Stiffness matrix for Darcy in 2D:edges" << std::endl;
#ifdef HAVE_MPI
  auto comm = Amanzi::getDefaultComm();
#else
  auto comm = Amanzi::getCommSelf();
#endif

  Teuchos::RCP<Teuchos::ParameterList> factory_plist = Teuchos::rcp(new Teuchos::ParameterList());
  factory_plist->set<bool>("request edges", true);
  factory_plist->set<bool>("request faces", true);
  MeshFactory meshfactory(comm, Teuchos::null, factory_plist);
  meshfactory.set_preference(Preference({ Framework::MSTK }));
  RCP<Mesh> mesh = meshfactory.create("test/one_pentagon.exo");

  Teuchos::ParameterList plist;
  plist.set<int>("method order", 1);
  MFD3D_CrouzeixRaviart mfd(plist, mesh);

  int nedges = 5, cell = 0;
  Tensor T(2, 1);
  T(0, 0) = 1;

  DenseMatrix A;
  for (int method = 0; method < 1; method++) {
    if (method == 0) { mfd.StiffnessMatrix(cell, T, A); }

    printf("Stiffness matrix for cell %3d\n", cell);
    PrintMatrix(A, "%8.4f ");

    // verify SPD propery
    for (int i = 0; i < nedges; i++) CHECK(A(i, i) > 0.0);

    // verify exact integration property
    auto edges = mesh->getCellEdges(cell);

    int d = mesh->getSpaceDimension();
    Point p(d);

    double xi, yi, xj;
    double vxx = 0.0, vxy = 0.0, volume = mesh->getCellVolume(cell);
    for (int i = 0; i < nedges; i++) {
      int e = edges[i];
      const AmanziGeometry::Point& xe = mesh->getEdgeCentroid(e);
      xi = xe[0];
      yi = xe[1];
      for (int j = 0; j < nedges; j++) {
        e = edges[j];
        const AmanziGeometry::Point& ye = mesh->getEdgeCentroid(e);
        xj = ye[0];
        vxx += A(i, j) * xi * xj;
        vxy += A(i, j) * yi * xj;
      }
    }
    CHECK_CLOSE(vxx, volume, 1e-10);
    CHECK_CLOSE(vxy, 0.0, 1e-10);
  }
}


/* **************************************************************** */
TEST(DARCY_STIFFNESS_3D)
{
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::WhetStone;

  std::cout << "\nTest: Stiffness matrix for Darcy in 3D" << std::endl;
#ifdef HAVE_MPI
  auto comm = Amanzi::getDefaultComm();
#else
  auto comm = Amanzi::getCommSelf();
#endif

  MeshFactory meshfactory(comm);
  meshfactory.set_preference(Preference({ Framework::MSTK }));
  RCP<Mesh> mesh = meshfactory.create("test/dodecahedron.exo");

  MFD3D_Diffusion mfd(mesh);

  int nnodes = 20, cell = 0;
  Tensor T(3, 1);
  T(0, 0) = 1.0;

  DenseMatrix A;
  mfd.StiffnessMatrixMMatrix(cell, T, A);

  printf("Stiffness matrix for cell %3d\n", cell);
  PrintMatrix(A, "%8.4f ");
  std::cout << "Number of simplex itrs=" << mfd.simplex_num_itrs() << std::endl;
  std::cout << "Functional value=" << mfd.simplex_functional() << std::endl;

  // verify SPD propery
  for (int i = 0; i < nnodes; i++) CHECK(A(i, i) > 0.0);

  // verify exact integration property
  AmanziMesh::Entity_ID_List nodes;
  nodes = mesh->getCellNodes(cell);

  int d = mesh->getSpaceDimension();
  Point p(d);

  double xi, yi, xj;
  double vxx = 0.0, vxy = 0.0, volume = mesh->getCellVolume(cell);
  for (int i = 0; i < nnodes; i++) {
    int v = nodes[i];
    p = mesh->getNodeCoordinate(v);
    xi = p[0];
    yi = p[1];
    for (int j = 0; j < nnodes; j++) {
      v = nodes[j];
      p = mesh->getNodeCoordinate(v);
      xj = p[0];
      vxx += A(i, j) * xi * xj;
      vxy += A(i, j) * yi * xj;
    }
  }
  CHECK_CLOSE(vxx, volume, 1e-10);
  CHECK_CLOSE(vxy, 0.0, 1e-10);
}


/* **************************************************************** */
TEST(RECOVER_GRADIENT_MIXED)
{
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::WhetStone;

  std::cout << "\nTest: Recover gradient from Darcy fluxes" << std::endl;
#ifdef HAVE_MPI
  auto comm = Amanzi::getDefaultComm();
#else
  auto comm = Amanzi::getCommSelf();
#endif

  MeshFactory meshfactory(comm);
  meshfactory.set_preference(Preference({ Framework::MSTK }));
  RCP<Mesh> mesh = meshfactory.create("test/one_trapezoid.exo");

  MFD3D_Diffusion mfd(mesh);

  // create Darcy fluxes
  int nfaces = 6, cell = 0;
  auto faces = mesh->getCellFaces(cell);

  Point flux(1.0, 2.0, 3.0);
  std::vector<Polynomial> solution(nfaces);

  for (int n = 0; n < nfaces; n++) {
    int f = faces[n];
    const Point& normal = mesh->getFaceNormal(f);
    solution[n].Reshape(3, 0);
    solution[n](0) = -normal * flux;
  }

  // gradient recovery
  Polynomial gradient(3, 1);
  mfd.L2Cell(cell, solution, solution, NULL, gradient);

  printf("Gradient %f %f %f\n", gradient(1), gradient(2), gradient(3));

  CHECK_CLOSE(gradient(1), 1.0, 1e-10);
  CHECK_CLOSE(gradient(2), 2.0, 1e-10);
  CHECK_CLOSE(gradient(3), 3.0, 1e-10);
}


/* **************************************************************** */
TEST(RECOVER_GRADIENT_NODAL)
{
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::WhetStone;

  std::cout << "\nTest: Recover gradient from nodal pressures" << std::endl;
#ifdef HAVE_MPI
  auto comm = Amanzi::getDefaultComm();
#else
  auto comm = Amanzi::getCommSelf();
#endif

  MeshFactory meshfactory(comm);
  meshfactory.set_preference(Preference({ Framework::MSTK }));
  RCP<Mesh> mesh = meshfactory.create("test/one_trapezoid.exo");

  Teuchos::ParameterList plist;
  plist.set<int>("method order", 1).set<bool>("use low-order scheme", true);
  MFD3D_Lagrange mfd(plist, mesh);

  // create pressure solution
  AmanziMesh::Entity_ID_List nodes;
  int nnodes = 8, cell = 0;
  nodes = mesh->getCellNodes(cell);

  Point slope(1.0, 2.0, 3.0);
  std::vector<Polynomial> solution(nnodes);
  Point xv(3);

  for (int n = 0; n < nnodes; n++) {
    int v = nodes[n];
    xv = mesh->getNodeCoordinate(v);
    solution[n].Reshape(3, 0);
    solution[n](0) = slope * xv;
  }

  // gradient recovery
  WhetStone::Polynomial gradient(3, 1);
  mfd.L2Cell(cell, solution, solution, NULL, gradient);

  printf("Gradient %f %f %f\n", gradient(1), gradient(2), gradient(3));

  CHECK_CLOSE(gradient(1), 1.0, 1e-10);
  CHECK_CLOSE(gradient(2), 2.0, 1e-10);
  CHECK_CLOSE(gradient(3), 3.0, 1e-10);
}


/* **************************************************************** */
TEST(DARCY_INVERSE_MASS_2D)
{
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::WhetStone;

  std::cout << "\nTest: Inverse mass matrix for Darcy, 2D" << std::endl;
#ifdef HAVE_MPI
  auto comm = Amanzi::getDefaultComm();
#else
  auto comm = Amanzi::getCommSelf();
#endif

  MeshFactory meshfactory(comm);
  meshfactory.set_preference(Preference({ Framework::MSTK }));
  RCP<Mesh> mesh = meshfactory.create("test/dodecahedron.exo");

  MFD3D_Diffusion mfd(mesh);

  int ok, nfaces = 12, cell = 0, dim = mesh->getSpaceDimension();
  Tensor T(dim, 2); // tensor of rank 1
  T(0, 0) = 1.0;
  T(1, 1) = 1.0;
  T(2, 2) = 1.0;
  T(0, 1) = T(1, 0) = 0.0;

  DenseMatrix W(nfaces, nfaces);
  for (int method = 0; method < 1; method++) {
    if (method == 0) {
      ok = mfd.MassMatrixInverseMMatrix(cell, T, W);
      std::cout << "Number of simplex itrs=" << mfd.simplex_num_itrs() << std::endl;
      std::cout << "Functional value=" << mfd.simplex_functional() << std::endl;
    }

    printf("Inverse of mass matrix for method=%d  ierr=%d\n", method, ok);
    PrintMatrix(W, "%8.4f ");

    // verify monotonicity propery
    for (int i = 0; i < nfaces; i++) {
      CHECK(W(i, i) > 0.0);
      for (int j = i + 1; j < nfaces; j++) CHECK(W(i, j) < 1e-10);
    }

    // verify exact integration property
    W.Inverse();
    T.Inverse();

    auto [faces, dirs] = mesh->getCellFacesAndDirections(cell);

    double xi, yi, xj;
    double vxx = 0.0, vxy = 0.0, volume = mesh->getCellVolume(cell);
    for (int i = 0; i < nfaces; i++) {
      int f = faces[i];
      xi = mesh->getFaceNormal(f)[0] * dirs[i];
      yi = mesh->getFaceNormal(f)[1] * dirs[i];
      for (int j = 0; j < nfaces; j++) {
        f = faces[j];
        xj = mesh->getFaceNormal(f)[0] * dirs[j];
        vxx += W(i, j) * xi * xj;
        vxy += W(i, j) * yi * xj;
      }
    }
    CHECK_CLOSE(volume * T(0, 0), vxx, 1e-10);
    CHECK_CLOSE(volume * T(0, 1), vxy, 1e-10);
    CHECK_CLOSE(mfd.simplex_functional(), 60.0, 1e-2);
  }
}
