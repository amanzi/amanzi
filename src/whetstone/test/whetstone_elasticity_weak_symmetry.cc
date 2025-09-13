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
#include "MFD3D_ElasticityWeakSym.hh"
#include "MFD3D_ElasticityWeakSymBdV.hh"
#include "NumericalIntegration.hh"
#include "Tensor.hh"


using namespace Teuchos;
using namespace Amanzi;
using namespace Amanzi::AmanziGeometry;
using namespace Amanzi::AmanziMesh;
using namespace Amanzi::WhetStone;


/* **************************************************************** */
TEST(ELASTICITY_WEAK_SYMMETRY_2D)
{
  std::cout << "\nTEST: 2D matrices for elasticity with weak symmetry" << std::endl;
  auto comm = Amanzi::getDefaultComm();

  MeshFactory meshfactory(comm);
  meshfactory.set_preference(Preference({ Framework::MSTK }));
  // RCP<Mesh> mesh = meshfactory.create(0.0, 0.0, 1.0, 1.0, 1, 1);
  RCP<Mesh> mesh = meshfactory.create("test/one_pentagon.exo");

  Teuchos::ParameterList plist;
  plist.set<std::string>("base", "cell");
  MFD3D_ElasticityWeakSym mfd(plist, mesh);

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

    printf("Mass matrix for cell 0, method=%d size=%d\n", method, M.NumRows());
    PrintMatrix(M, "%8.4f ");
    printf("Stiffness matrix for cell 0, method=%d size=%d\n", method, A.NumRows());
    PrintMatrix(A, "%10.6f ");

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
      vx[2 * i] = xf[0];
      vx[2 * i + 1] = xf[1];
    }

    const auto& xc = mesh->getCellCentroid(0);
    vx[2 * nfaces] = xc[0];
    vx[2 * nfaces + 1] = xc[1];

    double axx(0.0);
    for (int i = 0; i < mrows; i++) {
      for (int j = 0; j < mrows; j++) {
        axx += A(i, j) * vx[i] * vx[j];
      }
    }
    CHECK_CLOSE((4 * lambda + 2 * mu) * volume, axx, 1e-10);
  }
}


/* **************************************************************** */
TEST(ELASTICITY_WEAK_SYMMETRY_3D)
{
  std::cout << "\nTest: 3D matrices for elasticity with weak symmetry" << std::endl;
  auto comm = Amanzi::getDefaultComm();

  MeshFactory meshfactory(comm);
  meshfactory.set_preference(Preference({ Framework::MSTK }));
  RCP<Mesh> mesh = meshfactory.create(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1, 1, 1);

  Teuchos::ParameterList plist;
  plist.set<std::string>("base", "cell");
  MFD3D_ElasticityWeakSym mfd(plist, mesh);

  DenseMatrix M, G, A;
  for (int method = 0; method < 2; method++) {
    Tensor T;
    double lambda(0.0), mu(1.0);
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
      T(0, 0) = T(1, 1) = T(2, 2) = lambda + mu;
      T(0, 1) = T(1, 0) = T(0, 2) = T(2, 0) = T(1, 2) = T(2, 1) = lambda;
    }
    Tensor Tinv(T);
    Tinv.Inverse();

    mfd.MassMatrix(0, Tinv, M);
    mfd.RotationMatrix(0, G);
    mfd.StiffnessMatrix(0, Tinv, A);

    printf("Mass matrix for cell 0, method=%d size=%d\n", method, M.NumRows());
    PrintMatrix(M, "%8.4f ", 12);
    printf("Stiffness matrix for cell 0, method=%d size=%d\n", method, A.NumRows());
    PrintMatrix(A, "%10.6f ", 12);

    // verify SPD propery
    int nrows = M.NumRows();
    for (int i = 0; i < nrows; i++) CHECK(M(i, i) > 0.0);

    int mrows = A.NumRows();
    for (int i = 0; i < mrows; i++) CHECK(A(i, i) > 0.0);

    // verify exact integration property for mass matrix
    int orientation;
    int nfaces = mesh->getNumEntities(AmanziMesh::Entity_kind::FACE, Parallel_kind::ALL);
    std::vector<double> xx(nrows), yy(nrows), rr(nrows);

    Tensor Tx(3, 2), Ty(3, 2);
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
      xx[3 * i] = p1[0] / area;
      xx[3 * i + 1] = p1[1] / area;
      xx[3 * i + 2] = p1[2] / area;

      const auto p2 = T2 * normal;
      yy[3 * i] = p2[0] / area;
      yy[3 * i + 1] = p2[1] / area;
      yy[3 * i + 2] = p2[2] / area;

      const auto p3 = Ty * normal;
      rr[3 * i] = p3[0] / area;
      rr[3 * i + 1] = p3[1] / area;
      rr[3 * i + 2] = p3[2] / area;
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

    // verify exact integration property for stiffness matrix
    std::vector<double> vx(mrows);
    for (int i = 0; i < nfaces; i++) {
      int f = faces[i];
      const auto& xf = mesh->getFaceCentroid(f);
      vx[3 * i] = xf[0];
      vx[3 * i + 1] = xf[1];
      vx[3 * i + 2] = xf[2];
    }

    const auto& xc = mesh->getCellCentroid(0);
    vx[3 * nfaces] = xc[0];
    vx[3 * nfaces + 1] = xc[1];
    vx[3 * nfaces + 2] = xc[2];

    double axx(0.0);
    for (int i = 0; i < mrows; i++) {
      for (int j = 0; j < mrows; j++) {
        axx += A(i, j) * vx[i] * vx[j];
      }
    }
    CHECK_CLOSE((9 * lambda + 3 * mu) * volume, axx, 1e-10);
  }
}


/* **************************************************************** */
TEST(ELASTICITY_WEAK_SYMMETRY_BDV_2D)
{
  std::cout << "\nTEST: BdV matrices for elasticity with weak symmetry in 2D" << std::endl;
  auto comm = Amanzi::getDefaultComm();

  MeshFactory meshfactory(comm);
  meshfactory.set_preference(Preference({ Framework::MSTK }));
  // RCP<Mesh> mesh = meshfactory.create(0.0, 0.0, 0.1, 0.1, 1, 1);
  RCP<Mesh> mesh = meshfactory.create("test/one_pentagon.exo");
  // RCP<Mesh> mesh = meshfactory.create("test/one_quad.exo");

  Teuchos::ParameterList plist;
  plist.set<std::string>("base", "cell");
  MFD3D_ElasticityWeakSymBdV mfd(plist, mesh);

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

    printf("Mass matrix for cell 0, method=%d size=%d\n", method, M.NumRows());
    PrintMatrix(M, "%8.4f ", 6);

    // verify SPD propery
    int nrows = M.NumRows();
    for (int i = 0; i < nrows; i++) CHECK(M(i, i) > 0.0);

    int mrows = A.NumRows();
    for (int i = 0; i < mrows; i++) CHECK(A(i, i) > 0.0);

    // verify exact integration property for mass matrix
    int orientation;
    int nfaces = mesh->getCellNumFaces(0);
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
      xx[4 * i] = p1[0] / area;
      xx[4 * i + 1] = p1[1] / area;

      const auto p2 = T2 * normal;
      yy[4 * i] = p2[0] / area;
      yy[4 * i + 1] = p2[1] / area;

      const auto p3 = Ty * normal;
      rr[4 * i] = p3[0] / area;
      rr[4 * i + 1] = p3[1] / area;
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

    // verify exact integration property for stiffness matrix
    std::vector<double> vx(mrows);
    for (int i = 0; i < nfaces; i++) {
      int f = faces[i];
      double area = mesh->getFaceArea(f);
      const auto& xf = mesh->getFaceCentroid(f);
      const auto& normal = mesh->getFaceNormal(f);
      vx[4 * i] = 0.0;
      vx[4 * i + 1] = xf[1];
      vx[4 * i + 2] = 0.0;
      vx[4 * i + 3] = normal[0] / 12;
    }

    const auto& xc = mesh->getCellCentroid(0);
    vx[4 * nfaces] = 0.0;
    vx[4 * nfaces + 1] = xc[1];

    double axx(0.0);
    for (int i = 0; i < mrows; i++) {
      for (int j = 0; j < mrows; j++) {
        axx += A(i, j) * vx[i] * vx[j];
      }
    }

    CHECK_CLOSE((lambda + mu) * volume, axx, 1e-10);
  }
}


/* **************************************************************** */
void
RunWeakSymmetryBdV3D(const std::string& filename)
{
  std::cout << "\nTEST: 3D BdV matrices for elasticity with weak symmetry: " << filename
            << std::endl;
#ifdef HAVE_MPI
  auto comm = Amanzi::getDefaultComm();
#else
  auto comm = Amanzi::getCommSelf();
#endif

  auto fac_list = Teuchos::rcp(new ParameterList());
  fac_list->set<bool>("request edges", true);
  MeshFactory meshfactory(comm, Teuchos::null, fac_list);
  meshfactory.set_preference(Preference({ Framework::MSTK }));
  // RCP<Mesh> mesh = meshfactory.create(0.0, 0.0, 0.0, 0.1, 0.1, 0.1, 1, 1, 1);
  RCP<Mesh> mesh = meshfactory.create(filename);

  Teuchos::ParameterList plist;
  plist.set<std::string>("base", "cell");
  MFD3D_ElasticityWeakSymBdV mfd(plist, mesh);

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

    printf("Mass matrix for cell 0, method=%d\n", method);
    PrintMatrix(M, "%8.4f ", 12);

    // verify SPD propery
    int nrows = M.NumRows();
    for (int i = 0; i < nrows; i++) CHECK(M(i, i) > 0.0);

    int mrows = A.NumRows();
    for (int i = 0; i < mrows; i++) CHECK(A(i, i) > 0.0);

    // verify exact integration property for mass matrix
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
      for (int k = 0; k < 3; ++k) xx[9 * i + k] = p[k] / area;

      p = T2 * normal;
      for (int k = 0; k < 3; ++k) yy[9 * i + k] = p[k] / area;

      p = T3 * normal;
      for (int k = 0; k < 3; ++k) zz[9 * i + k] = p[k] / area;

      p = Ty * normal;
      for (int k = 0; k < 3; ++k) rr[9 * i + k] = p[k] / area;
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

    // verify exact integration property for stiffness matrix
    int index[3] = { 0, 0, 0 };
    NumericalIntegration numi(mesh);
    std::vector<const PolynomialBase*> polys(2);

    std::vector<double> vx(mrows);
    for (int i = 0; i < nfaces; i++) {
      int f = faces[i];
      const auto& xf = mesh->getFaceCentroid(f);
      double area = mesh->getFaceArea(f);

      AmanziGeometry::Point normal = mesh->getFaceNormal(f);
      auto coordsys = std::make_shared<AmanziGeometry::SurfaceCoordinateSystem>(xf, normal);

      vx[9 * i] = xf[0];
      vx[9 * i + 1] = 1.0;

      index[0] = 1;
      Polynomial cmono(3, index, 1.0);
      index[0] = 0;
      polys[0] = &cmono;

      index[0] = 1;
      Polynomial fmono1(2, index, 1.0);
      index[0] = 0;
      fmono1.InverseChangeCoordinates(xf, *coordsys->tau());
      polys[1] = &fmono1;
      vx[9 * i + 3] = numi.IntegratePolynomialsFace(f, polys) / area;

      index[1] = 1;
      Polynomial fmono2(2, index, 1.0);
      index[1] = 0;
      fmono2.InverseChangeCoordinates(xf, *coordsys->tau());
      polys[1] = &fmono2;
      vx[9 * i + 6] = numi.IntegratePolynomialsFace(f, polys) / area;
    }

    const auto& xc = mesh->getCellCentroid(0);
    vx[9 * nfaces] = xc[0];
    vx[9 * nfaces + 1] = 1.0;

    double axx(0.0);
    for (int i = 0; i < mrows; i++) {
      for (int j = 0; j < mrows; j++) {
        axx += A(i, j) * vx[i] * vx[j];
      }
    }
    CHECK_CLOSE((lambda + mu) * volume, axx, 1e-10);
  }
}

TEST(ELASTICITY_WEAK_SYMMETRY_3D_UNIT_CUBE)
{
  RunWeakSymmetryBdV3D("test/cube_unit.exo");
}

TEST(ELASTICITY_WEAK_SYMMETRY_3D_UNIT_CUBE_ROTATED)
{
  RunWeakSymmetryBdV3D("test/cube_unit_rotated.exo");
}

TEST(ELASTICITY_WEAK_SYMMETRY_3D_PARALLEPIPED_ROTATED)
{
  RunWeakSymmetryBdV3D("test/parallepiped_rotated.exo");
}
