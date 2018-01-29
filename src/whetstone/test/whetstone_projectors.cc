/*
  WhetStone

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <vector>

#include "Teuchos_RCP.hpp"
#include "UnitTest++.h"

#include "Mesh.hh"
#include "MeshFactory.hh"
#include "Point.hh"

#include "MFD3D_CrouzeixRaviart.hh"
#include "MFD3D_Diffusion.hh"
#include "MFD3D_Lagrange.hh"
#include "NumericalIntegration.hh"
#include "Tensor.hh"
#include "Polynomial.hh"
#include "Projector.hh"


/* **************************************************************** */
TEST(HARMONIC_PROJECTORS_SQUARE_CR) {
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::WhetStone;

  std::cout << "\nTest: HO nonconforming projectors for square (linear deformation)" << std::endl;
  Epetra_MpiComm *comm = new Epetra_MpiComm(MPI_COMM_WORLD);

  MeshFactory meshfactory(comm);
  meshfactory.preference(FrameworkPreference({MSTK}));
  Teuchos::RCP<const Amanzi::AmanziGeometry::GeometricModel> gm;
  Teuchos::RCP<Mesh> mesh = meshfactory(0.0, 0.0, 1.2, 1.1, 1, 1, gm, true, true); 
 
  int cell(0);
  VectorPolynomial uc;
  std::vector<VectorPolynomial> vf(4);

  // test zero cell deformation
  for (int n = 0; n < 4; ++n) {
    vf[n].resize(2);
    for (int i = 0; i < 2; ++i) {
      vf[n][i].Reshape(2, 1, true);
    }
  }

  Projector projector(mesh);
  auto moments = std::make_shared<WhetStone::DenseVector>();

  projector.HarmonicCell_CRk(cell, 2, vf, moments, uc);

  CHECK(uc[0].NormMax() < 1e-12 && uc[1].NormMax() < 1e-12);

  // test linear deformation
  for (int n = 0; n < 4; ++n) {
    for (int i = 0; i < 2; ++i) {
      vf[n][i](0, 0) = 1.0;
      vf[n][i](1, 0) = 2.0;
      vf[n][i](1, 1) = 3.0;
    }
  }
  
  projector.HarmonicCell_CR1(cell, vf, uc);
  std::cout << uc[0] << std::endl;

  uc[0] -= vf[0][0];
  uc[1] -= vf[0][1];
  CHECK(uc[0].NormMax() < 1e-12 && uc[1].NormMax() < 1e-12);

  for (int k = 2; k < 4; ++k) {
    projector.HarmonicCell_CRk(cell, k, vf, moments, uc);
    std::cout << uc[1] << std::endl;

    uc[0] -= vf[0][0];
    uc[1] -= vf[0][1];
    CHECK(uc[0].NormMax() < 1e-12 && uc[1].NormMax() < 1e-12);
  }

  // test re-location of the right-top corner to (2,3)
  std::cout << "Test: HO nonconforming projectors for square (bilinear deformation)" << std::endl;
  for (int n = 0; n < 4; ++n) vf[n].PutScalar(0.0);
  vf[1][0](1, 1) = 0.8 / 1.1; 
  vf[1][1](1, 1) = 1.9 / 1.1; 

  vf[2][0](1, 0) = 0.8 / 1.2; 
  vf[2][1](1, 0) = 1.9 / 1.2; 

  projector.HarmonicCell_CRk(cell, 2, vf, moments, uc);

  std::cout << uc[0] << std::endl;
  std::cout << uc[1] << std::endl;
  auto p = AmanziGeometry::Point(1.2, 1.1);
  CHECK(fabs(uc[0].Value(p) - 0.8) < 1e-12 &&
        fabs(uc[1].Value(p) - 1.9) < 1e-12);

  // add additive constant deformation 
  std::cout << "Test: HO nonconforming projectors for square (2nd bilinear deformation)" << std::endl;
  for (int n = 0; n < 4; ++n) vf[n][0](0, 0) = 1.0;
  projector.HarmonicCell_CRk(cell, 2, vf, moments, uc);
  std::cout << uc[0] << std::endl;

  projector.HarmonicCell_CRk(cell, 3, vf, moments, uc);
  std::cout << uc[0] << std::endl;

  // test harmonic functions (is it always true for square?)
  std::cout << "\nTest: HO nonconforming projectors for square (harmonic function)" << std::endl;
  for (int n = 0; n < 4; ++n) {
    for (int i = 0; i < 2; ++i) {
      vf[n][i].Reshape(2, 2, false);
      vf[n][i](2, 0) = 4.0;
      vf[n][i](2, 1) = 5.0;
      vf[n][i](2, 2) = 6.0;
    }
  }
  projector.HarmonicCell_CRk(cell, 2, vf, moments, uc);
  std::cout << uc[0] << std::endl;

  CHECK(fabs(uc[0](2, 0) + uc[0](2, 2)) < 1e-12 &&
        fabs(uc[1](2, 0) + uc[1](2, 2)) < 1e-12);

  delete comm;
}


/* **************************************************************** */
TEST(HARMONIC_PROJECTORS_POLYGON_CR) {
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::WhetStone;

  std::cout << "\nTest: HO nonconforming projectors for pentagon (linear deformation)" << std::endl;
  Epetra_MpiComm *comm = new Epetra_MpiComm(MPI_COMM_WORLD);

  MeshFactory meshfactory(comm);
  meshfactory.preference(FrameworkPreference({MSTK}));
  Teuchos::RCP<const Amanzi::AmanziGeometry::GeometricModel> gm;
  Teuchos::RCP<Mesh> mesh = meshfactory("test/one_pentagon.exo", gm, true, true);
  // Teuchos::RCP<Mesh> mesh = meshfactory("test/one_quad.exo", gm, true, true);
 
  int cell(0), nfaces(5);
  VectorPolynomial uc;
  std::vector<VectorPolynomial> vf(nfaces);

  // test linear deformation
  for (int n = 0; n < nfaces; ++n) {
    vf[n].resize(2);
    for (int i = 0; i < 2; ++i) {
      vf[n][i].Reshape(2, 1, true);
      vf[n][i](0, 0) = 1.0;
      vf[n][i](1, 0) = 2.0;
      vf[n][i](1, 1) = 3.0;
    }
  }
  
  // -- old scheme
  Projector projector(mesh);
  auto moments = std::make_shared<WhetStone::DenseVector>();

  projector.HarmonicCell_CR1(cell, vf, uc);
  std::cout << uc[0] << std::endl;

  uc[0] -= vf[0][0];
  uc[1] -= vf[0][1];
  CHECK(uc[0].NormMax() < 1e-12 && uc[1].NormMax() < 1e-12);

  // -- new scheme (k=1)
  for (int k = 1; k < 4; ++k) {
    projector.HarmonicCell_CRk(cell, k, vf, moments, uc);
    std::cout << uc[0] << std::endl;

    uc[0] -= vf[0][0];
    uc[1] -= vf[0][1];
    CHECK(uc[0].NormMax() < 1e-12 && uc[1].NormMax() < 1e-12);
  }

  // test quadratic deformation
  std::cout << "\nTest: HO nonconforming for pentagon (quadratic deformation)" << std::endl;
  for (int n = 0; n < nfaces; ++n) {
    for (int i = 0; i < 2; ++i) {
      vf[n][i].Reshape(2, 2, false);
      vf[n][i](2, 0) = 4.0;
      vf[n][i](2, 1) = 5.0;
      vf[n][i](2, 2) = -4.0;
    }
  }

  for (int k = 2; k < 4; ++k) {
    projector.HarmonicCell_CRk(cell, k, vf, moments, uc);
    std::cout << uc[0] << std::endl;
    uc[0] -= vf[0][0];
    uc[1] -= vf[0][1];
    CHECK(uc[0].NormMax() < 1e-10 && uc[1].NormMax() < 1e-10);
  }

  // test cubic deformation
  std::cout << "\nTest: HO nonconforming projectors for pentagon (cubic deformation)" << std::endl;
  for (int n = 0; n < nfaces; ++n) {
    for (int i = 0; i < 2; ++i) {
      vf[n][i].Reshape(2, 3, false);
      vf[n][i](3, 0) = 2.0;
      vf[n][i](3, 1) = -6.0;
      vf[n][i](3, 2) = -6.0;
      vf[n][i](3, 3) = 2.0;
    }
  }

  for (int k = 3; k < 4; ++k) {
    projector.HarmonicCell_CRk(cell, k, vf, moments, uc);
    std::cout << vf[0][0] << std::endl;
    std::cout << uc[0] << std::endl;
    uc[0] -= vf[0][0];
    uc[1] -= vf[0][1];
    CHECK(uc[0].NormMax() < 1e-11 && uc[1].NormMax() < 1e-11);
  }

  // test trace compatibility between function and its projecton (k < 3 only!)
  std::cout << "\nTest: HO non-conforming projectors for pentagon (harmonic function)" << std::endl;
  for (int n = 0; n < nfaces; ++n) {
    for (int i = 0; i < 2; ++i) {
      vf[n][i].Reshape(2, 2, false);
      vf[n][i](2, 0) = 4.0;
      vf[n][i](2, 1) = 5.0;
      vf[n][i](2, 2) = 6.0;
    }
  }
  projector.HarmonicCell_CRk(cell, 2, vf, moments, uc);
  std::cout << uc[0] << std::endl;

  int dir;
  double val1(0.0), valx(0.0);
  NumericalIntegration numi(mesh);

  for (int n = 0; n < nfaces; ++n) {
    const AmanziGeometry::Point& normal = mesh->face_normal(n, false, cell, &dir);
    double factor = normal[0] / mesh->face_area(n) * dir;

    std::vector<const Polynomial*> polys;

    Polynomial tmp = vf[n][0] - uc[0];
    polys.push_back(&tmp);
    val1 += factor * numi.IntegratePolynomialsFace(n, polys);

    Polynomial q(2, 1);
    q(1, 0) = 1.0;
    q(1, 1) = 2.0;
    polys.push_back(&q);
    valx += factor * numi.IntegratePolynomialsFace(n, polys);
  }
  std::cout << "values: " << val1 << " " << valx << std::endl;
  CHECK_CLOSE(val1, 0.0, 1e-12);
  CHECK_CLOSE(valx, 0.0, 1e-12);

  // preservation of moments (reusing previous boundary functions)
  std::cout << "\nTest: HO nonconforming projectors for pentagon (verify moments)" << std::endl;
  for (int k = 2; k < 4; ++k) {
    Polynomial vc(2, k - 2);
    int nk = vc.size();

    moments->Reshape(nk);
    for (int i = 0; i < nk; ++i) {
      (*moments)(i) = 1.0 + i;
    }

    projector.EllipticCell_CRk(cell, k, vf, moments, uc);

    for (auto it = vc.begin(); it.end() <= vc.end(); ++it) {
      Polynomial mono(2, it.multi_index(), 1.0);
      mono.set_origin(mesh->cell_centroid(cell));
   
      Polynomial poly(uc[0]);
      poly.ChangeOrigin(mesh->cell_centroid(cell));
      poly *= mono;

      double val = numi.IntegratePolynomialCell(cell, poly) / mesh->cell_volume(cell);
      int n = it.PolynomialPosition();
      if (n == 0) CHECK_CLOSE(1.0, val, 1e-12);
      if (n >= 1) CHECK(fabs(val - (1.0 + n)) > 0.05);
    }
  }

  delete comm;
}


/* **************************************************************** */
TEST(HARMONIC_PROJECTORS_SQUARE_PK) {
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::WhetStone;

  std::cout << "\nTest: HO Lagrange projectors for square (linear deformation)" << std::endl;
  Epetra_MpiComm *comm = new Epetra_MpiComm(MPI_COMM_WORLD);

  MeshFactory meshfactory(comm);
  meshfactory.preference(FrameworkPreference({MSTK}));
  Teuchos::RCP<const Amanzi::AmanziGeometry::GeometricModel> gm;
  Teuchos::RCP<Mesh> mesh = meshfactory(0.0, 0.0, 1.2, 1.1, 1, 1, gm, true, true); 
 
  int cell(0);
  VectorPolynomial uc, uc2;
  std::vector<VectorPolynomial> vf(4);

  // test zero cell deformation
  for (int n = 0; n < 4; ++n) {
    vf[n].resize(2);
    for (int i = 0; i < 2; ++i) {
      vf[n][i].Reshape(2, 1, true);
    }
  }

  Projector projector(mesh);
  auto moments = std::make_shared<WhetStone::DenseVector>();

  // test linear deformation
  for (int n = 0; n < 4; ++n) {
    for (int i = 0; i < 2; ++i) {
      vf[n][i](0, 0) = 1.0;
      vf[n][i](1, 0) = 2.0;
      vf[n][i](1, 1) = 3.0;
    }
  }
  
  for (int k = 1; k < 4; ++k) {
    projector.HarmonicCell_Pk(cell, k, vf, moments, uc);  
    for (int i = 0; i < 2; ++i) {
      uc[i] -= vf[0][i];
      CHECK(uc[i].NormMax() < 1e-12);
    }
  }

  // test re-location of the right-top corner to (2,3)
  // cross-check with the CR projectors
  std::cout << "Test: HO Lagrange projectors for square (bilinear deformation)" << std::endl;
  for (int n = 0; n < 4; ++n) vf[n].PutScalar(0.0);
  vf[1][0](1, 1) = 0.8 / 1.1; 
  vf[1][1](1, 1) = 1.9 / 1.1; 

  vf[2][0](1, 0) = 0.8 / 1.2; 
  vf[2][1](1, 0) = 1.9 / 1.2; 

  for (int k = 1; k < 3; ++k) { 
    projector.HarmonicCell_Pk(cell, k, vf, moments, uc);  
    projector.HarmonicCell_CRk(cell, k, vf, moments, uc2);
    for (int i = 0; i < 2; ++i) {
      uc2[i] -= uc[i];
      CHECK(uc2[i].NormMax() < 1e-12);
    }

    // Compare H1 and L2 projectors
    projector.L2Cell_Pk(cell, k, vf, moments, uc2);
    uc2[0] -= uc[0];
    CHECK(uc2[0].NormMax() < 1e-12);
  }

  auto p = AmanziGeometry::Point(1.2, 1.1);
  CHECK(fabs(uc[0].Value(p) - 0.8) < 1e-12 &&
        fabs(uc[1].Value(p) - 1.9) < 1e-12);

  delete comm;
}

/* **************************************************************** */
TEST(HARMONIC_PROJECTORS_POLYGON_PK) {
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::WhetStone;

  std::cout << "\nTest: HO Lagrange projectors for pentagon (linear deformation)" << std::endl;
  Epetra_MpiComm *comm = new Epetra_MpiComm(MPI_COMM_WORLD);

  MeshFactory meshfactory(comm);
  meshfactory.preference(FrameworkPreference({MSTK}));
  Teuchos::RCP<const Amanzi::AmanziGeometry::GeometricModel> gm;
  Teuchos::RCP<Mesh> mesh = meshfactory("test/one_pentagon.exo", gm, true, true);
  // Teuchos::RCP<Mesh> mesh = meshfactory("test/one_quad.exo", gm, true, true);
 
  int cell(0), nfaces(5);
  VectorPolynomial uc, uc2;
  std::vector<VectorPolynomial> vf(nfaces);

  Projector projector(mesh);
  auto moments = std::make_shared<WhetStone::DenseVector>();

  // test globally linear deformation
  for (int n = 0; n < nfaces; ++n) {
    vf[n].resize(2);
    for (int i = 0; i < 2; ++i) {
      vf[n][i].Reshape(2, 1, true);
      vf[n][i](0, 0) = 1.0;
      vf[n][i](1, 0) = 2.0;
      vf[n][i](1, 1) = 3.0;
    }
  }
  
  for (int k = 1; k < 4; ++k) {
    projector.HarmonicCell_Pk(cell, k, vf, moments, uc);
    std::cout << uc[0] << std::endl;

    uc[0] -= vf[0][0];
    uc[1] -= vf[0][1];
    CHECK(uc[0].NormMax() < 1e-12 && uc[1].NormMax() < 1e-12);
  }

  // test globally quadratic deformation
  std::cout << "\nTest: HO Lagrange for pentagon (quadratic deformation)" << std::endl;
  for (int n = 0; n < nfaces; ++n) {
    for (int i = 0; i < 2; ++i) {
      vf[n][i].Reshape(2, 2, false);
      vf[n][i](2, 0) = 4.0;
      vf[n][i](2, 1) = 5.0;
      vf[n][i](2, 2) = -4.0;
    }
  }

  for (int k = 2; k < 4; ++k) {
    projector.HarmonicCell_Pk(cell, k, vf, moments, uc);
    std::cout << uc[0] << std::endl;
    uc[0] -= vf[0][0];
    uc[1] -= vf[0][1];
    CHECK(uc[0].NormMax() < 1e-10 && uc[1].NormMax() < 1e-10);
  }

  // test trace compatibility between function and its projecton (k < 3 only!)
  std::cout << "\nTest: HO Lagrange projectors for pentagon (harmonic function)" << std::endl;
  for (int n = 0; n < nfaces; ++n) {
    for (int i = 0; i < 2; ++i) {
      vf[n][i].Reshape(2, 2, false);
      vf[n][i](2, 0) = 4.0;
      vf[n][i](2, 1) = 5.0;
      vf[n][i](2, 2) = 6.0;
    }
  }
  projector.HarmonicCell_Pk(cell, 2, vf, moments, uc);
  std::cout << uc[0] << std::endl;

  int dir;
  double val1(0.0), valx(0.0);
  NumericalIntegration numi(mesh);

  for (int n = 0; n < nfaces; ++n) {
    const AmanziGeometry::Point& normal = mesh->face_normal(n, false, cell, &dir);
    double factor = normal[0] / mesh->face_area(n) * dir;

    std::vector<const Polynomial*> polys;

    Polynomial tmp = vf[n][0] - uc[0];
    polys.push_back(&tmp);
    val1 += factor * numi.IntegratePolynomialsFace(n, polys);

    Polynomial q(2, 1);
    q(1, 0) = 1.0;
    q(1, 1) = 2.0;
    polys.push_back(&q);
    valx += factor * numi.IntegratePolynomialsFace(n, polys);
  }
  std::cout << "values: " << val1 << " " << valx << std::endl;
  CHECK_CLOSE(val1, 0.0, 1e-12);
  CHECK_CLOSE(valx, 0.0, 1e-12);

  // test piecewise linear deformation
  std::cout << "\nTest: HO Lagrange projectors for pentagon (piece-wice linear deformation)" << std::endl;
  std::vector<AmanziGeometry::Point> vv;
  vv.push_back(AmanziGeometry::Point( 0.0, 0.0));
  vv.push_back(AmanziGeometry::Point( 0.0,-0.1));
  vv.push_back(AmanziGeometry::Point( 0.1, 0.0));
  vv.push_back(AmanziGeometry::Point( 0.0, 0.1));
  vv.push_back(AmanziGeometry::Point(-0.1, 0.0));

  AmanziGeometry::Point x1(2), x2(2), tau(2);
  for (int n = 0; n < 5; ++n) {
    int m = (n + 1) % 5;
    mesh->node_get_coordinates(n, &x1);
    mesh->node_get_coordinates(m, &x2);
    tau = x2 - x1;
    tau /= AmanziGeometry::L22(tau);

    for (int i = 0; i < 2; ++i) {
      vf[n][i].Reshape(2, 1);
      vf[n][i](0, 0) = vv[n][i] * (x2 * tau) - vv[m][i] * (x1 * tau);

      vf[n][i](1, 0) = (vv[m][i] - vv[n][i]) * tau[0];
      vf[n][i](1, 1) = (vv[m][i] - vv[n][i]) * tau[1];
    }
  }
  
  for (int k = 1; k < 4; ++k) {
    projector.HarmonicCell_Pk(cell, k, vf, moments, uc);
    projector.HarmonicCell_CRk(cell, k, vf, moments, uc2);
    uc2 -= uc;
    if (k == 1) CHECK(uc2[0].NormMax() < 1e-12);
    if (k == 2) CHECK(uc2[0].NormMax() > 1e-3 && uc2[0].NormMax() < 2e-3);

    projector.L2Cell_Pk(cell, k, vf, moments, uc2);
    uc2 -= uc;
    if (k < 3) CHECK(uc2[0].NormMax() < 1e-12);
    if (k > 1) std::cout << "k=" << k << " moments:" << *moments;
  }

  // preservation of moments (reusing previous boundary functions)
  std::cout << "\nTest: HO Lagrange projectors for pentagon (verify moments)" << std::endl;
  for (int k = 2; k < 4; ++k) {
    moments->Reshape(k * (k - 1) / 2);
    moments->PutScalar(1.0);
    projector.EllipticCell_Pk(cell, k, vf, moments, uc);
    double tmp = numi.IntegratePolynomialCell(cell, uc[0]) / mesh->cell_volume(cell);
    CHECK_CLOSE(1.0, tmp, 1e-12);
  }

  delete comm;
}

/* **************************************************************** */
TEST(SERENDIPITY_PROJECTORS_POLYGON_PK) {
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::WhetStone;

  std::cout << "\nTest: HO Serendipity Lagrange projectors for pentagon (linear deformation)" << std::endl;
  Epetra_MpiComm *comm = new Epetra_MpiComm(MPI_COMM_WORLD);

  MeshFactory meshfactory(comm);
  meshfactory.preference(FrameworkPreference({MSTK}));
  Teuchos::RCP<const Amanzi::AmanziGeometry::GeometricModel> gm;
  Teuchos::RCP<Mesh> mesh = meshfactory("test/one_pentagon.exo", gm, true, true);
 
  int cell(0), nfaces(5);
  VectorPolynomial uc, uc2;
  std::vector<VectorPolynomial> vf(nfaces);

  Projector projector(mesh);
  auto moments = std::make_shared<WhetStone::DenseVector>();

  // test globally linear deformation
  for (int n = 0; n < nfaces; ++n) {
    vf[n].resize(2);
    for (int i = 0; i < 2; ++i) {
      vf[n][i].Reshape(2, 1, true);
      vf[n][i](0, 0) = 1.0;
      vf[n][i](1, 0) = 2.0;
      vf[n][i](1, 1) = 3.0;
    }
  }
  
  for (int k = 1; k < 4; ++k) {
    projector.L2Cell_SerendipityPk(cell, k, vf, moments, uc);
    std::cout << uc[0] << std::endl;

    uc[0] -= vf[0][0];
    uc[1] -= vf[0][1];
    CHECK(uc[0].NormMax() < 1e-12 && uc[1].NormMax() < 1e-12);
  }

  // test globally quadratic deformation
  std::cout << "\nTest: HO Serendipity Lagrange projectors for pentagon (quadratic deformation)" << std::endl;
  for (int n = 0; n < nfaces; ++n) {
    for (int i = 0; i < 2; ++i) {
      vf[n][i].Reshape(2, 2, false);
      vf[n][i](2, 0) = 4.0;
      vf[n][i](2, 1) = 5.0;
      vf[n][i](2, 2) = 6.0;
    }
  }

  for (int k = 2; k < 4; ++k) {
    projector.L2Cell_SerendipityPk(cell, k, vf, moments, uc);
    // std::cout << uc[0] << std::endl;
    uc[0] -= vf[0][0];
    uc[1] -= vf[0][1];
    CHECK(uc[0].NormMax() < 1e-10 && uc[1].NormMax() < 1e-10);
  }

  // test piecewise linear deformation
  // work in progress...
  std::cout << "\nTest: HO Sependipity Lagrange projectors for pentagon (piece-wice linear)" << std::endl;
  std::vector<AmanziGeometry::Point> vv;
  vv.push_back(AmanziGeometry::Point( 0.0, 0.0));
  vv.push_back(AmanziGeometry::Point( 0.05,-0.1));
  vv.push_back(AmanziGeometry::Point( 0.1, 0.0));
  vv.push_back(AmanziGeometry::Point( 0.0, 0.1));
  vv.push_back(AmanziGeometry::Point(-0.1, 0.0));

  AmanziGeometry::Point x1(2), x2(2), tau(2);
  for (int n = 0; n < 5; ++n) {
    int m = (n + 1) % 5;
    mesh->node_get_coordinates(n, &x1);
    mesh->node_get_coordinates(m, &x2);
    tau = x2 - x1;
    tau /= AmanziGeometry::L22(tau);

    for (int i = 0; i < 2; ++i) {
      vf[n][i].Reshape(2, 1);
      vf[n][i](0, 0) = vv[n][i] * (x2 * tau) - vv[m][i] * (x1 * tau);

      vf[n][i](1, 0) = (vv[m][i] - vv[n][i]) * tau[0];
      vf[n][i](1, 1) = (vv[m][i] - vv[n][i]) * tau[1];
    }
  }
  
  for (int k = 1; k < 4; ++k) {
    projector.L2Cell_SerendipityPk(cell, k, vf, moments, uc);
    uc[0].ChangeOrigin(mesh->cell_centroid(cell));
    std::cout << uc[0] << std::endl;
  }

  delete comm;
}

