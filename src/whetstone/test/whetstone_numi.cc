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

#include "NumericalIntegration.hh"
#include "Polynomial.hh"


/* **************************************************************** */
TEST(NUMI_CELL_2D_EULER_FORMULA) {
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::WhetStone;

  std::cout << "Test: Numerical integration: Euler's formula for polygon" << std::endl;
  Epetra_MpiComm *comm = new Epetra_MpiComm(MPI_COMM_WORLD);

  MeshFactory meshfactory(comm);
  meshfactory.preference(FrameworkPreference({MSTK}));
  Teuchos::RCP<const Amanzi::AmanziGeometry::GeometricModel> gm;
  Teuchos::RCP<Mesh> mesh = meshfactory("test/one_pentagon.exo", gm, true, true); 
 
  NumericalIntegration numi(mesh);

  int cell(0);
  double val;

  // 0th-order polynomial
  Polynomial poly(2, 0);
  poly(0, 0) = 1.0;
  val = numi.IntegratePolynomialCell(cell, poly);

  printf("order=0  value=%10.6g\n", val);
  CHECK_CLOSE(val, mesh->cell_volume(cell), 1e-10);
 
  // 1st-order polynomial
  poly.Reshape(2, 1);
  poly(1, 0) = 2.0;
  poly(1, 1) = 3.0;
  poly.set_origin(mesh->cell_centroid(cell));
  val = numi.IntegratePolynomialCell(cell, poly);

  printf("order=1  value=%10.6g\n", val);
  CHECK_CLOSE(val, mesh->cell_volume(cell), 1e-10);
 
  delete comm;
}


/* **************************************************************** */
TEST(NUMI_CELL_2D_QUADRATURE_POLYGON) {
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::WhetStone;

  std::cout << "Test: Numerical integration: quadrature rule" << std::endl;
  Epetra_MpiComm *comm = new Epetra_MpiComm(MPI_COMM_WORLD);

  MeshFactory meshfactory(comm);
  meshfactory.preference(FrameworkPreference({MSTK}));
  Teuchos::RCP<const Amanzi::AmanziGeometry::GeometricModel> gm;
  Teuchos::RCP<Mesh> mesh = meshfactory("test/one_pentagon.exo", gm, true, true); 
 
  NumericalIntegration numi(mesh);

  int cell(0);
  double val1, val2;

  // 1st-order polynomial
  Polynomial poly(2, 1);
  poly(1, 0) = 2.0;
  poly(1, 1) = 3.0;
  poly.set_origin(mesh->cell_centroid(cell));

  std::vector<const WhetStoneFunction*> polys(1);
  polys[0] = &poly;

  for (int order = 1; order < 10; ++order) {
    val1 = numi.IntegrateFunctionsTrianglatedCell(cell, polys, order);

    printf("order=%d  value=%10.6g\n", order, val1);
    CHECK_CLOSE(val1, poly.Value(mesh->cell_centroid(cell)), 1e-12);
  }
 
  // cross-comparison of integartors
  for (int order = 1; order < 10; ++order) {
    poly.Reshape(2, order, true);
    poly(order, 0) = 1.0;

    val1 = numi.IntegrateFunctionsTrianglatedCell(cell, polys, order);
    val2 = numi.IntegratePolynomialCell(cell, poly);

    printf("order=%d  values: %10.6g %10.6g\n", order, val1, val2);
    CHECK_CLOSE(val1, val2, 1e-12);
  }

  delete comm;
}


/* **************************************************************** */
TEST(NUMI_CELL_2D_QUADRATURE_SQUARE) {
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::WhetStone;

  std::cout << "Test: Numerical integration: quadrature rule for square" << std::endl;
  Epetra_MpiComm *comm = new Epetra_MpiComm(MPI_COMM_WORLD);

  MeshFactory meshfactory(comm);
  meshfactory.preference(FrameworkPreference({MSTK}));
  Teuchos::RCP<const Amanzi::AmanziGeometry::GeometricModel> gm;
  Teuchos::RCP<Mesh> mesh = meshfactory(-1.0, -1.0, 1.0, 1.0, 2, 2); 
 
  NumericalIntegration numi(mesh);

  int cell(3);
  double val, exact;

  for (int order = 1; order < 10; ++order) {
    Polynomial poly(2, order);
    poly(order, 0) = 1.0;

    std::vector<const WhetStoneFunction*> polys(1);
    polys[0] = &poly;

    val = numi.IntegrateFunctionsTrianglatedCell(cell, polys, poly.order());
    exact = 1.0 / (order + 1);

    printf("order=%d  value=%10.6g  exact=%10.6g\n", order, val, exact);
    CHECK_CLOSE(val, exact, 1e-12);
  }
 
  delete comm;
}

