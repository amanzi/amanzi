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

#include "Teuchos_RCP.hpp"
#include "UnitTest++.h"

#include "Mesh.hh"
#include "MeshFactory.hh"

#include "DenseMatrix.hh"
#include "dg.hh"


/* ****************************************************************
* Test of DG mass matrices
**************************************************************** */
TEST(DG_MASS_MATRIX) {
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::WhetStone;

  std::cout << "Test: DG mass matrices" << std::endl;
  Epetra_MpiComm *comm = new Epetra_MpiComm(MPI_COMM_WORLD);

  MeshFactory meshfactory(comm);
  meshfactory.preference(FrameworkPreference({MSTK}));
  Teuchos::RCP<Mesh> mesh = meshfactory(0.0, 0.0, 0.5, 0.5, 1, 1); 
 
  DG dg(mesh);

  for (int k = 0; k < 3; k++) {
    int nk = (k + 1) * (k + 2) / 2;
    DenseMatrix M(nk, nk);

    dg.TaylorMassMatrix(0, k, 1.0, M);

    printf("Mass matrix for order=%d\n", k);
    for (int i = 0; i < nk; i++) {
      for (int j = 0; j < nk; j++ ) printf("%9.5f ", M(i, j)); 
      printf("\n");
    }

    double area = mesh->cell_volume(0);
    CHECK_CLOSE(M(0, 0), area, 1e-12);
    if (k > 0) {
      CHECK_CLOSE(M(1, 1), area / 48, 1e-12);
    }
    if (k > 1) {
      CHECK_CLOSE(M(3, 3), area / 1280, 1e-12);
      // CHECK_CLOSE(M(4, 4), area / 144, 1e-12);
    }
  }

  delete comm;
}


/* ****************************************************************
* Test of DG advection matrices in a cell
**************************************************************** */
TEST(DG_ADVECTION_MATRIX_CELL) {
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::WhetStone;

  std::cout << "\nTest: DG advection matrices in cells" << std::endl;
  Epetra_MpiComm *comm = new Epetra_MpiComm(MPI_COMM_WORLD);

  MeshFactory meshfactory(comm);
  meshfactory.preference(FrameworkPreference({MSTK}));
  // Teuchos::RCP<Mesh> mesh = meshfactory(0.0, 0.0, 1.0, 1.0, 1, 1); 
  Teuchos::RCP<Mesh> mesh = meshfactory("test/one_cell2.exo");
 
  DG dg(mesh);

  AmanziGeometry::Point zero(0.0, 0.0);

  for (int k = 0; k < 3; k++) {
    int nk = (k + 1) * (k + 2) / 2;
    DenseMatrix A0(nk, nk), A1(nk, nk);

    std::vector<AmanziGeometry::Point> u;
    u.push_back(AmanziGeometry::Point(1.0, 2.0));

    dg.TaylorAdvectionMatrixCell(0, k, u, A0);

    printf("Advection matrix (cell-based) for order=%d\n", k);
    for (int i = 0; i < nk; i++) {
      for (int j = 0; j < nk; j++ ) printf("%8.4f ", A0(i, j)); 
      printf("\n");
    }

    // TEST1: accuracy
    DenseVector v(nk), av(nk);
    if (k > 1) {
      double tmp;
      const AmanziGeometry::Point& xc = mesh->cell_centroid(0);

      v.PutScalar(0.0);
      v(0) = xc[0] + 2 * xc[1];
      v(1) = 1.0;
      v(2) = 2.0;
    
      A0.Multiply(v, av, false);
      v.Dot(av, &tmp);
      CHECK_CLOSE(tmp, 5 * v(0) * mesh->cell_volume(0), 1e-12);
    }

    // TEST2: add zero linear component to constant u
    u.push_back(zero);
    u.push_back(zero);

    dg.TaylorAdvectionMatrixCell(0, k, u, A1);

    A1 -= A0;
    CHECK_CLOSE(0.0, A1.NormInf(), 1e-12);

    // TEST3: nonzero linear component of u
    u.clear();
    u.push_back(zero);
    u.push_back(zero);
    u.push_back(AmanziGeometry::Point(0.0, 1.0));

    dg.TaylorAdvectionMatrixCell(0, k, u, A1);

    printf("Advection matrix (cell-based) for order=%d u=(0,y)\n", k);
    for (int i = 0; i < nk; i++) {
      for (int j = 0; j < nk; j++ ) printf("%8.4f ", A1(i, j)); 
      printf("\n");
    }

    if (k > 1) {
      double mon_y2(0.384318693694);
      CHECK_CLOSE(2 * mon_y2, A1(0, 5), 1e-12);
      CHECK_CLOSE(mon_y2, A1(2, 2), 1e-12);
    }
  }

  delete comm;
}


/* ****************************************************************
* Test of DG advection matrices on a face
**************************************************************** */
TEST(DG_ADVECTION_MATRIX_FACE) {
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::WhetStone;

  std::cout << "\nTest: DG advection matrices on faces" << std::endl;
  Epetra_MpiComm *comm = new Epetra_MpiComm(MPI_COMM_WORLD);

  MeshFactory meshfactory(comm);
  meshfactory.preference(FrameworkPreference({MSTK}));
  Teuchos::RCP<Mesh> mesh = meshfactory(0.0, 0.0, 1.0, 1.0, 2, 2); 
 
  DG dg(mesh);

  AmanziGeometry::Point zero(0.0, 0.0);

  for (int k = 0; k < 2; k++) {
    int nk = (k + 1) * (k + 2);
    DenseMatrix A0(nk, nk), A1(nk, nk);

    std::vector<AmanziGeometry::Point> u;
    u.push_back(AmanziGeometry::Point(1.0, 2.0));

    // TEST1: constant u
    dg.TaylorAdvectionMatrixFace(1, k, u, A0);

    printf("Advection matrix (face-based) for order=%d  u=constant\n", k);
    for (int i = 0; i < nk; i++) {
      for (int j = 0; j < nk; j++ ) printf("%8.4f ", A0(i, j)); 
      printf("\n");
    }

    // TEST2: linear u with zero gradient
    u.push_back(zero);
    u.push_back(zero);

    dg.TaylorAdvectionMatrixFace(1, k, u, A1);

    A1 -= A0;
    CHECK_CLOSE(0.0, A1.NormInf(), 1e-12);

    // TEST3: nonzero linear component of u
    u.clear();
    u.push_back(zero);
    u.push_back(zero);
    u.push_back(AmanziGeometry::Point(1.0, 0.0));

    dg.TaylorAdvectionMatrixFace(1, k, u, A1);

    printf("Advection matrix (cell-based) for order=%d u=(y-y0,0)\n", k);
    for (int i = 0; i < nk; i++) {
      for (int j = 0; j < nk; j++ ) printf("%8.4f ", A1(i, j)); 
      printf("\n");
    }
  }

  delete comm;
}


/* ****************************************************************
* Test of polynomial approximation in cells
**************************************************************** */
TEST(DG_MAP_APPROXIMATION_CELL) {
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::WhetStone;

  std::cout << "\nTest: Polynomial approximation of map in cells." << std::endl;
  Epetra_MpiComm *comm = new Epetra_MpiComm(MPI_COMM_WORLD);

  MeshFactory meshfactory(comm);
  meshfactory.preference(FrameworkPreference({MSTK}));
  Teuchos::RCP<Mesh> mesh = meshfactory("test/one_cell2.exo");

  // extract polygon from the mesh
  Entity_ID_List nodes;
  AmanziGeometry::Point xv;
  std::vector<AmanziGeometry::Point> x1;

  mesh->cell_get_nodes(0, &nodes);

  for (int i = 0; i < nodes.size(); ++i) {
    mesh->node_get_coordinates(nodes[i], &xv);
    x1.push_back(xv);
  }

  // test identity map
  DG dg(mesh);
  std::vector<AmanziGeometry::Point> u;
  AmanziGeometry::Point ex(1.0, 0.0), ey(0.0, 1.0);

  dg.TaylorLeastSquareFit(1, x1, x1, u);
  CHECK_CLOSE(norm(u[0]), 0.0, 1e-12);
  CHECK_CLOSE(norm(u[1] - ex), 0.0, 1e-12);
  CHECK_CLOSE(norm(u[2] - ey), 0.0, 1e-12);

  // test linear map
  std::vector<AmanziGeometry::Point> x2(x1);
  AmanziGeometry::Point shift(0.1, 0.2);
  for (int i = 0; i < nodes.size(); ++i) {
    x2[i] += shift;
  }

  dg.TaylorLeastSquareFit(1, x1, x2, u);
  CHECK_CLOSE(norm(u[0] - shift), 0.0, 1e-12);
  CHECK_CLOSE(norm(u[1] - ex), 0.0, 1e-12);
  CHECK_CLOSE(norm(u[2] - ey), 0.0, 1e-12);

  // test rotation map
  double s(std::sin(0.3)), c(std::cos(0.3));
  for (int i = 0; i < nodes.size(); ++i) {
    x2[i][0] = c * x1[i][0] - s * x1[i][1];
    x2[i][1] = s * x1[i][0] + c * x1[i][1];
  }

  dg.TaylorLeastSquareFit(1, x1, x2, u);
  CHECK_CLOSE(norm(u[0]), 0.0, 1e-12);
  CHECK_CLOSE(norm(u[1] - AmanziGeometry::Point(c, s)), 0.0, 1e-12);
  CHECK_CLOSE(norm(u[2] - AmanziGeometry::Point(-s, c)), 0.0, 1e-12);

  // test non-linear deformation map
  x1.clear();
  x1.push_back(AmanziGeometry::Point(-0.5, -0.5));
  x1.push_back(AmanziGeometry::Point( 0.5, -0.5));
  x1.push_back(AmanziGeometry::Point(-0.5,  0.5));
  x1.push_back(AmanziGeometry::Point( 0.5,  0.5));

  x2 = x1;
  x2[3] += AmanziGeometry::Point(0.1, 0.1);

  dg.TaylorLeastSquareFit(1, x1, x2, u);

  for (int i = 0; i < u.size(); ++i) {
    printf("u[%d] = %8.4g %8.4g\n", i, u[i][0], u[i][1]); 
  }

  CHECK_CLOSE(0.025, u[0][0], 1e-12);

  delete comm;
}

