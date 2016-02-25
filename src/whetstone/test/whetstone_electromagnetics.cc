/*
  This is the mimetic discretization component of the Amanzi code. 

  Copyright 2010-2012 held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <vector>

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_LAPACK.hpp"
#include "UnitTest++.h"

#include "MeshFactory.hh"
#include "Mesh.hh"
#include "Point.hh"

#include "mfd3d_electromagnetics.hh"
#include "Tensor.hh"


/* **************************************************************** */
TEST(MASS_MATRIX_2D) {
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::WhetStone;

  std::cout << "\nTest: Mass matrix for edge elements in 2D" << std::endl;
#ifdef HAVE_MPI
  Epetra_MpiComm *comm = new Epetra_MpiComm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm *comm = new Epetra_SerialComm();
#endif

  FrameworkPreference pref;
  pref.clear();
  pref.push_back(MSTK);

  MeshFactory factory(comm);
  factory.preference(pref);

  bool request_faces(true), request_edges(true);
  // Teuchos::RCP<Mesh> mesh = factory(0.0, 0.0, 1.0, 1.0, 1, 1,NULL, true, true); 
  Teuchos::RCP<Mesh> mesh = factory("test/two_cell2.exo", Teuchos::null, request_faces, request_edges); 
 
  MFD3D_Electromagnetics mfd(mesh);

  int cell = 0;
  AmanziMesh::Entity_ID_List edges;
  mesh->cell_get_edges(cell, &edges);

  int nedges = edges.size();
  int nrows = nedges;

  Tensor T(2, 2);
  T(0, 0) = 2.0;
  T(1, 1) = 1.0;
  T(0, 1) = 1.0;
  T(1, 0) = 1.0;

  for (int method = 0; method < 4; method++) {
    DenseMatrix M(nrows, nrows);

    if (method == 0) {
      mfd.MassMatrix(cell, T, M);
    } else if (method == 1) {
      mfd.MassMatrixOptimized(cell, T, M);
    } else if (method == 2) {
      mfd.MassMatrixInverse(cell, T, M);
      M.Inverse();
    } else if (method == 3) {
      mfd.MassMatrixInverseOptimized(cell, T, M);
      M.Inverse();
    }

    printf("Mass matrix for cell %3d method=%d\n", cell, method);
    for (int i = 0; i < nrows; i++) {
      for (int j = 0; j < nrows; j++ ) printf("%8.4f ", M(i, j)); 
      printf("\n");
    }

    // verify SPD propery
    for (int i = 0; i < nrows; i++) CHECK(M(i, i) > 0.0);

    // verify exact integration property
    double xi, yi, xj;
    double vxx = 0.0, vxy = 0.0, volume = mesh->cell_volume(cell); 
    for (int i = 0; i < nedges; i++) {
      int e1 = edges[i];
      const AmanziGeometry::Point& t1 = mesh->edge_vector(e1);
      double a1 = mesh->edge_length(e1);

      xi = t1[0] / a1;
      yi = t1[1] / a1;
      for (int j = 0; j < nedges; j++) {
        int e2 = edges[j];
        const AmanziGeometry::Point& t2 = mesh->edge_vector(e2);
        double a2 = mesh->edge_length(e2);
        xj = t2[0] / a2;
        vxx += M(i, j) * xi * xj;
        vxy += M(i, j) * yi * xj;
      }
    }
    CHECK_CLOSE(volume, vxx, 1e-10);
    CHECK_CLOSE(-volume, vxy, 1e-10);
  }

  delete comm;
}


/* **************************************************************** */
TEST(MASS_MATRIX_3D) {
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::WhetStone;

  std::cout << "\nTest: Mass matrix for edge elements in 3D" << std::endl;
#ifdef HAVE_MPI
  Epetra_MpiComm *comm = new Epetra_MpiComm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm *comm = new Epetra_SerialComm();
#endif

  FrameworkPreference pref;
  pref.clear();
  pref.push_back(MSTK);

  MeshFactory meshfactory(comm);
  meshfactory.preference(pref);

  bool request_faces(true), request_edges(true);

  // RCP<Mesh> mesh = meshfactory("test/dodecahedron.exo", NULL, request_faces, request_edges); 
  RCP<Mesh> mesh = meshfactory("test/one_cell.exo", Teuchos::null, request_faces, request_edges); 
  // RCP<Mesh> mesh = meshfactory(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1, 2, 3, Teuchos::null, true, true); 
 
  MFD3D_Electromagnetics mfd(mesh);

  int cell = 0;
  AmanziMesh::Entity_ID_List edges;
  mesh->cell_get_edges(cell, &edges);

  int nedges = edges.size();
  int nrows = nedges;

  Tensor T(3, 2);
  T(0, 0) = 2.0;
  T(1, 1) = 1.0;
  T(0, 1) = 1.0;
  T(1, 0) = 1.0;
  T(2, 2) = 1.0;

  for (int method = 0; method < 4; method++) {
    DenseMatrix M(nrows, nrows);

    if (method == 0) {
      mfd.MassMatrix(cell, T, M);
    } else if (method == 1) {
      mfd.MassMatrixOptimized(cell, T, M);
    } else if (method == 2) {
      mfd.MassMatrixInverse(cell, T, M);
      M.Inverse();
    } else if (method == 3) {
      mfd.MassMatrixInverseOptimized(cell, T, M);
      M.Inverse();
    }

    printf("Mass matrix for cell %3d method=%d\n", cell, method);
    for (int i = 0; i < nrows; i++) {
      for (int j = 0; j < nrows; j++ ) printf("%8.4f ", M(i, j)); 
      printf("\n");
    }

    // verify SPD propery
    for (int i = 0; i < nrows; i++) CHECK(M(i, i) > 0.0);

    // verify exact integration property
    double xi, yi, xj;
    double vxx = 0.0, vxy = 0.0, volume = mesh->cell_volume(cell); 
    for (int i = 0; i < nedges; i++) {
      int e1 = edges[i];
      const AmanziGeometry::Point& t1 = mesh->edge_vector(e1);
      double a1 = mesh->edge_length(e1);

      xi = t1[0] / a1;
      yi = t1[1] / a1;
      for (int j = 0; j < nedges; j++) {
        int e2 = edges[j];
        const AmanziGeometry::Point& t2 = mesh->edge_vector(e2);
        double a2 = mesh->edge_length(e2);
        xj = t2[0] / a2;
        vxx += M(i, j) * xi * xj;
        vxy += M(i, j) * yi * xj;
      }
    }
    CHECK_CLOSE(volume, vxx, 1e-10);
    CHECK_CLOSE(-volume, vxy, 1e-10);
  }

  delete comm;
}


/* **************************************************************** */
TEST(STIFFNESS_MATRIX_2D) {
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::WhetStone;

  std::cout << "\nTest: Stiffness matrix for edge elements in 2D" << std::endl;
#ifdef HAVE_MPI
  Epetra_MpiComm *comm = new Epetra_MpiComm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm *comm = new Epetra_SerialComm();
#endif

  FrameworkPreference pref;
  pref.clear();
  pref.push_back(MSTK);

  MeshFactory factory(comm);
  factory.preference(pref);

  bool request_faces(true), request_edges(true);
  // Teuchos::RCP<Mesh> mesh = factory(0.0, 0.0, 1.0, 1.0, 1, 1, Teuchos::null, true, true); 
  Teuchos::RCP<Mesh> mesh = factory("test/two_cell2.exo", Teuchos::null, request_faces, request_edges); 
 
  MFD3D_Electromagnetics mfd(mesh);

  int cell = 0;
  AmanziMesh::Entity_ID_List faces;
  std::vector<int> dirs;
  mesh->cell_get_faces_and_dirs(cell, &faces, &dirs);

  int nfaces = faces.size();
  int nrows = nfaces;

  Tensor T(2, 1);
  T(0, 0) = 1.0;

  for (int method = 0; method < 2; method++) {
    DenseMatrix A(nrows, nrows);

    if (method == 0) {
      mfd.StiffnessMatrix(cell, T, A);
    } else if (method == 1) {
      mfd.StiffnessMatrixOptimized(cell, T, A);
    }

    printf("Stiffness matrix for cell %3d method=%d\n", cell, method);
    for (int i = 0; i < nrows; i++) {
      for (int j = 0; j < nrows; j++ ) printf("%8.4f ", A(i, j)); 
      printf("\n");
    }

    // verify SPD propery
    for (int i = 0; i < nrows; i++) CHECK(A(i, i) > 0.0);

    // verify exact integration property
    double xi, xj, vxx(0.0);
    AmanziGeometry::Point p1(2), p2(2);

    const AmanziGeometry::Point& xc = mesh->cell_centroid(cell);
    double volume = mesh->cell_volume(cell); 

    for (int i = 0; i < nrows; i++) {
      int f1 = faces[i];
      const AmanziGeometry::Point& n1 = mesh->face_normal(f1);
      const AmanziGeometry::Point& xf1 = mesh->face_centroid(f1);
      double a1 = mesh->face_area(f1);

      xi = ((xf1 - xc) * n1) * dirs[i] / a1;

      for (int j = 0; j < nrows; j++) {
        int f2 = faces[j];
        const AmanziGeometry::Point& n2 = mesh->face_normal(f2);
        const AmanziGeometry::Point& xf2 = mesh->face_centroid(f2);
        double a2 = mesh->face_area(f2);

        xj = ((xf2 - xc) * n2) * dirs[j] / a2;
        vxx += A(i, j) * xi * xj;
      }
    }
    CHECK_CLOSE(4 * volume, vxx, 1e-10);
  }

  delete comm;
}


/* **************************************************************** */
TEST(STIFFNESS_MATRIX_3D) {
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::WhetStone;

  std::cout << "\nTest: Stiffness matrix for edge elements in 3D" << std::endl;
#ifdef HAVE_MPI
  Epetra_MpiComm *comm = new Epetra_MpiComm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm *comm = new Epetra_SerialComm();
#endif

  FrameworkPreference pref;
  pref.clear();
  pref.push_back(MSTK);

  MeshFactory meshfactory(comm);
  meshfactory.preference(pref);

  bool request_faces(true), request_edges(true);

  // RCP<Mesh> mesh = meshfactory("test/dodecahedron.exo", Teuchos::null, request_faces, request_edges); 
  RCP<Mesh> mesh = meshfactory("test/one_cell.exo", Teuchos::null, request_faces, request_edges); 
  // RCP<Mesh> mesh = meshfactory(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1, 2, 3, Teuchos::null, true, true); 
 
  MFD3D_Electromagnetics mfd(mesh);

  int cell = 0;
  AmanziMesh::Entity_ID_List edges;
  mesh->cell_get_edges(cell, &edges);

  int nedges = edges.size();
  int nrows = nedges;

  Tensor T(3, 2);
  T(0, 0) = 2.0;
  T(1, 1) = 1.0;
  T(0, 1) = 1.0;
  T(1, 0) = 1.0;
  T(2, 2) = 1.0;

  for (int method = 0; method < 2; method++) {
    DenseMatrix A(nrows, nrows);

    if (method == 0) {
      mfd.StiffnessMatrix(cell, T, A);
    } else if (method == 1) {
      mfd.StiffnessMatrixOptimized(cell, T, A);
    }

    printf("Stiffness matrix for cell %3d\n", cell);
    for (int i = 0; i < nrows; i++) {
      for (int j = 0; j < nrows; j++ ) printf("%8.4f ", A(i, j)); 
      printf("\n");
    }

    // verify SPD propery
    for (int i = 0; i < nrows; i++) CHECK(A(i, i) > 0.0);

    // verify exact integration property
    int n1, n2;
    double xi, xj;
    double vxx = 0.0, volume = mesh->cell_volume(cell); 
    AmanziGeometry::Point p1(3), p2(3), v1(3);

    for (int i = 0; i < nedges; i++) {
      int e1 = edges[i];
      const AmanziGeometry::Point& t1 = mesh->edge_vector(e1);
      double a1 = mesh->edge_length(e1);

      mesh->edge_get_nodes(e1, &n1, &n2);
      mesh->node_get_coordinates(n1, &p1);
      mesh->node_get_coordinates(n2, &p2);
 
      v1 = (p1 + p2)^t1;
      xi = v1[0] / a1;

      for (int j = 0; j < nedges; j++) {
        int e2 = edges[j];
        const AmanziGeometry::Point& t2 = mesh->edge_vector(e2);
        double a2 = mesh->edge_length(e2);

        mesh->edge_get_nodes(e2, &n1, &n2);
        mesh->node_get_coordinates(n1, &p1);
        mesh->node_get_coordinates(n2, &p2);
 
        v1 = (p1 + p2)^t2;
        xj = v1[0] / a2;

        vxx += A(i, j) * xi * xj;
      }
    }
    CHECK_CLOSE(32 * volume, vxx, 1e-10);
  }

  delete comm;
}
