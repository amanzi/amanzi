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

#include "Teuchos_RCP.hpp"
#include "Teuchos_LAPACK.hpp"
#include "UnitTest++.h"

#include "MeshFactory.hh"
#include "Mesh.hh"
#include "Point.hh"

#include "MFD3D_Electromagnetics.hh"
#include "Tensor.hh"
#include "VEM_NedelecSerendipityType2.hh"


/* ******************************************************************
* Mass matrix in 2D
****************************************************************** */
TEST(MASS_MATRIX_2D) {
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::WhetStone;

  std::cout << "\nTest: Mass matrix for edge elements in 2D" << std::endl;
#ifdef HAVE_MPI
  auto comm = Amanzi::getDefaultComm();
#else
  auto comm = Amanzi::getCommSelf();
#endif

  MeshFactory meshfactory(comm);
  meshfactory.set_preference(Preference({Framework::MSTK}));

  bool request_faces(true), request_edges(true);
  // Teuchos::RCP<Mesh> mesh = meshfactory.create(0.0, 0.0, 1.0, 1.0, 20, 20, true, true); 
  Teuchos::RCP<Mesh> mesh = meshfactory.create("test/two_cell2.exo", request_faces, request_edges); 
 
  Teuchos::ParameterList plist;
  plist.set<int>("method order", 0);
  MFD3D_Electromagnetics mfd(plist, mesh);
  VEM_NedelecSerendipityType2 vem(plist, mesh);

  int cell = 0;
  AmanziMesh::Entity_ID_List edges;
  mesh->cell_get_edges(cell, &edges);

  int nedges = edges.size();

  Tensor T(2, 2);
  T(0, 0) = 2.0;
  T(1, 1) = 1.0;
  T(0, 1) = 1.0;
  T(1, 0) = 1.0;

  for (int method = 0; method < 6; method++) {
    DenseMatrix M;

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
    } else if (method == 4) {
      vem.set_order(0);
      vem.MassMatrix(cell, T, M);
    } else if (method == 5) {
      vem.set_order(1);
      vem.MassMatrix(cell, T, M);
    }

    printf("Mass matrix for cell %3d method=%d\n", cell, method);
    int nrows = M.NumRows();
    for (int i = 0; i < nrows; i++) {
      for (int j = 0; j < nrows; j++ ) printf("%9.5f ", M(i, j)); 
      printf("\n");
    }

    // verify SPD propery
    for (int i = 0; i < nrows; i++) CHECK(M(i, i) > 0.0);

    // verify exact integration property
    WhetStone::DenseVector u(nrows), v(nrows), w(nrows);
    u.PutScalar(0.0);
    v.PutScalar(0.0);

    int stride = (method == 5) ? 2 : 1;
    for (int k = 0, i = 0; i < nedges; ++i, k += stride) {
      int e = edges[i];
      const AmanziGeometry::Point& tau = mesh->edge_vector(e);
      u(k) = tau[0] / mesh->edge_length(e);
      v(k) = tau[1] / mesh->edge_length(e);
    }

    double volume = mesh->cell_volume(cell); 
    M.Multiply(u, w, false);
    double vxx = u * w;
    double vxy = v * w;

    CHECK_CLOSE(volume, vxx, 1e-10);
    CHECK_CLOSE(-volume, vxy, 1e-10);
  }
}


/* ******************************************************************
* Mass matrix in 3D
****************************************************************** */
void MassMatrix3D(std::string mesh_file, int max_row) {
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::WhetStone;

  std::cout << "\nTest: Mass matrix for edge elements in 3D: " << mesh_file << std::endl;
#ifdef HAVE_MPI
  auto comm = Amanzi::getDefaultComm();
#else
  auto comm = Amanzi::getCommSelf();
#endif

  MeshFactory meshfactory(comm);
  meshfactory.set_preference(Preference({Framework::MSTK}));

  RCP<Mesh> mesh;
  if (mesh_file == "")
    mesh = meshfactory.create(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1, 1, 1, true, true); 
  else 
    mesh = meshfactory.create(mesh_file, true, true); 
 
  Teuchos::ParameterList plist;
  plist.set<int>("method order", 0);

  MFD3D_Electromagnetics mfd(plist, mesh);
  VEM_NedelecSerendipityType2 vem(plist, mesh);

  int cell = 0;
  AmanziMesh::Entity_ID_List edges;
  mesh->cell_get_edges(cell, &edges);

  int nedges = edges.size();

  Tensor T(3, 2);
  T(0, 0) = 2.0;
  T(1, 1) = 1.0;
  T(0, 1) = 1.0;
  T(1, 0) = 1.0;
  T(2, 2) = 1.0;

  for (int method = 0; method < 6; method++) {
    DenseMatrix M;

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
    } else if (method == 4) {
      vem.set_order(0);
      vem.MassMatrix(cell, T, M);
    } else if (method == 5) {
      vem.set_order(1);
      vem.MassMatrix(cell, T, M);
    }

    int nrows = M.NumRows();
    int m = std::min(nrows, max_row);
    printf("Mass matrix: method=%d  edges=%d  size=%d  submatrix=%dx%d\n", method, nedges, nrows, m, m);

    for (int i = 0; i < m; i++) {
      for (int j = 0; j < m; j++ ) printf("%8.4f ", M(i, j)); 
      printf("\n");
    }

    // verify SPD propery
    for (int i = 0; i < nrows; i++) CHECK(M(i, i) > 0.0);

    // verify exact integration property
    WhetStone::DenseVector u(nrows), v(nrows), w(nrows);
    u.PutScalar(0.0);
    v.PutScalar(0.0);

    int stride = (method == 5) ? 2 : 1;
    for (int k = 0, i = 0; i < nedges; ++i, k += stride) {
      int e = edges[i];
      const AmanziGeometry::Point& tau = mesh->edge_vector(e);
      u(k) = tau[0] / mesh->edge_length(e);
      v(k) = tau[1] / mesh->edge_length(e);
    }

    double volume = mesh->cell_volume(cell); 
    M.Multiply(u, w, false);
    double vxx = u * w;
    double vxy = v * w;

    CHECK_CLOSE(volume, vxx, 1e-10);
    CHECK_CLOSE(-volume, vxy, 1e-10);
  }
}

TEST(MASS_MATRIX_3D_CUBE) {
  MassMatrix3D("", 12);
exit(0);
}

/*
TEST(MASS_MATRIX_3D_CUBE_ROTATED) {
  MassMatrix3D("test/cube_unit_rotated.exo", 12);
}

TEST(MASS_MATRIX_3D_HEX) {
  MassMatrix3D("test/one_trapezoid.exo", 12);
}

TEST(MASS_MATRIX_3D_24SIDED) {
  MassMatrix3D("test/cube_triangulated.exo", 10);
}

TEST(MASS_MATRIX_3D_DODECAHEDRON) {
  MassMatrix3D("test/dodecahedron.exo", 10);
}
*/


/* ******************************************************************
* Stiffness matrix in 2D
****************************************************************** */
TEST(STIFFNESS_MATRIX_2D) {
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::WhetStone;

  std::cout << "\nTest: Stiffness matrix for edge elements in 2D" << std::endl;
#ifdef HAVE_MPI
  auto comm = Amanzi::getDefaultComm();
#else
  auto comm = Amanzi::getCommSelf();
#endif

  MeshFactory meshfactory(comm);
  meshfactory.set_preference(Preference({Framework::MSTK}));

  bool request_faces(true), request_edges(true);
  // Teuchos::RCP<Mesh> mesh = meshfactory.create(0.0, 0.0, 1.0, 1.0, 1, 1, true, true); 
  Teuchos::RCP<Mesh> mesh = meshfactory.create("test/two_cell2.exo", request_faces, request_edges); 
 
  Teuchos::ParameterList plist;
  MFD3D_Electromagnetics mfd(plist, mesh);

  int cell = 0;
  AmanziMesh::Entity_ID_List faces;
  std::vector<int> dirs;
  mesh->cell_get_faces_and_dirs(cell, &faces, &dirs);

  int nfaces = faces.size();
  int nrows = nfaces;

  Tensor T(2, 1);
  T(0, 0) = 1.0;

  for (int method = 1; method < 2; method++) {
    DenseMatrix A(nrows, nrows);

    if (method == 0) {
      mfd.StiffnessMatrix(cell, T, A);
    } else if (method == 1) {
      mfd.StiffnessMatrix_GradCorrection(cell, T, A);
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
}


/* ******************************************************************
* Stiffness matrix in 3D
****************************************************************** */
void StiffnessMatrix3D(std::string mesh_file, int max_row) {
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::WhetStone;

  std::cout << "\nTest: Stiffness matrix for edge elements in 3D: " << mesh_file << std::endl;
#ifdef HAVE_MPI
  auto comm = Amanzi::getDefaultComm();
#else
  auto comm = Amanzi::getCommSelf();
#endif

  MeshFactory meshfactory(comm);
  meshfactory.set_preference(Preference({Framework::MSTK}));

  RCP<Mesh> mesh;
  if (mesh_file == "")
    mesh = meshfactory.create(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1, 1, 1, true, true); 
  else 
    mesh = meshfactory.create(mesh_file, true, true); 
 
  Teuchos::ParameterList plist;
  MFD3D_Electromagnetics mfd(plist, mesh);

  int cell = 0;
  AmanziMesh::Entity_ID_List edges;
  mesh->cell_get_edges(cell, &edges);

  int nedges = edges.size();
  int nrows = nedges;

  Tensor T(3, 2);
  T(0, 0) = 2.0;
  T(1, 1) = 1.0;
  T(2, 2) = 3.0;
  T(0, 1) = 1.0;
  T(1, 0) = 1.0;

  for (int method = 0; method < 2; method++) {
    DenseMatrix A(nrows, nrows);

    if (method == 0) {
      mfd.StiffnessMatrix(cell, T, A);
    } else if (method == 1) {
      mfd.StiffnessMatrix_GradCorrection(cell, T, A);
    }

    int m = std::min(nrows, max_row);
    printf("Stiffness matrix: method=%d  edges=%d  submatrix=%dx%d\n", method, nedges, m, m);

    for (int i = 0; i < m; i++) {
      for (int j = 0; j < m; j++ ) printf("%9.5f ", A(i, j)); 
      printf("\n");
    }

    // verify SPD propery
    for (int i = 0; i < nrows; i++) CHECK(A(i, i) > 0.0);

    // verify exact integration property
    int n1, n2;
    double xi, xj, yj;
    double vxx(0.0), vxy(0.0), volume = mesh->cell_volume(cell); 
    AmanziGeometry::Point v1(3);
    const AmanziGeometry::Point& xc = mesh->cell_centroid(cell); 

    for (int i = 0; i < nedges; i++) {
      int e1 = edges[i];
      const AmanziGeometry::Point& xe = mesh->edge_centroid(e1);
      const AmanziGeometry::Point& t1 = mesh->edge_vector(e1);
      double a1 = mesh->edge_length(e1);

      v1 = xe ^ t1;
      xi = v1[0] / a1;

      for (int j = 0; j < nedges; j++) {
        int e2 = edges[j];
        const AmanziGeometry::Point& ye = mesh->edge_centroid(e2);
        const AmanziGeometry::Point& t2 = mesh->edge_vector(e2);
        double a2 = mesh->edge_length(e2);

        v1 = ye ^ t2;
        xj = v1[0] / a2;
        yj = v1[1] / a2;

        vxx += A(i, j) * xi * xj;
        vxy += A(i, j) * xi * yj;
      }
    }
    double tol = vxx * 1e-10;
    CHECK_CLOSE(4 * volume * T(0,0), vxx, tol);
    CHECK_CLOSE(4 * volume * T(0,1), vxy, tol);
  }
}

/*
TEST(STIFFNESS_MATRIX_3D_CUBE) {
  StiffnessMatrix3D("", 12);
}

TEST(STIFFNESS_MATRIX_3D_HEX) {
  StiffnessMatrix3D("test/one_trapezoid.exo", 12);
}

TEST(STIFFNESS_MATRIX_3D_DODECAHEDRON) {
  StiffnessMatrix3D("test/dodecahedron.exo", 10);
}

TEST(STIFFNESS_MATRIX_3D_24SIDES) {
  StiffnessMatrix3D("test/cube_triangulated.exo", 10);
} 
*/


