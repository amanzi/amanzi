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

  int c = 0;
  AmanziMesh::Entity_ID_List edges;
  mesh->cell_get_edges(c, &edges);

  int nedges = edges.size();

  Tensor T(2, 2);
  T(0, 0) = 2.0;
  T(1, 1) = 1.0;
  T(0, 1) = 1.0;
  T(1, 0) = 1.0;

  for (int method = 0; method < 6; method++) {
    DenseMatrix M;

    if (method == 0) {
      mfd.MassMatrix(c, T, M);
    } else if (method == 1) {
      mfd.MassMatrixOptimized(c, T, M);
    } else if (method == 2) {
      mfd.MassMatrixInverse(c, T, M);
      M.Inverse();
    } else if (method == 3) {
      mfd.MassMatrixInverseOptimized(c, T, M);
      M.Inverse();
    } else if (method == 4) {
      vem.set_order(0);
      vem.MassMatrix(c, T, M);
    } else if (method == 5) {
      vem.set_order(1);
      vem.MassMatrix(c, T, M);
    }

    printf("Mass matrix for cell %3d method=%d\n", c, method);
    int nrows = M.NumRows();
    for (int i = 0; i < nrows; i++) {
      for (int j = 0; j < nrows; j++ ) printf("%9.5f ", M(i, j)); 
      printf("\n");
    }

    // verify SPD propery
    for (int i = 0; i < nrows; i++) CHECK(M(i, i) > 0.0);

    // verify exact integration property
    int k = (method == 5) ? 1 : 0;
    std::vector<VectorPolynomial> uf(nedges), vf(nedges), wf(nedges);
    for (int i = 0; i < nedges; ++i) {
      uf[i].Reshape(2, 2, 0);
      uf[i][0](0) = 1.0;

      vf[i].Reshape(2, 2, 0);
      vf[i][1](0) = 1.0;

      wf[i].Reshape(2, 2, k);
      wf[i][0](k) = 1.0;
    }

    WhetStone::DenseVector u(nrows), v(nrows), w(nrows), a(nrows);
    vem.CalculateDOFsOnBoundary<AmanziMesh::Mesh>(mesh, c, uf, uf, u);
    vem.CalculateDOFsOnBoundary<AmanziMesh::Mesh>(mesh, c, vf, vf, v);
    vem.CalculateDOFsOnBoundary<AmanziMesh::Mesh>(mesh, c, wf, wf, w);

    M.Multiply(w, a, false);
    double vxx = u * a;
    double vxy = v * a;

    double volume = mesh->cell_volume(c); 
    const AmanziGeometry::Point& xc = mesh->cell_centroid(c);
    if (method != 5) {
      CHECK_CLOSE(volume, vxx, 1e-10);
      CHECK_CLOSE(-volume, vxy, 1e-10);
    } else {
      CHECK_CLOSE(volume * xc[0], vxx, 1e-10);
      CHECK_CLOSE(-volume * xc[0], vxy, 1e-10);
    }
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
    mesh = meshfactory.create(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 2, 2, 2, true, true); 
  else 
    mesh = meshfactory.create(mesh_file, true, true); 
 
  // loop over cells will test orientation of edges
  int ncells = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  for (int c = 0; c < ncells; ++c) {
    Teuchos::ParameterList plist;
    plist.set<int>("method order", 0);

    MFD3D_Electromagnetics mfd(plist, mesh);
    VEM_NedelecSerendipityType2 vem(plist, mesh);

    AmanziMesh::Entity_ID_List edges;
    mesh->cell_get_edges(c, &edges);
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
        mfd.MassMatrix(c, T, M);
      } else if (method == 1) {
        mfd.MassMatrixOptimized(c, T, M);
      } else if (method == 2) {
        mfd.MassMatrixInverse(c, T, M);
        M.Inverse();
      } else if (method == 3) {
        mfd.MassMatrixInverseOptimized(c, T, M);
        M.Inverse();
      } else if (method == 4) {
        vem.set_order(0);
        vem.MassMatrix(c, T, M);
      } else if (method == 5) {
        vem.set_order(1);
        vem.MassMatrix(c, T, M);
      }

      int nrows = M.NumRows();
      int m = std::min(nrows, max_row);
      printf("Mass matrix: cell=%d method=%d  edges=%d  size=%d  submatrix=%dx%d\n", c, method, nedges, nrows, m, m);

      for (int i = 0; i < m; i++) {
        for (int j = 0; j < m; j++ ) printf("%8.4f ", M(i, j)); 
        printf("\n");
      }

      // verify SPD propery
      for (int i = 0; i < nrows; i++) CHECK(M(i, i) > 0.0);

      // verify exact integration property
      std::vector<VectorPolynomial> uf(nedges), vf(nedges), wf(nedges);
      const AmanziGeometry::Point& xc = mesh->cell_centroid(c);

      int k = (method == 5) ? 1 : 0;
      for (int i = 0; i < nedges; ++i) {
        uf[i].Reshape(3, 3, 0);
        uf[i][0](0) = 1.0;
        uf[i].set_origin(xc);

        vf[i].Reshape(3, 3, 0);
        vf[i][1](0) = 1.0;
        vf[i].set_origin(xc);

        wf[i].Reshape(3, 3, k);
        wf[i][0](k) = 1.0;
        wf[i].set_origin(xc);
      }

      WhetStone::DenseVector u(nrows), v(nrows), w(nrows), a(nrows);
      vem.CalculateDOFsOnBoundary<AmanziMesh::Mesh>(mesh, c, uf, uf, u);
      vem.CalculateDOFsOnBoundary<AmanziMesh::Mesh>(mesh, c, vf, vf, v);
      vem.CalculateDOFsOnBoundary<AmanziMesh::Mesh>(mesh, c, wf, wf, w);

      M.Multiply(w, a, false);
      double vx1 = u * a;
      double vy1 = v * a;
      double vxx = w * a;

      double volume = mesh->cell_volume(c); 
      if (method != 5) {
        CHECK_CLOSE(volume, vx1, 1e-10);
        CHECK_CLOSE(-volume, vy1, 1e-10);
      } else {
        CHECK_CLOSE(0.0, vx1, 1e-10);
        CHECK_CLOSE(0.0, vy1, 1e-10);
        CHECK_CLOSE(vem.integrals().poly()(4), vxx, 1e-10 * vxx);
      }
    }
  }
}

TEST(MASS_MATRIX_3D_CUBE) {
  MassMatrix3D("", 12);
}

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

  // loop over cells will test orientation of edges
  int c = 0;
  AmanziMesh::Entity_ID_List faces;
  std::vector<int> dirs;
  mesh->cell_get_faces_and_dirs(c, &faces, &dirs);

  int nfaces = faces.size();
  int nrows = nfaces;

  Tensor T(2, 1);
  T(0, 0) = 1.0;

  for (int method = 1; method < 2; method++) {
    DenseMatrix A(nrows, nrows);

    if (method == 0) {
      mfd.StiffnessMatrix(c, T, A);
    } else if (method == 1) {
      mfd.StiffnessMatrix_GradCorrection(c, T, A);
    }

    printf("Stiffness matrix for cell %3d method=%d\n", c, method);
    for (int i = 0; i < nrows; i++) {
      for (int j = 0; j < nrows; j++ ) printf("%8.4f ", A(i, j)); 
        printf("\n");
    }

    // verify SPD propery
    for (int i = 0; i < nrows; i++) CHECK(A(i, i) > 0.0);

    // verify exact integration property
    double xi, xj, vxx(0.0);
    AmanziGeometry::Point p1(2), p2(2);

    const AmanziGeometry::Point& xc = mesh->cell_centroid(c);
    double volume = mesh->cell_volume(c); 

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
    mesh = meshfactory.create(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 2, 2, 2, true, true); 
  else 
    mesh = meshfactory.create(mesh_file, true, true); 
 
  int ncells = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  for (int c = 0; c < 1; ++c) {
  // for (int c = 0; c < ncells; ++c) {
    Teuchos::ParameterList plist;
    plist.set<int>("method order", 0);

    MFD3D_Electromagnetics mfd(plist, mesh);
    VEM_NedelecSerendipityType2 vem(plist, mesh);

    AmanziMesh::Entity_ID_List edges;
    mesh->cell_get_edges(c, &edges);

    int nedges = edges.size();
    int nrows = nedges;

    Tensor T(3, 2);
    T(0, 0) = 2.0;
    T(1, 1) = 1.0;
    T(2, 2) = 3.0;
    T(0, 1) = 1.0;
    T(1, 0) = 1.0;

    for (int method = 0; method < 3; method++) {
      DenseMatrix A(nrows, nrows);

      if (method == 0) {
        mfd.StiffnessMatrix(c, T, A);
      } else if (method == 1) {
        mfd.StiffnessMatrix_GradCorrection(c, T, A);
      } else if (method == 2) {
        vem.set_order(0);
        vem.StiffnessMatrix(c, T, A);
      }

      int m = std::min(nrows, max_row);
      printf("Stiffness matrix: cell=%d  method=%d  edges=%d  submatrix=%dx%d\n", c, method, nedges, m, m);

      for (int i = 0; i < m; i++) {
        for (int j = 0; j < m; j++ ) printf("%9.5f ", A(i, j)); 
        printf("\n");
      }

      // verify SPD propery
      for (int i = 0; i < nrows; i++) CHECK(A(i, i) > 0.0);

      // verify exact integration property
      double xi, xj, yj;
      double vxx(0.0), vxy(0.0), volume = mesh->cell_volume(c); 
      AmanziGeometry::Point v1(3);

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
}

TEST(STIFFNESS_MATRIX_3D_CUBE) {
  StiffnessMatrix3D("", 12);
}

/*
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


