/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  This is the mimetic discretization component of the Amanzi code.

*/

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <vector>

#include "Teuchos_RCP.hpp"
#include "UnitTest++.h"

#include "MeshFactory.hh"
#include "Mesh.hh"
#include "Point.hh"

#include "MFD3D_Electromagnetics.hh"
#include "Tensor.hh"

#if 0
/* ******************************************************************
* Mass matrix in 2D
****************************************************************** */
TEST(MASS_MATRIX_2D)
{
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::WhetStone;

  std::cout << "\nTest: Mass matrix for edge elements in 2D" << std::endl;
#  ifdef HAVE_MPI
  auto comm = Amanzi::getDefaultComm();
#  else
  auto comm = Amanzi::getCommSelf();
#  endif
  Teuchos::RCP<Teuchos::ParameterList> factory_plist = Teuchos::rcp(new Teuchos::ParameterList()); 
  factory_plist->set<bool>("request edges", true); 
  factory_plist->set<bool>("request faces", true); 
  MeshFactory meshfactory(comm,Teuchos::null,factory_plist);
  meshfactory.set_preference(Preference({ Framework::MSTK }));

  Teuchos::RCP<Mesh> mesh = meshfactory.create("test/two_cell2.exo");

  Teuchos::ParameterList plist;
  MFD3D_Electromagnetics mfd(plist, mesh);

  int cell = 0;
  auto edges = mesh->getCellEdges(cell);

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
      for (int j = 0; j < nrows; j++) printf("%9.5f ", M(i, j));
      printf("\n");
    }

    // verify SPD propery
    for (int i = 0; i < nrows; i++) CHECK(M(i, i) > 0.0);

    // verify exact integration property
    double xi, yi, xj;
    double vxx = 0.0, vxy = 0.0, volume = mesh->getCellVolume(cell);
    for (int i = 0; i < nedges; i++) {
      int e1 = edges[i];
      const AmanziGeometry::Point& t1 = mesh->getEdgeVector(e1);
      double a1 = mesh->getEdgeLength(e1);

      xi = t1[0] / a1;
      yi = t1[1] / a1;
      for (int j = 0; j < nedges; j++) {
        int e2 = edges[j];
        const AmanziGeometry::Point& t2 = mesh->getEdgeVector(e2);
        double a2 = mesh->getEdgeLength(e2);
        xj = t2[0] / a2;
        vxx += M(i, j) * xi * xj;
        vxy += M(i, j) * yi * xj;
      }
    }
    CHECK_CLOSE(volume, vxx, 1e-10);
    CHECK_CLOSE(-volume, vxy, 1e-10);
  }
}

#endif

/* ******************************************************************
* Mass matrix in 3D
****************************************************************** */
void
MassMatrix3D(std::string mesh_file, int max_row)
{
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

  Teuchos::RCP<Teuchos::ParameterList> factory_plist = Teuchos::rcp(new Teuchos::ParameterList());
  factory_plist->set<bool>("request edges", true);
  factory_plist->set<bool>("request faces", true);
  MeshFactory meshfactory(comm, Teuchos::null, factory_plist);
  meshfactory.set_preference(Preference({ Framework::MSTK }));

  RCP<Mesh> mesh;
  if (mesh_file == "")
    mesh = meshfactory.create(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1, 1, 1);
  else
    mesh = meshfactory.create(mesh_file);

  Teuchos::ParameterList plist;
  MFD3D_Electromagnetics mfd(plist, mesh);

  int cell = 0;
  auto edges = mesh->getCellEdges(cell);

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

    int m = std::min(nrows, max_row);
    printf("Mass matrix: method=%d  edges=%d  submatrix=%dx%d\n", method, nedges, m, m);

    for (int i = 0; i < m; i++) {
      for (int j = 0; j < m; j++) printf("%8.4f ", M(i, j));
      printf("\n");
    }

    // verify SPD propery
    for (int i = 0; i < nrows; i++) CHECK(M(i, i) > 0.0);

    // verify exact integration property
    double xi, yi, xj;
    double vxx = 0.0, vxy = 0.0, volume = mesh->getCellVolume(cell);
    for (int i = 0; i < nedges; i++) {
      int e1 = edges[i];
      const AmanziGeometry::Point& t1 = mesh->getEdgeVector(e1);
      double a1 = mesh->getEdgeLength(e1);

      xi = t1[0] / a1;
      yi = t1[1] / a1;
      for (int j = 0; j < nedges; j++) {
        int e2 = edges[j];
        const AmanziGeometry::Point& t2 = mesh->getEdgeVector(e2);
        double a2 = mesh->getEdgeLength(e2);
        xj = t2[0] / a2;
        vxx += M(i, j) * xi * xj;
        vxy += M(i, j) * yi * xj;
      }
    }
    CHECK_CLOSE(volume, vxx, 1e-10);
    CHECK_CLOSE(-volume, vxy, 1e-10);
  }
}

TEST(MASS_MATRIX_3D_CUBE)
{
  MassMatrix3D("", 12);
}

#if 0 

TEST(MASS_MATRIX_3D_HEX)
{
  MassMatrix3D("test/one_trapezoid.exo", 12);
}

TEST(MASS_MATRIX_3D_24SIDED)
{
  MassMatrix3D("test/cube_triangulated.exo", 10);
}

TEST(MASS_MATRIX_3D_DODECAHEDRON)
{
  MassMatrix3D("test/dodecahedron.exo", 10);
}


/* ******************************************************************
* Stiffness matrix in 2D
****************************************************************** */
TEST(STIFFNESS_MATRIX_2D)
{
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::WhetStone;

  std::cout << "\nTest: Stiffness matrix for edge elements in 2D" << std::endl;
#  ifdef HAVE_MPI
  auto comm = Amanzi::getDefaultComm();
#  else
  auto comm = Amanzi::getCommSelf();
#  endif

  Teuchos::RCP<Teuchos::ParameterList> factory_plist = Teuchos::rcp(new Teuchos::ParameterList()); 
  factory_plist->set<bool>("request edges", true); 
  factory_plist->set<bool>("request faces", true); 
  MeshFactory meshfactory(comm,Teuchos::null,factory_plist);
  meshfactory.set_preference(Preference({ Framework::MSTK }));

  Teuchos::RCP<Mesh> mesh = meshfactory.create("test/two_cell2.exo");

  Teuchos::ParameterList plist;
  MFD3D_Electromagnetics mfd(plist, mesh);

  int cell = 0;
  auto [faces, dirs] = mesh->getCellFacesAndDirections(cell);

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
      for (int j = 0; j < nrows; j++) printf("%8.4f ", A(i, j));
      printf("\n");
    }

    // verify SPD propery
    for (int i = 0; i < nrows; i++) CHECK(A(i, i) > 0.0);

    // verify exact integration property
    double xi, xj, vxx(0.0);
    AmanziGeometry::Point p1(2), p2(2);

    const AmanziGeometry::Point& xc = mesh->getCellCentroid(cell);
    double volume = mesh->getCellVolume(cell);

    for (int i = 0; i < nrows; i++) {
      int f1 = faces[i];
      const AmanziGeometry::Point& n1 = mesh->getFaceNormal(f1);
      const AmanziGeometry::Point& xf1 = mesh->getFaceCentroid(f1);
      double a1 = mesh->getFaceArea(f1);

      xi = ((xf1 - xc) * n1) * dirs[i] / a1;

      for (int j = 0; j < nrows; j++) {
        int f2 = faces[j];
        const AmanziGeometry::Point& n2 = mesh->getFaceNormal(f2);
        const AmanziGeometry::Point& xf2 = mesh->getFaceCentroid(f2);
        double a2 = mesh->getFaceArea(f2);

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
void
StiffnessMatrix3D(std::string mesh_file, int max_row)
{
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::WhetStone;

  std::cout << "\nTest: Stiffness matrix for edge elements in 3D: " << mesh_file << std::endl;
#  ifdef HAVE_MPI
  auto comm = Amanzi::getDefaultComm();
#  else
  auto comm = Amanzi::getCommSelf();
#  endif

  Teuchos::RCP<Teuchos::ParameterList> factory_plist = Teuchos::rcp(new Teuchos::ParameterList()); 
  factory_plist->set<bool>("request edges", true); 
  factory_plist->set<bool>("request faces", true); 
  MeshFactory meshfactory(comm,Teuchos::null,factory_plist);
  meshfactory.set_preference(Preference({ Framework::MSTK }));

  RCP<Mesh> mesh;
  if (mesh_file == "")
    mesh = meshfactory.create(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1, 1, 1);
  else
    mesh = meshfactory.create(mesh_file);

  Teuchos::ParameterList plist;
  MFD3D_Electromagnetics mfd(plist, mesh);

  int cell = 0;
  auto edges = mesh->getCellEdges(cell);

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
      for (int j = 0; j < m; j++) printf("%9.5f ", A(i, j));
      printf("\n");
    }

    // verify SPD propery
    for (int i = 0; i < nrows; i++) CHECK(A(i, i) > 0.0);

    // verify exact integration property
    double xi, xj, yj;
    double vxx(0.0), vxy(0.0), volume = mesh->getCellVolume(cell);
    AmanziGeometry::Point v1(3);

    for (int i = 0; i < nedges; i++) {
      int e1 = edges[i];
      const AmanziGeometry::Point& xe = mesh->getEdgeCentroid(e1);
      const AmanziGeometry::Point& t1 = mesh->getEdgeVector(e1);
      double a1 = mesh->getEdgeLength(e1);

      v1 = xe ^ t1;
      xi = v1[0] / a1;

      for (int j = 0; j < nedges; j++) {
        int e2 = edges[j];
        const AmanziGeometry::Point& ye = mesh->getEdgeCentroid(e2);
        const AmanziGeometry::Point& t2 = mesh->getEdgeVector(e2);
        double a2 = mesh->getEdgeLength(e2);

        v1 = ye ^ t2;
        xj = v1[0] / a2;
        yj = v1[1] / a2;

        vxx += A(i, j) * xi * xj;
        vxy += A(i, j) * xi * yj;
      }
    }
    double tol = vxx * 1e-10;
    CHECK_CLOSE(4 * volume * T(0, 0), vxx, tol);
    CHECK_CLOSE(4 * volume * T(0, 1), vxy, tol);
  }
}

TEST(STIFFNESS_MATRIX_3D_CUBE)
{
  StiffnessMatrix3D("", 12);
}

TEST(STIFFNESS_MATRIX_3D_HEX)
{
  StiffnessMatrix3D("test/one_trapezoid.exo", 12);
}

TEST(STIFFNESS_MATRIX_3D_DODECAHEDRON)
{
  StiffnessMatrix3D("test/dodecahedron.exo", 10);
}

TEST(STIFFNESS_MATRIX_3D_24SIDES)
{
  StiffnessMatrix3D("test/cube_triangulated.exo", 10);
}
#endif