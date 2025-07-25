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

#include "MeshFactory.hh"
#include "MeshAudit.hh"

#include "Mesh.hh"
#include "Point.hh"

#include "MFD3D_BernardiRaugel.hh"
#include "Tensor.hh"

/* ******************************************************************
* DoFs: vectors at nodes, normal component on faces.
****************************************************************** */
TEST(STIFFNESS_STOKES_2D)
{
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::WhetStone;

  std::cout << "\nTest: Stiffness matrix for Stokes in 2D" << std::endl;
  auto comm = Amanzi::getDefaultComm();

  MeshFactory meshfactory(comm);
  meshfactory.set_preference(Preference({ Framework::MSTK }));
  // RCP<Mesh> mesh = meshfactory.create(0.0, 0.0, 1.0, 1.0, 1, 1);
  RCP<Mesh> mesh = meshfactory.create("test/one_pentagon.exo");

  Teuchos::ParameterList plist;
  MFD3D_BernardiRaugel mfd(plist, mesh);

  // extract single cell
  int cell(0);
  auto nodes = mesh->getCellNodes(cell);
  int nnodes = nodes.size();

  auto [faces, dirs] = mesh->getCellFacesAndDirections(cell);
  int nfaces = faces.size();

  // calcualte stiffness matrix
  Tensor T(2, 1);
  T(0, 0) = 1.0;

  int nrows = 2 * nnodes + nfaces;
  DenseMatrix A(nrows, nrows);

  mfd.StiffnessMatrix(cell, T, A);

  printf("Stiffness matrix for cell %3d\n", cell);
  for (int i = 0; i < nrows; i++) {
    for (int j = 0; j < nrows; j++) printf("%7.4f ", A(i, j));
    printf("\n");
  }

  // verify SPD propery
  for (int i = 0; i < nrows; i++) CHECK(A(i, i) > 0.0);

  // verify exact integration property
  int d = mesh->getSpaceDimension();
  Point p(d);

  DenseVector ax(nrows), ay(nrows);
  ax.PutScalar(0.0);
  ay.PutScalar(0.0);
  for (int i = 0; i < nnodes; i++) {
    int v = nodes[i];
    p = mesh->getNodeCoordinate(v);
    ax(2 * i) = p[0];
    ay(2 * i + 1) = p[1];
  }
  for (int i = 0; i < nfaces; i++) {
    int f = faces[i];
    const AmanziGeometry::Point& normal = mesh->getFaceNormal(f);
    const AmanziGeometry::Point& xf = mesh->getFaceCentroid(f);
    double area = mesh->getFaceArea(f);

    ax(2 * nnodes + i) = xf[0] * normal[0] / area;
    ay(2 * nnodes + i) = xf[1] * normal[1] / area;
  }

  double vxx(0.0), vxy(0.0);
  double volume = mesh->getCellVolume(cell);

  for (int i = 0; i < nrows; i++) {
    for (int j = 0; j < nrows; j++) {
      vxx += A(i, j) * ax(i) * ax(j);
      vxy += A(i, j) * ay(i) * ax(j);
    }
  }
  CHECK_CLOSE(vxx, volume, 1e-10);
  CHECK_CLOSE(vxy, 0.0, 1e-10);
}


/* ******************************************************************
* DoFs: vectors at nodes, normal component on faces.
****************************************************************** */
void
runStiffness3D(const std::string& filename)
{
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::WhetStone;

  std::cout << "\nTest: Stiffness matrix for Stokes in 3D" << std::endl;
  auto comm = Amanzi::getDefaultComm();

  MeshFactory meshfactory(comm);
  meshfactory.set_preference(Preference({ Framework::MSTK }));
  RCP<Mesh> mesh;
  if (filename == "") mesh = meshfactory.create(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1, 2, 3);
  else mesh = meshfactory.create(filename);

  Teuchos::ParameterList plist;
  MFD3D_BernardiRaugel mfd(plist, mesh);

  // extract single cell
  int c(0);
  auto nodes = mesh->getCellNodes(c);
  int nnodes = nodes.size();

  auto [faces, dirs] = mesh->getCellFacesAndDirections(c);
  int nfaces = faces.size();

  // calcualte stiffness matrix
  Tensor T(3, 1);
  T(0, 0) = 1.0;

  int nrows = 3 * nnodes + nfaces;
  DenseMatrix A(nrows, nrows);

  mfd.StiffnessMatrix(c, T, A);

  printf("Stiffness matrix for cell %3d\n", c);
  PrintMatrix(A, "%11.6f", 12);

  // verify SPD propery
  for (int i = 0; i < nrows; i++) CHECK(A(i, i) > 0.0);

  // verify exact integration property
  int d = mesh->getSpaceDimension();
  Point p(d);

  DenseVector ax(nrows), ay(nrows), az(nrows);
  ax.PutScalar(0.0);
  ay.PutScalar(0.0);
  az.PutScalar(0.0);
  for (int i = 0; i < nnodes; i++) {
    int v = nodes[i];
    p = mesh->getNodeCoordinate(v);
    ax(3 * i) = p[0];
    ay(3 * i + 1) = p[1];
    az(3 * i + 2) = p[2];
  }
  for (int i = 0; i < nfaces; i++) {
    int f = faces[i];
    const AmanziGeometry::Point& normal = mesh->getFaceNormal(f);
    const AmanziGeometry::Point& xf = mesh->getFaceCentroid(f);
    double area = mesh->getFaceArea(f);

    ax(3 * nnodes + i) = xf[0] * normal[0] / area;
    ay(3 * nnodes + i) = xf[1] * normal[1] / area;
    az(3 * nnodes + i) = xf[2] * normal[2] / area;
  }

  double vxx(0.0), vxy(0.0), vzz(0.0);
  double volume = mesh->getCellVolume(c);

  for (int i = 0; i < nrows; i++) {
    for (int j = 0; j < nrows; j++) {
      vxx += A(i, j) * ax(i) * ax(j);
      vxy += A(i, j) * ay(i) * ax(j);
      vzz += A(i, j) * az(i) * az(j);
    }
  }
  CHECK_CLOSE(vxx, volume, 1e-10);
  CHECK_CLOSE(vxy, 0.0, 1e-10);
  CHECK_CLOSE(vzz, volume, 1e-10);
}


TEST(STIFFNESS_STOKES_3D_BOXMESH)
{
  runStiffness3D("");
}
TEST(STIFFNESS_STOKES_3D_POLYTOPE)
{
  runStiffness3D("test/dodecahedron.exo");
}


/* ******************************************************************
* Dofs: vectors at nodes, normal component on faces.
****************************************************************** */
TEST(ADVECTION_NAVIER_STOKES_2D)
{
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::WhetStone;

  std::cout << "\nTest: Advection matrix for Navier-Stokes in 2D" << std::endl;
  auto comm = Amanzi::getDefaultComm();

  MeshFactory meshfactory(comm);
  meshfactory.set_preference(Preference({ Framework::MSTK }));
  // RCP<Mesh> mesh = meshfactory.create(0.0, 0.0, 1.0, 1.0, 1, 1);
  RCP<Mesh> mesh = meshfactory.create("test/one_pentagon.exo");

  Teuchos::ParameterList plist;
  MFD3D_BernardiRaugel mfd(plist, mesh);

  // extract single cell
  int cell(0);

  auto nodes = mesh->getCellNodes(cell);
  int nnodes = nodes.size();

  auto [faces, dirs] = mesh->getCellFacesAndDirections(cell);

  // setup velocity
  AmanziMesh::Point_List u(nnodes);
  for (int i = 0; i < nnodes; ++i) {
    u[i] = AmanziGeometry::Point(1.0, 2.0);
  }

  // calculate advection matrix
  DenseMatrix A;
  mfd.AdvectionMatrix(cell, u, A);

  printf("Advection matrix for cell %3d\n", cell);
  for (int i = 0; i < 2 * nnodes; i++) {
    for (int j = 0; j < 2 * nnodes; j++) printf("%8.5f ", A(i, j));
    printf("\n");
  }
}
