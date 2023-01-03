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
#include "Tensor.hh"


/* **************************************************************** */
TEST(ELASTICITY_STIFFNESS_2D)
{
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::WhetStone;

  std::cout << "\nTest: Stiffness matrix for elasticity in 2D" << std::endl;
#ifdef HAVE_MPI
  auto comm = Amanzi::getDefaultComm();
#else
  auto comm = Amanzi::getCommSelf();
#endif

  MeshFactory meshfactory(comm);
  meshfactory.set_preference(Preference({ Framework::MSTK }));
  // RCP<Mesh> mesh = meshfactory.create(0.0, 0.0, 1.0, 1.0, 1, 1);
  RCP<Mesh> mesh = meshfactory.create("test/one_pentagon.exo");

  Teuchos::ParameterList plist;
  plist.set<std::string>("base", "cell");
  MFD3D_Elasticity mfd(plist, mesh);

  int nnodes = 5, nrows = nnodes * 2, cell = 0;
  DenseMatrix A(nrows, nrows);

  for (int method = 0; method < 2; method++) {
    Tensor T;
    double lambda(1.0), mu(0.0);
    if (method == 0) {
      T.Init(2, 1);
      T(0, 0) = 1.0;
    } else if (method == 1) {
      mu = 3.0;
      lambda = 0.0;
      T.Init(2, 4);
      T(0, 0) = T(1, 1) = lambda + 2 * mu;
      T(0, 1) = T(1, 0) = lambda;
      T(2, 2) = T(3, 3) = mu;
    }

    // mfd.StiffnessMatrix(cell, T, A);
    mfd.StiffnessMatrixOptimized(cell, T, A);
    // mfd.StiffnessMatrixMMatrix(cell, T, A);
    // std::cout << "Number of simplex itrs=" << mfd.simplex_num_itrs() << std::endl;
    // std::cout << "Functional value=" << mfd.simplex_functional() << std::endl;

    printf("Stiffness matrix for cell %3d\n", cell);
    PrintMatrix(A, "%8.4f ");

    // verify SPD propery
    for (int i = 0; i < nrows; i++) CHECK(A(i, i) > 0.0);

    // verify exact integration property
    AmanziMesh::Entity_ID_List nodes;
    nodes = mesh->getCellNodes(cell);

    int d = mesh->getSpaceDimension();
    Point p(d);

    double xx[nrows], yy[nrows];
    for (int i = 0; i < nnodes; i++) {
      int v = nodes[i];
      p = mesh->getNodeCoordinate(v);
      xx[2 * i] = p[0];
      xx[2 * i + 1] = 0.0;

      yy[2 * i] = 0.0;
      yy[2 * i + 1] = p[1];
    }

    double vxx = 0.0, vxy = 0.0, volume = mesh->getCellVolume(cell);
    for (int i = 0; i < nrows; i++) {
      for (int j = 0; j < nrows; j++) {
        vxx += A(i, j) * xx[i] * xx[j];
        vxy += A(i, j) * yy[i] * xx[j];
      }
    }
    CHECK_CLOSE((2 * mu + lambda) * volume, vxx, 1e-10);
    CHECK_CLOSE(0.0, vxy, 1e-10);
  }
}


/* **************************************************************** */
TEST(ELASTICITY_STIFFNESS_3D)
{
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::WhetStone;

  std::cout << "\nTest: Stiffness matrix for Elasticity in 3D" << std::endl;
  auto comm = Amanzi::getDefaultComm();

  MeshFactory meshfactory(comm);
  meshfactory.set_preference(Preference({ Framework::MSTK }));
  // RCP<Mesh> mesh = meshfactory.create(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1, 1, 1);
  RCP<Mesh> mesh = meshfactory.create("test/one_trapezoid.exo");

  Teuchos::ParameterList plist;
  plist.set<std::string>("base", "cell");
  MFD3D_Elasticity mfd(plist, mesh);

  int nrows = 24, nnodes = 8, cell = 0;
  Tensor T(3, 1);
  T(0, 0) = 1;

  DenseMatrix A(nrows, nrows);
  mfd.StiffnessMatrix(cell, T, A);

  printf("Stiffness matrix for cell %3d\n", cell);
  PrintMatrix(A, "%8.4f ");

  // verify SPD propery
  for (int i = 0; i < nrows; i++) CHECK(A(i, i) > 0.0);

  // verify exact integration property
  AmanziMesh::Entity_ID_List nodes;
  nodes = mesh->getCellNodes(cell);

  int d = mesh->getSpaceDimension();
  Point p(d);

  double xx[nrows], yy[nrows];
  for (int i = 0; i < nnodes; i++) {
    int v = nodes[i];
    p = mesh->getNodeCoordinate(v);
    xx[3 * i] = p[0];
    xx[3 * i + 1] = 0.0;
    xx[3 * i + 2] = 0.0;

    yy[3 * i] = 0.0;
    yy[3 * i + 1] = p[1];
    yy[3 * i + 2] = 0.0;
  }

  double vxx = 0.0, vxy = 0.0, volume = mesh->getCellVolume(cell);
  for (int i = 0; i < nrows; i++) {
    for (int j = 0; j < nrows; j++) {
      vxx += A(i, j) * xx[i] * xx[j];
      vxy += A(i, j) * xx[i] * yy[j];
    }
  }
  CHECK_CLOSE(vxx, volume, 1e-10);
  CHECK_CLOSE(vxy, 0.0, 1e-10);
}
