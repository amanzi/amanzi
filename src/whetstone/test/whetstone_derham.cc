/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  WhetStone, Version 2.2
  Release name: naka-to.

  The mimetic finite difference method.
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

#include "DeRham_Node.hh"
#include "Tensor.hh"


/* ****************************************************************
* DeRham complex for nodes
**************************************************************** */
TEST(DERHAM_COMPLEX_NODE)
{
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::WhetStone;

  std::cout << "Test: Mass matrix for nodes" << std::endl;
  auto comm = Amanzi::getDefaultComm();

  MeshFactory meshfactory(comm);
  meshfactory.set_preference(Preference({ Framework::MSTK }));
  Teuchos::RCP<Mesh> mesh = meshfactory.create(0.0, 0.0, 0.5, 1.0, 1, 1);

  DeRham_Node drc(mesh);

  int nnodes(4), cell(0);
  DenseMatrix M(nnodes, nnodes);

  Tensor T(2, 1);
  T(0, 0) = 1.0;

  drc.MassMatrix(cell, T, M);

  printf("Mass matrix for cell %3d\n", cell);
  PrintMatrix(M, "%8.4f ");

  // verify SPD propery
  for (int i = 0; i < nnodes; ++i) CHECK(M(i, i) > 0.0);

  // verify exact integration property
  auto nodes = mesh->getCellNodes(cell);

  double xi, yi, xj;
  double vxx = 0.0, vxy = 0.0, volume = mesh->getCellVolume(cell);
  AmanziGeometry::Point p1(2), p2(2);

  for (int i = 0; i < nnodes; i++) {
    int v1 = nodes[i];
    p1 = mesh->getNodeCoordinate(v1);
    for (int j = 0; j < nnodes; j++) {
      xi = 1.0; // p1[0];
      yi = 1.0; // p1[1];
      xj = 1.0; // p2[0];

      vxx += M(i, j) * xi * xj;
      vxy += M(i, j) * yi * xj;
    }
  }

  CHECK_CLOSE(T(0, 0) * volume, vxx, 1e-10);
  CHECK_CLOSE(T(0, 0) * volume, vxy, 1e-10);
}
