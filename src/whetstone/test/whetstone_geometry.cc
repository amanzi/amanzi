/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  WhetStone, version 2.0
  Release name: naka-to.

*/

#include <cstdlib>
#include <cmath>
#include <iostream>

#include "UnitTest++.h"

#include "Mesh.hh"
#include "MeshFactory.hh"
#include "WhetStoneMeshUtils.hh"

/* ****************************************************************
* Test of face centroids
**************************************************************** */
TEST(FACE_CENTROIDS)
{
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::WhetStone;

  std::cout << "Test: calculation of face centroids" << std::endl;

  Preference pref;
  pref.clear();
  pref.push_back(Framework::MSTK);

  auto comm = Amanzi::getDefaultComm();
  MeshFactory meshfactory(comm);
  meshfactory.set_preference(pref);
  // Teuchos::RCP<Mesh> mesh = meshfactory.create(0.0, 0.0, 1.0, 1.0, 1, 1);
  Teuchos::RCP<Mesh> mesh = meshfactory.create("test/one_pentagon.exo");

  AmanziGeometry::Point p(2), xc(2);
  Entity_ID_List nodes;
  std::vector<double> weights;

  double area = mesh->getCellVolume(0);
  nodes = mesh->getCellNodes(0);
  PolygonCentroidWeights(*mesh, nodes, area, weights);

  // verification
  for (int n = 0; n < nodes.size(); ++n) {
    p = mesh->getNodeCoordinate(nodes[n]);
    xc += weights[n] * p;
  }
  const AmanziGeometry::Point& xm = mesh->getCellCentroid(0);
  std::cout << xc << " = " << xm << std::endl;
  CHECK_CLOSE(0.0, norm(xc - xm), 1e-10);
}
