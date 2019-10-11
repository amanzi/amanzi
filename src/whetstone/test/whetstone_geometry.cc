/*
  WhetStone, version 2.0
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <cstdlib>
#include <cmath>
#include <iostream>

#include "UnitTest++.h"

#include "AmanziComm.hh"
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
  Kokkos::View<Entity_ID*> nodes;
  std::vector<double> weights;

  double area = mesh->cell_volume(0, false);
  mesh->cell_get_nodes(0, nodes);
  PolygonCentroidWeights(*mesh, nodes, area, weights);

  // verification
  for (int n = 0; n < nodes.extent(0); ++n) {
    mesh->node_get_coordinates(nodes(n), &p);
    xc += weights[n] * p;
  }
  const AmanziGeometry::Point& xm = mesh->cell_centroid(0);
  std::cout << xc << " = " << xm << std::endl;
  CHECK_CLOSE(0.0, norm(xc - xm), 1e-10);
}
