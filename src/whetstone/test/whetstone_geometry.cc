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

#include "Mesh.hh"
#include "MeshFactory.hh"

#include "MFD3D_Diffusion.hh"

/* ****************************************************************
* Test of face centroids
**************************************************************** */
TEST(FACE_CENTROIDS) {
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::WhetStone;

  std::cout << "Test: calculation of face centroids" << std::endl;
 
  FrameworkPreference pref;
  pref.clear();
  pref.push_back(MSTK);

  Epetra_MpiComm *comm = new Epetra_MpiComm(MPI_COMM_WORLD);
  MeshFactory factory(comm);
  factory.preference(pref);
  // Teuchos::RCP<Mesh> mesh = factory(0.0, 0.0, 1.0, 1.0, 1, 1); 
  Teuchos::RCP<Mesh> mesh = factory("test/one_cell2.exo");
 
  MFD3D_Diffusion mfd(mesh);
 
  AmanziGeometry::Point p(2), xc(2);
  Entity_ID_List nodes;
  std::vector<double> weights;

  double area = mesh->cell_volume(0);
  mesh->cell_get_nodes(0, &nodes); 
  mfd.PolygonCentroidWeights(nodes, area, weights);
  
  // verification
  for (int n = 0; n < nodes.size(); ++n) {
    mesh->node_get_coordinates(nodes[n], &p);
    xc += weights[n] * p;
  }
  const AmanziGeometry::Point& xm = mesh->cell_centroid(0);
  std::cout << xc << " = " << xm << std::endl;
  CHECK_CLOSE(0.0, norm(xc - xm), 1e-10);

  delete comm;
}


