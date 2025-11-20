/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#include <iostream>
#include <vector>

#include "Teuchos_RCP.hpp"
#include "UnitTest++.h"

#include "AmanziComm.hh"
#include "MeshDefs.hh"
#include "MeshVirtual.hh"
#include "MeshLogicalAudit.hh"

#include "demo_mesh.hh"

using namespace Amanzi;
using namespace Amanzi::AmanziMesh;
using namespace Amanzi::AmanziGeometry;

TEST(VIRTUAL_MESH_REGULAR)
{
  std::cout << "TEST: Mesh virtual regular" << std::endl;

  auto comm = Amanzi::getDefaultComm();

  // Create the mesh:
  Double_List cell_volumes(4, 0.25);

  std::vector<std::vector<int>> face_cells = {
    { 0, 1 }, { 0, 2 }, { 1, 3 }, { 2, 3},
    { 0 }, { 1 }, { 1 }, { 3 }, { 3 }, { 2 }, { 2 }, { 0 } 
  };
                                             
  std::vector<Point> cell_centroids = {
     Point(0.25, 0.25),
     Point(0.75, 0.25),
     Point(0.25, 0.75),
     Point(0.75, 0.75)
  };

  std::vector<Point> face_centroids = {
    Point(0.5, 0.25), Point(0.25, 0.5), Point(0.75, 0.5), Point(0.5, 0.75),
    Point(0.25, 0.0), Point(0.75, 0.0),
    Point(1.0, 0.25), Point(1.0, 0.75),
    Point(0.75, 1.0), Point(0.25, 1.0),
    Point(0.0, 0.75), Point(0.0, 0.25)
  };

  std::vector<Point> face_normals = { 
    Point(0.5, 0.0), Point(0.5, 0.0), Point(0.0, 0.5), Point(0.5, 0.0),
    Point(0.0, -0.5), Point(0.0, -0.5),
    Point(0.5, 0.0), Point(0.5, 0.0),
    Point(0.0, 0.5), Point(0.0, 0.5),
    Point(-0.5, 0.0), Point(-0.5, 0.0)
  };

  Teuchos::RCP<Teuchos::ParameterList> plist;
  auto mesh_fw = Teuchos::rcp(new MeshVirtual(comm,
                                              plist,
                                              face_cells,
                                              cell_centroids,
                                              cell_volumes,
                                              face_centroids,
                                              face_normals));

  auto mesh = Teuchos::rcp(new Amanzi::AmanziMesh::MeshCache<MemSpace_kind::HOST>(
    mesh_fw, Teuchos::rcp(new Amanzi::AmanziMesh::MeshVirtualAlgorithms()), Teuchos::null));

  MeshLogicalAudit audit(mesh, std::cout);
  CHECK(!audit.Verify());
}
