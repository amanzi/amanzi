/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon
*/
#include "UnitTest++.h"

#include "AmanziComm.hh"
#include "RegionBox.hh"
#include "MeshFactory.hh"
#include "Patch.hh"

using namespace Amanzi;
using namespace Amanzi::AmanziGeometry;
using namespace Amanzi::AmanziMesh;

struct test_mesh {
  Comm_ptr_type comm;
  Teuchos::RCP<GeometricModel> gm;
  Teuchos::RCP<Mesh> mesh;

  test_mesh()
  {
    comm = getDefaultComm();
    gm = Teuchos::rcp(new GeometricModel(2));
    gm->AddRegion(Teuchos::rcp(new RegionBox("box1", 0, Point(0.,0.),
            Point(2.0, 4.0))));
    gm->AddRegion(Teuchos::rcp(new RegionBox("box2", 1, Point(2.,0.),
            Point(4.0, 4.0))));
    
    MeshFactory meshfactory(comm, gm);
    mesh = meshfactory.create(0.0, 0.0, 4.0, 4.0, 2, 2);
  }
  ~test_mesh() {}
};


SUITE(PATCH) {

TEST_FIXTURE(test_mesh, PATCH_CREATE) {
  PatchSpace space{mesh, false, "box1", CELL, 2, -1};
  CHECK_EQUAL(2, space.size());
}

TEST_FIXTURE(test_mesh, MULTIPATCH_CREATE) {
  MultiPatchSpace space;
  space.mesh = mesh;
  space.ghosted = false;

  space.AddPatch("box1", CELL, 2);
  space.AddPatch("box2", CELL, 2);
  auto mp = space.Create();

  auto p1 = (*mp)[0];
  CHECK_EQUAL(2, p1.data.extent(0));
  CHECK_EQUAL(2, p1.data.extent(1));

  auto p2 = (*mp)[1];
  CHECK_EQUAL(2, p2.data.extent(0));
  CHECK_EQUAL(2, p2.data.extent(1));

}  

}
