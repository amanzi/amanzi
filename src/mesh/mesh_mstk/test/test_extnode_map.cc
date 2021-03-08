#include <UnitTest++.h>

#include <fstream>

#include "../Mesh_MSTK.hh"
#include "MeshAudit.hh"

#include "Epetra_Map.h"
#include "AmanziComm.hh"
#include "Epetra_Import.h"
#include "Epetra_Vector.h"


TEST(MSTK_EXTNODE_MAP_4P)
{
  using namespace Amanzi;

  int nbnd;
  AmanziGeometry::Point xp(3);

  auto comm = Amanzi::getDefaultComm();
  int size = comm->NumProc();

  Teuchos::RCP<AmanziMesh::MeshFramework> mesh;
  if (size == 4) {
    mesh = Teuchos::rcp(new AmanziMesh::Mesh_MSTK("test/hex_3x3x3_sets.exo", comm));
    nbnd = 56;
  } else {
    mesh = Teuchos::rcp(new AmanziMesh::Mesh_MSTK(0.0, 0.0, 1.0, 1.0, 3, 3, comm));
    nbnd = 12;
  }

  int d = mesh->get_space_dimension();

  Epetra_Map node_map(mesh->node_map(true));
  Epetra_Map extnode_map(mesh->exterior_node_map(false));
  Epetra_Map extnode_map_wghost(mesh->exterior_node_map(true));

  CHECK(extnode_map.NumGlobalElements() == nbnd);

  for (int v = extnode_map.MinLID(); v <= extnode_map.MaxLID(); ++v) {
    int gid = extnode_map.GID(v);
    int v2 = node_map.LID(gid);

    CHECK_EQUAL(node_map.GID(v2), gid);
  }

  // Check if ghostmap contains only boundary faces
  for (int v = extnode_map_wghost.MinLID(); v <= extnode_map_wghost.MaxLID(); ++v) {
    int gid = extnode_map_wghost.GID(v);
    xp = mesh->getNodeCoordinate(node_map.LID(gid));
    CHECK(fabs(xp[0]) < 1e-6 || fabs(xp[0] - 1.0) < 1e-6 ||
          fabs(xp[1]) < 1e-6 || fabs(xp[1] - 1.0) < 1e-6 ||
          fabs(xp[d-1]) < 1e-6 || fabs(xp[d-1] - 1.0) < 1e-6);
  }
}

