#include <UnitTest++.h>

#include <fstream>

#include "MeshAlgorithms_impl.hh"
#include "MeshAlgorithms_decl.hh"
#include "../Mesh_MSTK.hh"
#include "MeshAudit.hh"

#include "AmanziComm.hh"

// Test edge functions in 2D

TEST(MSTK_EDGES_2D)
{
  auto comm = Amanzi::getDefaultComm();
  int size = comm->NumProc();
  CHECK_EQUAL(4,size);

  int DebugWait = 0;
  while (DebugWait);

  // Generate a 4x4 quad mesh distributed over four processors
  auto plist = Teuchos::rcp(new Teuchos::ParameterList("my mesh"));
  plist->set("request edges", true);
  Teuchos::RCP<Amanzi::AmanziMesh::MeshFramework>
    mesh(new Amanzi::AmanziMesh::Mesh_MSTK(0.0,0.0,2.0,1.0,4,4,comm,
            Teuchos::null, plist));

  // Check that we get the expected number of edges
  int ne_owned = mesh->getNumEntities(Amanzi::AmanziMesh::Entity_kind::EDGE,
          Amanzi::AmanziMesh::Parallel_type::OWNED);
  int ne_all = mesh->getNumEntities(Amanzi::AmanziMesh::Entity_kind::EDGE,
				  Amanzi::AmanziMesh::Parallel_type::ALL);
  CHECK(ne_owned <= ne_all);

  // This assumes a symmetric partitioning - not always the case with
  // ZOLTAN graph partitioning

  //  CHECK_EQUAL(24,ne_all);

  // In 2D, faces and edges are the same - so face global IDs and edge
  // global IDs for a cell must match

  int nc_owned = mesh->getNumEntities(Amanzi::AmanziMesh::Entity_kind::CELL,
          Amanzi::AmanziMesh::Parallel_type::OWNED);

  for (int c = 0; c < nc_owned; ++c) {
    Amanzi::AmanziMesh::Entity_ID_List cedges, cfaces, fedges;
    Amanzi::AmanziMesh::Entity_Direction_List cfdirs, fedirs;

    mesh->getCellEdges(c,cedges);
    mesh->getCellFacesAndDirs(c,cfaces,&cfdirs);

    bool pass = true;
    for (int e = 0; e < cedges.size(); ++e) {
      CHECK_EQUAL(mesh->getEntityGID(Amanzi::AmanziMesh::Entity_kind::EDGE, cedges[e]),
                  mesh->getEntityGID(Amanzi::AmanziMesh::Entity_kind::FACE, cfaces[e]));

      Amanzi::AmanziGeometry::Point evec(2), fnormal(2), ftangent(2);
      evec = mesh->getEdgeVector(cfaces[e]);
      fnormal = mesh->getFaceNormal(cfaces[e]);
      ftangent.set(-fnormal[1],fnormal[0]);
      pass &= (evec == ftangent);
      CHECK_EQUAL(evec[0],ftangent[0]);
      CHECK_EQUAL(evec[1],ftangent[1]);
      CHECK(evec == ftangent);

      // check directions is correct for star-convex
      auto cc = mesh->getCellCentroid(c);
      auto fc = mesh->getFaceCentroid(cfaces[e]);
      auto ec = mesh->getEdgeCentroid(cedges[e]);
      CHECK((fc - cc) * (cfdirs[e]*fnormal) > 0);
      pass &= ((fc - cc) * (cfdirs[e]*fnormal) > 0);

      // face and edge centroid are the same
      CHECK(fc == ec);
      pass &= (fc == ec);

      // normal with respect to cell
      int dir2;
      auto fnormal2 = mesh->getFaceNormal(cfaces[e], c, &dir2);
      CHECK_EQUAL(cfdirs[e], dir2);
      pass &= (cfdirs[e] == dir2);
      CHECK(fnormal2 == (cfdirs[e] * fnormal));
      pass &= (fnormal2 == (cfdirs[e] * fnormal));
    }
    if (!pass) {
      Amanzi::AmanziMesh::MeshAlgorithms::debugCell(*mesh, c);
    }

    for (int f = 0; f < cfaces.size(); ++f) {
      mesh->getFaceEdgesAndDirs(cfaces[f],fedges,&fedirs);
      CHECK_EQUAL(1,fedges.size()); // face is same as edge in 2D
      CHECK_EQUAL(1,fedirs[0]); // direction is always 1
    }

    // owing to how we constructed the mesh, the length of horizontal edges
    // should be 0.5 and vertical edges 0.25
    for (int e = 0; e < ne_owned; ++e) {
      Amanzi::AmanziGeometry::Point evec(2);
      double elen;

      evec = mesh->getEdgeVector(e);
      elen = mesh->getEdgeLength(e);
      if (evec[1] == 0.0) {
        CHECK_EQUAL(0.5,elen);
        CHECK_EQUAL(elen,norm(evec));
      }
      else if (evec[0] == 0.0) {
        CHECK_EQUAL(0.25,elen);
        CHECK_EQUAL(elen,norm(evec));
      }
    }
  }
}


// Test edge functions in 3D
TEST(MSTK_EDGES_3D)
{
  auto comm = Amanzi::getDefaultComm();
  int size = comm->NumProc();
  CHECK_EQUAL(4,size);

  //  if (rank == 0) {
  int DebugWait = 0;
  while (DebugWait);
  //  }

  // Generate a 4x4x4 quad mesh distributed over four processors
  auto plist = Teuchos::rcp(new Teuchos::ParameterList("my mesh"));
  plist->set("request edges", true);
  Teuchos::RCP<Amanzi::AmanziMesh::MeshFramework>
    mesh(new Amanzi::AmanziMesh::Mesh_MSTK(0.0,0.0,0.0,2.0,1.0,4.0,4,4,4,
            comm,Teuchos::null,plist));

  // How many owned and used edges are there?
  int ne_owned = mesh->getNumEntities(Amanzi::AmanziMesh::Entity_kind::EDGE,
          Amanzi::AmanziMesh::Parallel_type::OWNED);
  int ne_all = mesh->getNumEntities(Amanzi::AmanziMesh::Entity_kind::EDGE,
          Amanzi::AmanziMesh::Parallel_type::ALL);

  // Check that we got a non-zero number
  CHECK(ne_owned != 0);
  CHECK(ne_all != 0);

  // Go through the cells and retrieve their edges to make sure it
  // works correctly. Also, get the faces of the cells and the edges
  // of these faces and do additional checks
  int nc_owned = mesh->getNumEntities(Amanzi::AmanziMesh::Entity_kind::CELL,
          Amanzi::AmanziMesh::Parallel_type::OWNED);

  for (int c = 0; c < nc_owned; ++c) {
    Amanzi::AmanziMesh::Entity_ID_List cedges, cfaces, fedges;
    std::vector<int> cfdirs, fedirs;

    mesh->getCellEdges(c,cedges);
    mesh->getCellFacesAndDirs(c,cfaces,&cfdirs);

    for (int f = 0; f < cfaces.size(); ++f) {
      mesh->getFaceEdgesAndDirs(cfaces[f],fedges,&fedirs);
    }
  }

  // owing to how we constructed the mesh, the length of x-direction
  // should be 0.5, y-direction edges should
  // 0.25 and z-direction edges should be 1.0
  for (int e = 0; e < ne_owned; ++e) {
    Amanzi::AmanziGeometry::Point evec(2);
    double elen;

    evec = mesh->getEdgeVector(e);
    elen = mesh->getEdgeLength(e);
    if (evec[0] != 0.0 && evec[1] == 0.0) {
      CHECK_EQUAL(0.5,elen);
      CHECK_EQUAL(elen,norm(evec));
    }
    else if (evec[0] == 0.0 && evec[1] != 0.0) {
      CHECK_EQUAL(0.25,elen);
      CHECK_EQUAL(elen,norm(evec));
    }
    else  {
      CHECK_EQUAL(1.0,elen);
      CHECK_EQUAL(elen,norm(evec));
    }
  }
}

