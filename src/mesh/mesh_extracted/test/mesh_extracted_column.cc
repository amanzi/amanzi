/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Rao Garimella
*/

/*
  Note, this differs from src/mesh/test/mesh_column_extraction -- that test
  builds a 3D framework mesh from the extracted column, while this test builds
  a 1D MeshColumn object.
*/

#include <UnitTest++.h>

#include <fstream>

#include "AmanziComm.hh"
#include "Geometry.hh"
#include "MeshCache.hh"
#include "MeshFrameworkColumn.hh"
#include "RegionBox.hh"
#include "GeometricModel.hh"

#include "Mesh.hh"
#include "Mesh_MSTK.hh"

using namespace Amanzi;

TEST(COLUMN_MESH_3D)
{
  auto comm = getDefaultComm();

  int nx = 4, ny = 4, nz = 4;
  double lx = 4, ly = 4, lz = 4;

  // create a geometric model with regions
  Teuchos::RCP<AmanziGeometry::GeometricModel> gm =
      Teuchos::rcp(new AmanziGeometry::GeometricModel(3));

  AmanziGeometry::Point p0(3), p1(3);
  p0[2] = 2.5;
  p1[0] = 4.; p1[1] = 4.; p1[2] = 5.;
  auto r0 = Teuchos::rcp(new AmanziGeometry::RegionBox("myregion", 0, p0, p1));
  gm->AddRegion(r0);

  AmanziGeometry::Point p2(0.,0.,4.), p3(4.,4.,4.);
  auto r1 = Teuchos::rcp(new AmanziGeometry::RegionBox("surface", 0, p2, p3));
  gm->AddRegion(r1);

  // Create the mesh
  Teuchos::RCP<AmanziMesh::MeshFramework> mesh_fw =
    Teuchos::rcp(new AmanziMesh::Mesh_MSTK(0.0,0.0,0.0,
              lx,ly,lz, nx,ny,nz, comm,gm));
  auto mesh = Teuchos::rcp(new AmanziMesh::Mesh(mesh_fw, Teuchos::rcp(new AmanziMesh::MeshFrameworkAlgorithms()), Teuchos::null));
  mesh->buildColumns();

  // Perturb the nodes above the base layer just a bit
  int nnodes = mesh->getNumEntities(AmanziMesh::Entity_kind::NODE,
          AmanziMesh::Parallel_kind::OWNED);

  for (int n = 0; n < nnodes; n++) {
    AmanziGeometry::Point xyz(3);
    xyz = mesh->getNodeCoordinate(n);
    xyz[2] += 0.005*xyz[0]*xyz[1]*xyz[2];
    mesh->setNodeCoordinate(n,xyz);
  }

  // verify in-going topology
  CHECK_EQUAL(16,mesh->columns.num_columns_owned);
  CHECK_EQUAL(4, mesh->columns.cells_.size<MemSpace_kind::HOST>(10));
  CHECK_EQUAL(5, mesh->columns.faces_.size<MemSpace_kind::HOST>(10));

  // Create a column mesh from one of the columns
  auto colmesh_ext = Teuchos::rcp(new AmanziMesh::Mesh_MSTK(mesh->getMeshFramework(),
          mesh->columns.cells_.getRow<MemSpace_kind::HOST>(10), AmanziMesh::Entity_kind::CELL,
          false, getCommSelf(), gm, Teuchos::null));

  // Create the MeshColumn object
  Teuchos::RCP<AmanziMesh::MeshFramework> colmesh_fw =
    Teuchos::rcp(new AmanziMesh::MeshFrameworkColumn(colmesh_ext, Teuchos::null));
  AmanziMesh::Mesh colmesh(colmesh_fw, Teuchos::rcp(new AmanziMesh::MeshFrameworkColumnAlgorithms()), Teuchos::null);

  // Verify column mesh topology
  int ncells = colmesh.getNumEntities(AmanziMesh::Entity_kind::CELL,
          AmanziMesh::Parallel_kind::OWNED);
  CHECK_EQUAL(4,ncells);

  int nfaces = colmesh.getNumEntities(AmanziMesh::Entity_kind::FACE,
          AmanziMesh::Parallel_kind::OWNED);
  CHECK_EQUAL(5,nfaces);

  nnodes = colmesh.getNumEntities(AmanziMesh::Entity_kind::NODE,
          AmanziMesh::Parallel_kind::OWNED);
  CHECK_EQUAL(20,nnodes);

  for (int j = 0; j < ncells; j++) {
    AmanziMesh::cEntity_ID_View cfaces;
    AmanziMesh::cEntity_Direction_View cfdirs;
    colmesh.getCellFacesAndDirs(j,cfaces,&cfdirs);
    CHECK_EQUAL(2,cfaces.size());
  }

  for (int j = 0; j < nfaces; j++) {
    AmanziMesh::cEntity_ID_View fcells;
    colmesh.getFaceCells(j,AmanziMesh::Parallel_kind::OWNED,fcells);

    if (j == 0) {
      CHECK_EQUAL(1,fcells.size());
    } else if (j == nfaces-1) {
      CHECK_EQUAL(1,fcells.size());
    } else {
      CHECK_EQUAL(2,fcells.size());
    }
  }

  // Verify column mesh geometry
  // centroid of top face
  AmanziGeometry::Point fcenbase(3);
  fcenbase = colmesh.getFaceCentroid(0);
  CHECK_CLOSE(2.5, fcenbase[0], 1.e-10);
  CHECK_CLOSE(2.5, fcenbase[1], 1.e-10);

  // area of base face
  double fareabase = colmesh.getFaceArea(0);
  CHECK_CLOSE(1.0, fareabase, 1.e-10);

  // Make sure centroids of other faces are vertically stacked
  for (int j = 1; j < nfaces; j++) {
    AmanziGeometry::Point fcen(3);
    fcen = colmesh.getFaceCentroid(j);
    CHECK_EQUAL(fcenbase[0],fcen[0]);
    CHECK_EQUAL(fcenbase[1],fcen[1]);
  }

  // Make sure the normals of the faces are have only a Z component
  for (int j = 0; j < nfaces; j++) {
    auto normal = colmesh.getFaceNormal(j);
    CHECK_CLOSE(0.0,normal[0],1.e-10);
    CHECK_CLOSE(0.0,normal[1],1.e-10);
    CHECK_CLOSE(1.0,fabs(normal[2]),1.e-10);
  }

  // Make sure centroids of cells are stacked up
  // exactly above that of the base face and that
  // their z value is exactly between the z-values
  // of their corresponding bottom and top faces
  for (int i = 0; i < ncells; i++) {
    auto ccen = colmesh.getCellCentroid(i);
    auto fcen0 = colmesh.getFaceCentroid(i);
    auto fcen1 = colmesh.getFaceCentroid(i+1);

    CHECK_CLOSE(fcenbase[0],ccen[0],1.e-10);
    CHECK_CLOSE(fcenbase[1],ccen[1],1.e-10);
    CHECK_CLOSE((fcen0[2]+fcen1[2])/2.0,ccen[2],1.e-10);
  }

  // Verify the volume of cells is computed correctly in spite of
  // the topology being special (lateral faces of the cell are
  // omitted). The volume of each cell should be the area of the
  // base face times the distance between the centroids of the lower
  // face of the cell and the upper face of the cell
  for (int j = 0; j < ncells; j++) {
    const auto& cfaces = colmesh.getCellFaces(j);

    auto locen = colmesh.getFaceCentroid(cfaces[0]);
    auto hicen = colmesh.getFaceCentroid(cfaces[1]);

    double height = norm(hicen-locen);
    double expvolume = fareabase*height;
    CHECK_CLOSE(expvolume, colmesh.getCellVolume(j), 1.0e-08);
  }


  // verify that the regions have made it through
  AmanziMesh::Entity_ID_View myregion;
  myregion = colmesh.getSetEntities("myregion", AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::ALL);
  CHECK_EQUAL(2, myregion.size());
  CHECK(colmesh.getCellCentroid(myregion[0])[2] >= 2.5);
  CHECK(colmesh.getCellCentroid(myregion[1])[2] >= 2.5);
}


TEST(COLUMN_MESH_3D_FROM_SURFACE)
{
  auto comm = getDefaultComm();

  int nx = 4, ny = 4, nz = 4;
  double lx = 4, ly = 4, lz = 4;

  // create a geometric model with regions
  Teuchos::RCP<AmanziGeometry::GeometricModel> gm =
      Teuchos::rcp(new AmanziGeometry::GeometricModel(3));

  AmanziGeometry::Point p0(3), p1(3);
  p0[2] = 2.5;
  p1[0] = 4.; p1[1] = 4.; p1[2] = 5.;
  auto r0 = Teuchos::rcp(new AmanziGeometry::RegionBox("myregion", 0, p0, p1));
  gm->AddRegion(r0);

  AmanziGeometry::Point p2(0.,0.,4.), p3(4.,4.,4.);
  auto r1 = Teuchos::rcp(new AmanziGeometry::RegionBox("surface", 0, p2, p3));
  gm->AddRegion(r1);

  // Create the mesh
  Teuchos::RCP<AmanziMesh::MeshFramework> mesh_fw =
    Teuchos::rcp(new AmanziMesh::Mesh_MSTK(0.0,0.0,0.0,
            lx,ly,lz,nx,ny,nz,
            comm, gm));
  auto mesh = Teuchos::rcp(new AmanziMesh::Mesh(mesh_fw, Teuchos::rcp(new AmanziMesh::MeshFrameworkAlgorithms()), Teuchos::null));
  mesh->buildColumns({"surface"});

  int nnodes = mesh->getNumEntities(AmanziMesh::Entity_kind::NODE,
          AmanziMesh::Parallel_kind::OWNED);

  // verify in-going topology
  CHECK_EQUAL(16,mesh->columns.num_columns_owned);
  CHECK_EQUAL(4, mesh->columns.cells_.size<MemSpace_kind::HOST>(10));
  CHECK_EQUAL(5, mesh->columns.faces_.size<MemSpace_kind::HOST>(10));

  // Create a column mesh from one of the columns
  // Create a column mesh from one of the columns
  auto colmesh_ext = Teuchos::rcp(new AmanziMesh::Mesh_MSTK(mesh->getMeshFramework(),
          mesh->columns.cells_.getRow<MemSpace_kind::HOST>(10), AmanziMesh::Entity_kind::CELL,
          false, getCommSelf(), gm, Teuchos::null));
  // Create the MeshColumn object
  Teuchos::RCP<AmanziMesh::MeshFramework> colmesh_fw =
    Teuchos::rcp(new AmanziMesh::MeshFrameworkColumn(colmesh_ext, Teuchos::null));
  AmanziMesh::Mesh colmesh(colmesh_fw, Teuchos::rcp(new AmanziMesh::MeshFrameworkColumnAlgorithms()), Teuchos::null);

  // Verify column mesh topology
  int ncells = colmesh.getNumEntities(AmanziMesh::Entity_kind::CELL,
          AmanziMesh::Parallel_kind::OWNED);
  CHECK_EQUAL(4,ncells);

  int nfaces = colmesh.getNumEntities(AmanziMesh::Entity_kind::FACE,
          AmanziMesh::Parallel_kind::OWNED);
  CHECK_EQUAL(5,nfaces);

  nnodes = colmesh.getNumEntities(AmanziMesh::Entity_kind::NODE,
          AmanziMesh::Parallel_kind::OWNED);
  CHECK_EQUAL(20,nnodes);

  for (int j = 0; j < ncells; j++) {
    AmanziMesh::cEntity_ID_View cfaces;
    AmanziMesh::cEntity_Direction_View cfdirs;
    colmesh.getCellFacesAndDirs(j,cfaces,&cfdirs);

    CHECK_EQUAL(2,cfaces.size());
  }

  for (int j = 0; j < nfaces; j++) {
    AmanziMesh::cEntity_ID_View fcells;
    colmesh.getFaceCells(j,AmanziMesh::Parallel_kind::OWNED,fcells);

    if (j == 0) {
      CHECK_EQUAL(1,fcells.size());
    }
    else if (j == nfaces-1) {
      CHECK_EQUAL(1,fcells.size());
    }
    else {
      CHECK_EQUAL(2,fcells.size());
    }
  }

  // Verify column mesh geometry
  // centroid of base face
  AmanziGeometry::Point fcenbase(3);
  fcenbase = colmesh.getFaceCentroid(0);

  // area of base face
  double fareabase = colmesh.getFaceArea(0);

  // Make sure centroids of other faces are stacked up
  // exactly above that of the base face
  for (int j = 1; j < nfaces; j++) {
    AmanziGeometry::Point fcen(3);
    fcen = colmesh.getFaceCentroid(j);
    CHECK_CLOSE(fcenbase[0],fcen[0],1.e-10);
    CHECK_CLOSE(fcenbase[1],fcen[1],1.e-10);
  }

  // Make sure the normals of the faces are have only a Z component
  for (int j = 0; j < nfaces; j++) {
    AmanziGeometry::Point normal(3);
    normal = colmesh.getFaceNormal(j);
    CHECK_CLOSE(0.0,normal[0],1.e-10);
    CHECK_CLOSE(0.0,normal[1],1.e-10);
    CHECK_CLOSE(1.0,fabs(normal[2]),1.e-10);
  }

  // Make sure centroids of cells are stacked up
  // exactly above that of the base face and that
  // their z value is exactly between the z-values
  // of their corresponding bottom and top faces
  for (int i = 0; i < ncells; i++) {
    AmanziGeometry::Point ccen(3), fcen0(3), fcen1(3);

    ccen = colmesh.getCellCentroid(i);

    fcen0 = colmesh.getFaceCentroid(i);
    fcen1 = colmesh.getFaceCentroid(i+1);

    CHECK_CLOSE(fcenbase[0],ccen[0],1.e-10);
    CHECK_CLOSE(fcenbase[1],ccen[1],1.e-10);

    CHECK_CLOSE((fcen0[2]+fcen1[2])/2.0,ccen[2],1.e-10);
  }

  // Verify the volume of cells is computed correctly in spite of
  // the topology being special (lateral faces of the cell are
  // omitted). The volume of each cell should be the area of the
  // base face times the distance between the centroids of the lower
  // face of the cell and the upper face of the cell
  for (int j = 0; j < ncells; j++) {
    AmanziMesh::cEntity_ID_View cfaces;
    colmesh.getCellFaces(j,cfaces);

    AmanziGeometry::Point locen(3), hicen(3);
    locen = colmesh.getFaceCentroid(cfaces[0]);
    hicen = colmesh.getFaceCentroid(cfaces[1]);

    double height = norm(hicen-locen);

    double expvolume = fareabase*height;

    double volume = colmesh.getCellVolume(j);

    CHECK_CLOSE(expvolume,volume,1.0e-08);
  }

  // verify that the regions have made it through
  AmanziMesh::Entity_ID_View myregion;
  myregion = colmesh.getSetEntities("myregion", AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::ALL);
  CHECK_EQUAL(2, myregion.size());
  CHECK(colmesh.getCellCentroid(myregion[0])[2] >= 2.5);
  CHECK(colmesh.getCellCentroid(myregion[1])[2] >= 2.5);
}

