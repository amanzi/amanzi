/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Rao Garimella
*/

/*
  Tests the MeshSurfaceCell object, which is a single surface cell for use as
  the "surface" mesh corresponding to a MeshColumn object.
*/

#include <UnitTest++.h>

#include <mpi.h>
#include <fstream>

#include "AmanziComm.hh"

#include "Geometry.hh"
#include "Mesh_MSTK.hh"
#include "MeshFrameworkColumn.hh"
#include "MeshSurfaceCell.hh"
#include "RegionBox.hh"
#include "RegionLabeledSet.hh"
#include "GeometricModel.hh"
#include "MeshCache.hh"

using namespace Amanzi;

TEST(SURFACE_COLUMN_MESH_3D)
{
  auto comm = getDefaultComm();

  int nx = 4, ny = 4, nz = 4;
  double lx = 4, ly = 4, lz = 4;

  // create a geometric model with surface region
  Teuchos::RCP<AmanziGeometry::GeometricModel> gm =
      Teuchos::rcp(new AmanziGeometry::GeometricModel(3));

  // Create the mesh
  Teuchos::RCP<AmanziMesh::MeshFramework> mesh_fw =
      Teuchos::rcp(new AmanziMesh::Mesh_MSTK(0.0,0.0,0.0,
              lx,ly,lz,nx,ny,nz, comm, gm));
  Teuchos::RCP<AmanziMesh::Mesh> mesh = Teuchos::rcp(new AmanziMesh::Mesh(mesh_fw));
  mesh->buildColumns();

  // Perturb the nodes above the base layer just a bit
  int nnodes = mesh->getNumEntities(AmanziMesh::Entity_kind::NODE,
          AmanziMesh::Parallel_type::OWNED);

  for (int n = 0; n < nnodes; n++) {
    AmanziGeometry::Point xyz(3);
    xyz = mesh->getNodeCoordinate(n);
    xyz[2] += 0.005*xyz[0]*xyz[1]*xyz[2];
    mesh->setNodeCoordinate(n,xyz);
  }

  // Create a column mesh from one of the columns
  Teuchos::RCP<AmanziMesh::MeshFramework> colmesh_ext =
    Teuchos::rcp(new AmanziMesh::Mesh_MSTK(mesh_fw,
          mesh->columns.cells_.getRow<MemSpace_type::HOST>(10), AmanziMesh::Entity_kind::CELL,
          false, getCommSelf(), gm, Teuchos::null));

  // Create the MeshColumn object
  Teuchos::RCP<AmanziMesh::MeshFramework> colmesh =
    Teuchos::rcp(new AmanziMesh::MeshFrameworkColumn(colmesh_ext, Teuchos::null));

  // Extract the surface from this column
  AmanziMesh::MeshSurfaceCell col_surf(colmesh);

  // -- check basic mesh structure
  CHECK_EQUAL(1, col_surf.getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_type::OWNED));
  CHECK_EQUAL(4, col_surf.getNumEntities(AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_type::OWNED));
  CHECK_EQUAL(4, col_surf.getNumEntities(AmanziMesh::Entity_kind::NODE, AmanziMesh::Parallel_type::OWNED));

  // -- check flattened
  AmanziGeometry::Point node;
  node = col_surf.getNodeCoordinate(0);
  CHECK_EQUAL(2, node.dim());

  // -- check volumes
  CHECK_CLOSE(1.0, col_surf.getCellVolume(0), 1.e-9);
  CHECK_CLOSE(1.0, col_surf.getFaceArea(3), 1.e-9);

  // -- check parent entity is the right face
  CHECK(colmesh->getFaceCentroid(col_surf.getEntityParent(AmanziMesh::Entity_kind::CELL, 0))[2] > 4.);

}
