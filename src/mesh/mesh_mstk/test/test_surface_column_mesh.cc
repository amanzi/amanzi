/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:

*/

//!

// -------------------------------------------------------------
/**
 * @file   test_column_mesh.cc
 * @author Rao V. Garimella
 * @date   Thu Mar 18, 2015
 *
 * @brief
 *
 *
 */
// -------------------------------------------------------------
// -------------------------------------------------------------

#include <UnitTest++.h>

#include <mpi.h>
#include <fstream>

#include "AmanziComm.hh"

#include "Geometry.hh"
#include "../Mesh_MSTK.hh"
#include "../MeshColumn.hh"
#include "../MeshSurfaceCell.hh"
#include "RegionBox.hh"
#include "RegionLabeledSet.hh"
#include "GeometricModel.hh"


TEST(SURFACE_COLUMN_MESH_3D)
{
  auto comm = Amanzi::getDefaultComm();
  const int nproc(comm->getSize());
  const int me(comm->getRank());


  int nx = 4, ny = 4, nz = 4;
  double lx = 4, ly = 4, lz = 4;
  int dx = 1.0, dy = 1.0, dz = 1.0;

  // create a geometric model with surface region
  Teuchos::RCP<Amanzi::AmanziGeometry::GeometricModel> gm =
    Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(3));

  Amanzi::AmanziGeometry::Point p0(3), p1(3);
  p0[2] = 3.999;
  p1[0] = 4.;
  p1[1] = 4.;
  p1[2] = 5.;
  auto r0 =
    Teuchos::rcp(new Amanzi::AmanziGeometry::RegionBox("surface", -1, p0, p1));
  gm->AddRegion(r0);

  Amanzi::AmanziGeometry::Point p2(2), p3(2);
  p3[0] = 4.;
  p3[1] = 4.;
  auto r1 = Teuchos::rcp(
    new Amanzi::AmanziGeometry::RegionBox("surface_domain", -1, p2, p3));
  gm->AddRegion(r1);

  // Create the mesh
  Teuchos::RCP<Amanzi::AmanziMesh::Mesh> mesh =
    Teuchos::rcp(new Amanzi::AmanziMesh::Mesh_MSTK(
      0.0, 0.0, 0.0, lx, ly, lz, nx, ny, nz, comm, gm));

  CHECK_EQUAL(1, mesh->build_columns());

  // Perturb the nodes above the base layer just a bit
  int nnodes = mesh->num_entities(Amanzi::AmanziMesh::NODE,
                                  Amanzi::AmanziMesh::Parallel_type::OWNED);

  for (int n = 0; n < nnodes; n++) {
    Amanzi::AmanziGeometry::Point xyz(3);
    mesh->node_get_coordinates(n, &xyz);
    xyz[2] += 0.005 * xyz[0] * xyz[1] * xyz[2];
    mesh->node_set_coordinates(n, xyz);
  }

  // Create a column mesh from one of the columns
  auto colmesh = Teuchos::rcp(new Amanzi::AmanziMesh::MeshColumn(mesh, 10));

  std::cout << "Column mesh created" << std::endl;

  // Extract the surface from this column
  Amanzi::AmanziMesh::MeshSurfaceCell col_surf(colmesh, "surface");

  std::cout << "Extracted" << std::endl;

  // -- check basic mesh structure
  CHECK_EQUAL(1,
              col_surf.num_entities(Amanzi::AmanziMesh::CELL,
                                    Amanzi::AmanziMesh::Parallel_type::OWNED));
  CHECK_EQUAL(4,
              col_surf.num_entities(Amanzi::AmanziMesh::FACE,
                                    Amanzi::AmanziMesh::Parallel_type::OWNED));
  CHECK_EQUAL(4,
              col_surf.num_entities(Amanzi::AmanziMesh::NODE,
                                    Amanzi::AmanziMesh::Parallel_type::OWNED));

  // -- check flattened
  Amanzi::AmanziGeometry::Point node;
  col_surf.node_get_coordinates(0, &node);
  CHECK_EQUAL(2, node.dim());

  // -- check sets
  Kokkos::View<Amanzi::AmanziMesh::Entity_ID*> cells_in_surf;
  col_surf.get_set_entities_and_vofs("surface",
                                     Amanzi::AmanziMesh::CELL,
                                     Amanzi::AmanziMesh::Parallel_type::OWNED,
                                     cells_in_surf,
                                     NULL);
  CHECK_EQUAL(1, cells_in_surf.extent(0));
  CHECK_EQUAL(0, cells_in_surf(0));

  Kokkos::View<Amanzi::AmanziMesh::Entity_ID*> cells_in_surf_2D;
  col_surf.get_set_entities_and_vofs("surface_domain",
                                     Amanzi::AmanziMesh::CELL,
                                     Amanzi::AmanziMesh::Parallel_type::OWNED,
                                     cells_in_surf_2D,
                                     NULL);
  CHECK_EQUAL(1, cells_in_surf_2D.extent(0));
  CHECK_EQUAL(0, cells_in_surf_2D(0));

  // -- check volumes
  CHECK_CLOSE(1.0, col_surf.cell_volume(0, false), 1.e-9);
  CHECK_CLOSE(1.0, col_surf.face_area(3), 1.e-9);
}

TEST(SURFACE_COLUMN_MESH_3D_UNSTRUCTURED)
{
  auto comm = Amanzi::getDefaultComm();
  const int nproc(comm->getSize());
  const int me(comm->getRank());

  // create a geometric model with surface region
  Teuchos::RCP<Amanzi::AmanziGeometry::GeometricModel> gm =
    Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(3));

  Teuchos::RCP<Amanzi::AmanziGeometry::RegionLabeledSet> r0 =
    Teuchos::rcp(new Amanzi::AmanziGeometry::RegionLabeledSet(
      "surface", -1, "FACE", "test/slab-0.05-5x4x25.exo", "Exodus II", "1"));
  gm->AddRegion(r0);

  Amanzi::AmanziGeometry::Point p2(2), p3(2);
  p3[0] = 4.;
  p3[1] = 4.;
  auto r1 = Teuchos::rcp(
    new Amanzi::AmanziGeometry::RegionBox("surface_domain", -1, p2, p3));
  gm->AddRegion(r1);

  // Create the mesh
  Teuchos::RCP<Amanzi::AmanziMesh::Mesh> mesh = Teuchos::rcp(
    new Amanzi::AmanziMesh::Mesh_MSTK("test/slab-0.05-5x4x25.exo", comm, gm));

  CHECK_EQUAL(20,
              mesh->get_set_size("surface",
                                 Amanzi::AmanziMesh::FACE,
                                 Amanzi::AmanziMesh::Parallel_type::ALL));

  // Build columns in the mesh
  CHECK_EQUAL(1, mesh->build_columns());

  // Create a column mesh from one of the columns
  auto colmesh = Teuchos::rcp(new Amanzi::AmanziMesh::MeshColumn(mesh, 10));
  CHECK_EQUAL(1,
              colmesh->get_set_size("surface",
                                    Amanzi::AmanziMesh::FACE,
                                    Amanzi::AmanziMesh::Parallel_type::ALL));

  // Extract the surface from this column
  Amanzi::AmanziMesh::MeshSurfaceCell col_surf(colmesh, "surface");

  CHECK_EQUAL(1,
              col_surf.num_entities(Amanzi::AmanziMesh::CELL,
                                    Amanzi::AmanziMesh::Parallel_type::OWNED));
  CHECK_EQUAL(4,
              col_surf.num_entities(Amanzi::AmanziMesh::FACE,
                                    Amanzi::AmanziMesh::Parallel_type::OWNED));
  CHECK_EQUAL(4,
              col_surf.num_entities(Amanzi::AmanziMesh::NODE,
                                    Amanzi::AmanziMesh::Parallel_type::OWNED));

  Kokkos::View<Amanzi::AmanziMesh::Entity_ID*> cells_in_surf;
  col_surf.get_set_entities_and_vofs("surface",
                                     Amanzi::AmanziMesh::CELL,
                                     Amanzi::AmanziMesh::Parallel_type::OWNED,
                                     cells_in_surf,
                                     NULL);
  CHECK_EQUAL(1, cells_in_surf.extent(0));
  CHECK_EQUAL(0, cells_in_surf(0));

  Kokkos::View<Amanzi::AmanziMesh::Entity_ID*> cells_in_surf_2D;
  col_surf.get_set_entities_and_vofs("surface_domain",
                                     Amanzi::AmanziMesh::CELL,
                                     Amanzi::AmanziMesh::Parallel_type::OWNED,
                                     cells_in_surf_2D,
                                     NULL);
  CHECK_EQUAL(1, cells_in_surf.extent(0));
  CHECK_EQUAL(0, cells_in_surf(0));

  CHECK_CLOSE(6400.0, col_surf.cell_volume(0, false), 1.e-9);
  CHECK_CLOSE(80.0, col_surf.face_area(3), 1.e-9);
}
