/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

#pragma once

#include "MeshFrameworkTraits.hh"
#include "MeshAudit_impl.hh"

using namespace Amanzi;
using namespace Amanzi::AmanziMesh;

//
// Helper functions for testing geometry
//

//
// This is a helper function -- simply runs MeshAudit
//
template <class MeshAudit_type, class Mesh_type>
bool
testMeshAuditHost(const Teuchos::RCP<const Mesh_type>& mesh)
{
  // run MeshAudit
  // note on host, passing POINTER to the object
  Mesh_type const * mc = mesh.get();

  MeshAudit_type audit(mesh, mc);
  int status = audit.Verify();
  CHECK_EQUAL(0, status);
  return status;
}

template <class MeshAudit_type>
bool
testMeshAuditDevice(const Teuchos::RCP<const Mesh>& mesh)
{
  // run MeshAudit
  // note on device, passing the cache
  MeshAudit_type audit(mesh, mesh->getCache());
  int status = audit.Verify();
  CHECK_EQUAL(0, status);
  return status;
}


//
// Sums a scalar over all processes and compares the result to exp
//
// Used to check that partial counts add up to global count in parallel tests.
//
template <typename T>
void
CHECK_CLOSE_SUMALL(T exp, T contrib, const Amanzi::Comm_type& comm, T tol = 0)
{
  // MPI-based CHECK_CLOSE using SumAll
  T global;
  Teuchos::reduceAll(comm, Teuchos::REDUCE_SUM, 1, &contrib, &global);
  CHECK_CLOSE(exp, global, tol);
  if (std::abs(exp - global) > tol) {
    // for debugging
    // double diff = std::abs(exp - global);
  }
}


template<MemSpace_kind MEM>
const auto&
getMeshForKokkos(const Mesh& mesh) {
  if constexpr(MEM == MemSpace_kind::HOST) {
    return mesh;
  } else {
    return mesh.getCache();
  }
}

//
// Tests geometry given expected values
//
template <class MeshOrCache_type>
bool
testMeshGeometry_(const Mesh& mesh,
                  const MeshOrCache_type& m,
                  const typename std::remove_pointer<MeshOrCache_type>::type::cPoint_View& exp_cell_centroids,
                  const typename std::remove_pointer<MeshOrCache_type>::type::cDouble_View& exp_cell_volumes,
                  const typename std::remove_pointer<MeshOrCache_type>::type::cPoint_View& exp_face_centroids,
                  const typename std::remove_pointer<MeshOrCache_type>::type::cDouble_View& exp_face_areas,
                  const typename std::remove_pointer<MeshOrCache_type>::type::cPoint_View& exp_face_normals,
                  const typename std::remove_pointer<MeshOrCache_type>::type::cPoint_View& exp_node_coordinates)
{
  bool test_error = false;
  using MeshOrCache_noptr_type = typename std::remove_pointer<MeshOrCache_type>::type;

  // test cell-based quantities
  {
    std::cout << "Checking cell geometry ..." << std::endl;
    int ncells = mesh.getNumEntities(Entity_kind::CELL, Parallel_kind::OWNED);
    typename MeshOrCache_noptr_type::Entity_ID_View test1("testMeshGeometry cells 1", ncells);
    typename MeshOrCache_noptr_type::Entity_ID_View test2("testMeshGeometry cells 2", ncells);
    typename MeshOrCache_noptr_type::Entity_ID_View test3("testMeshGeometry cells 3", ncells);

    using Range_type =
      Kokkos::RangePolicy<LO, typename MeshOrCache_noptr_type::cEntity_ID_View::execution_space>;

    Kokkos::parallel_for(
      "GeometryHarnesses::testMeshGeometry",
      Range_type{ 0, ncells },
      KOKKOS_LAMBDA(const LO c) {
        test1(c) = 1;
        auto centroid = AmanziMesh::Impl::get(m).getCellCentroid(c);

        // Search for a cell with the same centroid in the expected list of centroids
        int j = 0;
        for (; j != exp_cell_centroids.size(); ++j) {
          auto diff = exp_cell_centroids[j] - centroid;
          if (AmanziGeometry::norm(diff) < 1.0e-10) {
            test1(c) = 0;
            break;
          }
        }

        // check cell volume matches
        if (fabs(exp_cell_volumes[j] - AmanziMesh::Impl::get(m).getCellVolume(c)) > 1e-10) { test2(c) = 1; }

        // check that the outward normals sum to 0
        auto cfaces = AmanziMesh::Impl::get(m).getCellFaces(c);
        AmanziGeometry::Point normal_sum(AmanziMesh::Impl::get(m).getSpaceDimension());
        for (int k = 0; k < cfaces.size(); k++) {
          auto normal = AmanziMesh::Impl::get(m).getFaceNormal(cfaces[k], c);
          normal_sum = normal_sum + normal;
        }
        double val = AmanziGeometry::norm(normal_sum);
        if (val > 1.e-10) {
          test3(c) = 1;
          // for (int k = 0; k < cfaces.size(); k++) {
          //   auto normal = AmanziMesh::Impl::get(m).getFaceNormal(cfaces[k], c);
          //   std::cout << "fail cell " << i << " normal (" << k << ") " << cfaces[k] << " = " << normal
          //             << std::endl;
          // }
        }
      });

    bool error = AmanziMesh::Impl::checkErrorList(
      test1, "GeometryHarnesses::testMeshGeometry cannot find matching cell centroid", std::cerr);
    error |= AmanziMesh::Impl::checkErrorList(
      test2, "GeometryHarnesses::testMeshGeometry bad cell volume", std::cerr);
    error |= AmanziMesh::Impl::checkErrorList(
      test3, "GeometryHarnesses::testMeshGeometry cell normals don't sum to 0", std::cerr);
    test_error = AmanziMesh::Impl::globalAny(*mesh.getComm(), error);
  }

  // test face-based quantities
  {
    std::cout << "Checking face geometry ..." << std::endl;
    int nfaces = mesh.getNumEntities(Entity_kind::FACE, Parallel_kind::OWNED);
    typename MeshOrCache_noptr_type::Entity_ID_View test1("testMeshGeometry faces 1", nfaces);
    typename MeshOrCache_noptr_type::Entity_ID_View test2("testMeshGeometry faces 2", nfaces);
    typename MeshOrCache_noptr_type::Entity_ID_View test3("testMeshGeometry faces 3", nfaces);
    typename MeshOrCache_noptr_type::Entity_ID_View test4("testMeshGeometry faces 4", nfaces);

    using Range_type =
      Kokkos::RangePolicy<LO, typename MeshOrCache_noptr_type::cEntity_ID_View::execution_space>;
    Kokkos::parallel_for(
      "GeometryHarnesses::testMeshGeometry", Range_type{ 0, nfaces }, KOKKOS_LAMBDA(const LO f) {
        test1(f) = 1;
        AmanziGeometry::Point centroid = AmanziMesh::Impl::get(m).getFaceCentroid(f);

        int j = 0;
        for (; j < exp_face_centroids.size(); ++j) {
          auto diff = exp_face_centroids[j] - centroid;
          if (AmanziGeometry::norm(diff) < 1.0e-10) {
            test1(f) = 0;
            break;
          }
        }

        // compare face areas
        if (fabs(exp_face_areas[j] - AmanziMesh::Impl::get(m).getFaceArea(f)) > 1.e-10) { test2(f) = 1; }

        // Natural normal is well-posed
        AmanziGeometry::Point natural_normal = AmanziMesh::Impl::get(m).getFaceNormal(f);

        // Check the normal with respect to each connected cell is given as the
        // natural times the orientation.
        auto cellids = AmanziMesh::Impl::get(m).getFaceCells(f);

        for (int k = 0; k < cellids.size(); k++) {
          int orientation = 0;
          auto normal_wrt_cell = AmanziMesh::Impl::get(m).getFaceNormal(f, cellids[k], &orientation);
          if (natural_normal * orientation != normal_wrt_cell) {
            test3(f) = 1;
            //   std::cout << "Fail face: " << i << " wrt cell " << cellids[k] << std::endl
            //             << "  nat_normal = " << natural_normal << std::endl
            //             << "  normal_wrt = " << normal_wrt_cell << std::endl
            //             << "  orientation = " << orientation << std::endl;
          }

          // check the cell's outward normal is indeed outward (assumes star-convex)
          AmanziGeometry::Point cellcentroid = AmanziMesh::Impl::get(m).getCellCentroid(cellids[k]);
          AmanziGeometry::Point facecentroid = AmanziMesh::Impl::get(m).getFaceCentroid(f);
          AmanziGeometry::Point outvec = facecentroid - cellcentroid;

          double dp = outvec * normal_wrt_cell;
          dp /= (AmanziGeometry::norm(outvec) * AmanziGeometry::norm(normal_wrt_cell));
          if (fabs(dp - 1) > 1.e-10) { test4(f) = 1; }
        }
      });


    bool error = AmanziMesh::Impl::checkErrorList(
      test1, "GeometryHarnesses::testMeshGeometry cannot find matching face centroid", std::cerr);
    error |= AmanziMesh::Impl::checkErrorList(
      test2, "GeometryHarnesses::testMeshGeometry bad face area", std::cerr);
    error |= AmanziMesh::Impl::checkErrorList(
      test3,
      "GeometryHarnesses::testMeshGeometry face normals wrt cell not consistent with orientation",
      std::cerr);
    error |= AmanziMesh::Impl::checkErrorList(
      test4, "GeometryHarnesses::testMeshGeometry face normal wrt cell not outward", std::cerr);
    test_error |= AmanziMesh::Impl::globalAny(*mesh.getComm(), error);
  }

  // test the node-based quantities
  {
    std::cout << "Checking node geometry ..." << std::endl;
    int nnodes = mesh.getNumEntities(Entity_kind::NODE, Parallel_kind::OWNED);
    typename MeshOrCache_noptr_type::Entity_ID_View test1("testMeshGeometry nodes 1", nnodes);

    using Range_type =
      Kokkos::RangePolicy<LO, typename MeshOrCache_noptr_type::cEntity_ID_View::execution_space>;
    Kokkos::parallel_for(
      "GeometryHarnesses::testMeshGeometry", Range_type{ 0, nnodes }, KOKKOS_LAMBDA(const LO n) {
        test1(n) = 1;

        AmanziGeometry::Point centroid = AmanziMesh::Impl::get(m).getNodeCoordinate(n);
        int j = 0;
        for (; j < exp_node_coordinates.size(); ++j) {
          auto diff = exp_node_coordinates[j] - centroid;
          if (AmanziGeometry::norm(diff) < 1.0e-10) {
            test1(n) = 0;
            break;
          }
        }
      });

    bool error = AmanziMesh::Impl::checkErrorList(
      test1, "GeometryHarnesses::testMeshGeometry cannot find matching node coordinate", std::cerr);
    test_error |= AmanziMesh::Impl::globalAny(*mesh.getComm(), error);
  }
  return test_error;
}


template <class MeshOrCache_type>
bool
testMeshGeometry(const Mesh& mesh,
                 const MeshOrCache_type& mesh_or_cache,
                 const Point_List& exp_cell_centroids,
                 const Double_List& exp_cell_volumes,
                 const Point_List& exp_face_centroids,
                 const Double_List& exp_face_areas,
                 const Point_List& exp_face_normals,
                 const Point_List& exp_node_coordinates)
{
  // run the testing
  return testMeshGeometry_(mesh,
                           mesh_or_cache,
                           toNonOwningView(exp_cell_centroids),
                           toNonOwningView(exp_cell_volumes),
                           toNonOwningView(exp_face_centroids),
                           toNonOwningView(exp_face_areas),
                           toNonOwningView(exp_face_normals),
                           toNonOwningView(exp_node_coordinates));
}

template <>
inline bool
testMeshGeometry<AmanziMesh::MeshCache>(const AmanziMesh::Mesh& mesh,
                                   const AmanziMesh::MeshCache& mesh_or_cache,
                                   const Point_List& exp_cell_centroids,
                                   const Double_List& exp_cell_volumes,
                                   const Point_List& exp_face_centroids,
                                   const Double_List& exp_face_areas,
                                   const Point_List& exp_face_normals,
                                   const Point_List& exp_node_coordinates)
{
  AmanziMesh::MeshCache::Point_View exp_cell_centroids_d("expected cell centroids",
                                                    exp_cell_centroids.size());
  Kokkos::deep_copy(exp_cell_centroids_d, toNonOwningView(exp_cell_centroids));
  AmanziMesh::MeshCache::Double_View exp_cell_volumes_d("expected cell volumes",
                                                   exp_cell_volumes.size());
  Kokkos::deep_copy(exp_cell_volumes_d, toNonOwningView(exp_cell_volumes));
  AmanziMesh::MeshCache::Point_View exp_face_centroids_d("expected face centroids",
                                                    exp_face_centroids.size());
  Kokkos::deep_copy(exp_face_centroids_d, toNonOwningView(exp_face_centroids));
  AmanziMesh::MeshCache::Double_View exp_face_areas_d("expected face areas", exp_face_areas.size());
  Kokkos::deep_copy(exp_face_areas_d, toNonOwningView(exp_face_areas));
  AmanziMesh::MeshCache::Point_View exp_face_normals_d("expected face normals", exp_face_normals.size());
  Kokkos::deep_copy(exp_face_normals_d, toNonOwningView(exp_face_normals));
  AmanziMesh::MeshCache::Point_View exp_node_coordinates_d("expected node coordinates",
                                                      exp_node_coordinates.size());
  Kokkos::deep_copy(exp_node_coordinates_d, toNonOwningView(exp_node_coordinates));

  return testMeshGeometry_(mesh,
                           mesh_or_cache,
                           exp_cell_centroids_d,
                           exp_cell_volumes_d,
                           exp_face_centroids_d,
                           exp_face_areas_d,
                           exp_face_normals_d,
                           exp_node_coordinates_d);
}


template <class Mesh_type>
void
testGeometryQuadBasics(const Mesh_type& mesh,
                       int nx, int ny)
{
  // test the basic dimensionality
  CHECK_EQUAL(2, mesh.getSpaceDimension());
  int ncells = mesh.getNumEntities(Entity_kind::CELL, Parallel_kind::OWNED);
  int ncells_test = nx * ny;
  CHECK_CLOSE_SUMALL(ncells_test, ncells, *mesh.getComm());

  int nfaces = mesh.getNumEntities(Entity_kind::FACE, Parallel_kind::OWNED);
  int nfaces_test = ny * (nx + 1) + nx * (ny + 1);
  CHECK_CLOSE_SUMALL(nfaces_test, nfaces, *mesh.getComm());

  int nnodes = mesh.getNumEntities(Entity_kind::NODE, Parallel_kind::OWNED);
  int nnodes_test = (ny + 1) * (nx + 1);
  CHECK_CLOSE_SUMALL(nnodes_test, nnodes, *mesh.getComm());
}

//
// Form the expected values and call testGeometry for a 2D box
//
template <class MeshOrCache_type>
void
testGeometryQuad(const Mesh& mesh,
                 const MeshOrCache_type& mesh_or_cache,
                 int nx, int ny)
{
  testGeometryQuadBasics(mesh, nx, ny);

  int ncells_test = nx * ny;

  // construct expected cell volumes, centroids
  Double_List exp_cell_volumes(ncells_test, 1. / nx * 1. / ny);
  Point_List exp_cell_centroids;
  for (int i = 0; i != nx; ++i) {
    for (int j = 0; j != ny; ++j) {
      exp_cell_centroids.emplace_back(AmanziGeometry::Point((i + .5) / nx, (j + .5) / ny));
    }
  }

  // construct expected face centroids, areas, and normals
  Point_List exp_face_centroids;
  Double_List exp_face_areas;
  Point_List exp_face_normals;
  for (int i = 0; i != nx + 1; ++i) {
    for (int j = 0; j != ny; ++j) {
      exp_face_centroids.emplace_back(AmanziGeometry::Point(((double)i) / nx, (j + .5) / ny));
      exp_face_areas.emplace_back(1.0 / ny);
      exp_face_normals.emplace_back(AmanziGeometry::Point(1.0, 0.0));
    }
  }
  for (int i = 0; i != nx; ++i) {
    for (int j = 0; j != ny + 1; ++j) {
      exp_face_centroids.emplace_back(AmanziGeometry::Point((i + .5) / nx, ((double)j) / ny));
      exp_face_areas.emplace_back(1.0 / nx);
      exp_face_normals.emplace_back(AmanziGeometry::Point(0.0, 1.0));
    }
  }

  // construct expected nodal locations
  Point_List exp_node_coordinates;
  for (int i = 0; i != nx + 1; ++i) {
    for (int j = 0; j != ny + 1; ++j) {
      exp_node_coordinates.emplace_back(AmanziGeometry::Point(((double)i) / nx, ((double)j) / ny));
    }
  }

  bool result = testMeshGeometry(mesh,
                                 mesh_or_cache,
                                 exp_cell_centroids,
                                 exp_cell_volumes,
                                 exp_face_centroids,
                                 exp_face_areas,
                                 exp_face_normals,
                                 exp_node_coordinates);
  CHECK(!result);
}


template <class Mesh_type>
void
testGeometryCubeBasics(const Mesh_type& mesh, int nx, int ny, int nz)
{
  // test the basic dimensionality
  CHECK_EQUAL(3, mesh.getSpaceDimension());
  int ncells = mesh.getNumEntities(Entity_kind::CELL, Parallel_kind::OWNED);
  int ncells_test = nx * ny * nz;
  CHECK_CLOSE_SUMALL(ncells_test, ncells, *mesh.getComm());
  int nfaces = mesh.getNumEntities(Entity_kind::FACE, Parallel_kind::OWNED);
  int nfaces_test = nx * ny * (nz + 1) + nx * (ny + 1) * nz + (nx + 1) * ny * nz;
  CHECK_CLOSE_SUMALL(nfaces_test, nfaces, *mesh.getComm());
  int nnodes = mesh.getNumEntities(Entity_kind::NODE, Parallel_kind::OWNED);
  int nnodes_test = (nx + 1) * (ny + 1) * (nz + 1);
  CHECK_CLOSE_SUMALL(nnodes_test, nnodes, *mesh.getComm());
}


//
// Form the expected values and call testGeometry for a 3D cube
//
template <class Mesh_type>
void
testGeometryCube(const Mesh& mesh, const Mesh_type& mesh_or_cache,
                 int nx, int ny, int nz)
{
  testGeometryCubeBasics(mesh, nx, ny, nz);

  int ncells_test = nx * ny * nz;
  // construct expected cell volumes, centroids
  Double_List exp_cell_volumes(ncells_test, 1. / nx * 1. / ny * 1. / nz);
  Point_List exp_cell_centroids;
  for (int i = 0; i != nx; ++i) {
    for (int j = 0; j != ny; ++j) {
      for (int k = 0; k != nz; ++k) {
        exp_cell_centroids.emplace_back(
          AmanziGeometry::Point((i + .5) / nx, (j + .5) / ny, (k + .5) / nz));
      }
    }
  }

  // construct expected face centroids, areas, and normals
  Point_List exp_face_centroids;
  Double_List exp_face_areas;
  Point_List exp_face_normals;
  for (int i = 0; i != (nx + 1); ++i) {
    for (int j = 0; j != ny; ++j) {
      for (int k = 0; k != nz; ++k) {
        exp_face_centroids.emplace_back(
          AmanziGeometry::Point(((double)i) / nx, (j + .5) / ny, (k + .5) / nz));
        exp_face_areas.emplace_back(1.0 / ny * 1.0 / nz);
        exp_face_normals.emplace_back(AmanziGeometry::Point(1.0, 0.0, 0.0));
      }
    }
  }
  for (int i = 0; i != nx; ++i) {
    for (int j = 0; j != (ny + 1); ++j) {
      for (int k = 0; k != nz; ++k) {
        exp_face_centroids.emplace_back(
          AmanziGeometry::Point((i + .5) / nx, ((double)j) / ny, (k + .5) / nz));
        exp_face_areas.emplace_back(1.0 / nx * 1.0 / nz);
        exp_face_normals.emplace_back(AmanziGeometry::Point(0.0, 1.0, 0.0));
      }
    }
  }
  for (int i = 0; i != nx; ++i) {
    for (int j = 0; j != ny; ++j) {
      for (int k = 0; k != (nz + 1); ++k) {
        exp_face_centroids.emplace_back(
          AmanziGeometry::Point((i + .5) / nx, (j + .5) / ny, ((double)k) / nz));
        exp_face_areas.emplace_back(1.0 / nx * 1.0 / ny);
        exp_face_normals.emplace_back(AmanziGeometry::Point(0.0, 0.0, 1.0));
      }
    }
  }

  // construct expected nodal locations
  Point_List exp_node_coordinates;
  for (int i = 0; i != nx + 1; ++i) {
    for (int j = 0; j != ny + 1; ++j) {
      for (int k = 0; k != nz + 1; ++k) {
        exp_node_coordinates.emplace_back(
          AmanziGeometry::Point(((double)i) / nx, ((double)j) / ny, ((double)k) / nz));
      }
    }
  }

  // run the testing
  bool res = testMeshGeometry(mesh,
                              mesh_or_cache,
                              exp_cell_centroids,
                              exp_cell_volumes,
                              exp_face_centroids,
                              exp_face_areas,
                              exp_face_normals,
                              exp_node_coordinates);
  CHECK(!res);
}


//
// Test the exterior face and node maps for a unit (2D quad or 3D hex) box.
//
// consistency checks on maps for cells/faces/nodes are done in MeshAudit.
// Here are are simply checking that boundary faces and nodes are actually on
// the boundary.
//
// Note this is valid on 2D quads too with default value nz = -1
//
template <class MeshOrCache_type>
void
testExteriorMapsUnitBox(const Mesh& mesh, const MeshOrCache_type& m, int nx, int ny, int nz = -1)
{
  // check faces are on the boundary
  int nbfaces = mesh.getMap(AmanziMesh::Entity_kind::BOUNDARY_FACE, false)->getGlobalNumElements();
  int nbfaces_test;
  if (nz < 0) {
    nbfaces_test = 2 * nx + 2 * ny;
  } else {
    nbfaces_test = nx * ny * 2 + nx * nz * 2 + ny * nz * 2;
  }
  CHECK_EQUAL(nbfaces_test, nbfaces);

  const auto& bfaces = *mesh.getMap(AmanziMesh::Entity_kind::BOUNDARY_FACE, true);
  const auto& faces = *mesh.getMap(AmanziMesh::Entity_kind::FACE, true);

  Kokkos::parallel_for(
    "for: testExteriorMapsUnitBox", bfaces.getLocalNumElements(), KOKKOS_LAMBDA(const int j) {
      auto bf = faces.getLocalElement(bfaces.getGlobalElement(j));
      auto f_centroid = AmanziMesh::Impl::get(m).getFaceCentroid(bf);
      bool found = false;
      for (int i = 0; i != AmanziMesh::Impl::get(m).getManifoldDimension(); ++i) {
        if (std::abs(f_centroid[i]) < 1e-10 || std::abs(f_centroid[i] - 1) < 1e-10) {
          found = true;
          break;
        }
      }

      if (!found) {
        printf("not found: %d at %g,%g,%g \n", bf, f_centroid[0], f_centroid[1], f_centroid[2]);
        assert(found);
      }
      //CHECK(found);
    });

  // check nodes are on the boundary
  //
  // NOTE: this appears broken in current master, see #583
  int nbnodes = mesh.getMap(AmanziMesh::Entity_kind::BOUNDARY_NODE, false)->getGlobalNumElements();
  int nbnodes_test;
  if (nz < 0) {
    nbnodes_test = 2 * (nx - 1) + 2 * (ny - 1) + 4; // don't double count the corners
  } else {
    nbnodes_test = 2 * (nx - 1) * (ny - 1) + 2 * (nx - 1) * (nz - 1) + 2 * (ny - 1) * (nz - 1) +
                   4 * (nx - 1) + 4 * (ny - 1) + 4 * (nz - 1) + 8;
  }
  CHECK_EQUAL(nbnodes_test, nbnodes);

  const auto& bnodes = *mesh.getMap(AmanziMesh::Entity_kind::BOUNDARY_NODE, true);
  const auto& nodes = *mesh.getMap(AmanziMesh::Entity_kind::NODE, true);
  Kokkos::parallel_for(
    "for: testExteriorMapUnitBox", bnodes.getLocalNumElements(), KOKKOS_LAMBDA(const int j) {
      printf(" bnode %d GID %d LID %d\n",
             j,
             bnodes.getGlobalElement(j),
             nodes.getLocalElement(bnodes.getGlobalElement(j)));

      auto bn = nodes.getLocalElement(bnodes.getGlobalElement(j));
      AmanziGeometry::Point nc;
      nc = AmanziMesh::Impl::get(m).getNodeCoordinate(bn);
      bool found = false;
      if (std::abs(nc[0]) < 1e-10 || std::abs(nc[0] - 1) < 1e-10 || std::abs(nc[1]) < 1e-10 ||
          std::abs(nc[1] - 1) < 1e-10 || std::abs(nc[2]) < 1e-10 || std::abs(nc[2] - 1) < 1e-10) {
        found = true;
      }
      assert(found);
      //CHECK(found);
    });
}


//
// Test a columnar system
//
inline void
testColumnsUniformDz(const Mesh& mesh, double dz)
{
  // tests the columnar structure of cells
  int n_columns = mesh.columns->num_columns_all;
  CHECK(n_columns > 0);

  // also tests that cols with ghost entities are listed first
  int ncells_owned =
    mesh.getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);
  bool owned = true;

  for (int col = 0; col != n_columns; ++col) {
    const auto& cells = mesh.columns->cells_.template getRowUnmanaged<MemSpace_kind::HOST>(col);

    // check all owned cells first, then all ghosted
    if (owned) {
      if (cells[0] < ncells_owned) {
        for (const auto& c : cells) { CHECK(c < ncells_owned); }
      } else {
        owned = false;
      }
    }
    if (!owned) {
      for (const auto& c : cells) { CHECK(c >= ncells_owned); }
    }

    // check geometry
    const auto& faces = mesh.columns->faces_.template getRowUnmanaged<MemSpace_kind::HOST>(col);
    CHECK(faces.size() == (cells.size() + 1));
    for (int i = 0; i != cells.size(); ++i) {
      Entity_ID c = cells[i];
      AmanziGeometry::Point cc = mesh.getCellCentroid(c);
      std::cout << "col " << col << " cell(" << i << ") is " << c << " with cent = " << cc
                << std::endl;
      AmanziGeometry::Point fd = mesh.getFaceCentroid(faces[i + 1]);
      AmanziGeometry::Point fu = mesh.getFaceCentroid(faces[i]);
      CHECK_CLOSE(cc[0], fu[0], 1e-10);
      CHECK_CLOSE(cc[1], fu[1], 1e-10);
      CHECK_CLOSE(cc[0], fd[0], 1e-10);
      CHECK_CLOSE(cc[1], fd[1], 1e-10);
      CHECK_CLOSE(cc[2], fd[2] + dz / 2.0, 1e-10);
      CHECK_CLOSE(cc[2], fu[2] - dz / 2.0, 1e-10);

      if (i + 1 != cells.size()) {
        AmanziGeometry::Point cc_d = mesh.getCellCentroid(cells[i + 1]);
        CHECK_CLOSE(cc[0], cc_d[0], 1.e-10);
        CHECK_CLOSE(cc[1], cc_d[1], 1.e-10);
        CHECK_CLOSE(cc[2], cc_d[2] + dz, 1.e-10);
      }
    }
  }

  // // test the columnar structure of nodes
  // int nnodes = mesh.getNumEntities(AmanziMesh::Entity_kind::NODE, AmanziMesh::Parallel_kind::ALL);
  // for (int n=0; n!=nnodes; ++n) {
  //   AmanziGeometry::Point nc;
  //   nc = mesh.getNodeCoordinate(n);

  //   int nu = mesh.node_get_node_above(n);
  //   if (nu >= 0) {
  //     AmanziGeometry::Point nuc;
  //     nuc = mesh.getNodeCoordinate(nu);
  //     CHECK_CLOSE(nc[0], nuc[0], 1.e-10);
  //     CHECK_CLOSE(nc[1], nuc[1], 1.e-10);
  //     CHECK_CLOSE(nc[2], nuc[2]-dz, 1.e-10);
  //   }
  // }
}
