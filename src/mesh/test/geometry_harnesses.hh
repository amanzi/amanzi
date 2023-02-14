/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

#pragma once

#include "MeshFrameworkTraits.hh"

using namespace Amanzi;
using namespace Amanzi::AmanziMesh;

using Entity_Direction_View = std::vector<int>;

//
// Helper functions for testing geometry
//

//
// This is a helper function -- simply runs MeshAudit
//
template<class MeshAudit_type, class Mesh_type>
void
testMeshAudit(const Teuchos::RCP<Mesh_type>& mesh) {
  // run MeshAudit
  MeshAudit_type audit(mesh);
  int status = audit.Verify();
  CHECK_EQUAL(0, status);
}

//
// Sums a scalar over all processes and compares the result to exp
//
// Used to check that partial counts add up to global count in parallel tests.
//
template<typename T>
void
CHECK_CLOSE_SUMALL(T exp, T contrib, const Amanzi::Comm_type& comm, T tol=0) {
  // MPI-based CHECK_CLOSE using SumAll
  T global;
  comm.SumAll(&contrib, &global, 1);
  CHECK_CLOSE(exp, global, tol);
  if (std::abs(exp - global) > tol) {
    // for debugging
    double diff = std::abs(exp - global);
  }
}

//
// Checks that contrib has a nonzero entry in every i of the array on at least
// one proc.  Used in parallel tests to ensure that the (globally expected
// array) is found and checked on at least one owning rank.
//
void
CHECK_MPI_ALL(std::vector<int>& contrib, const Amanzi::Comm_type& comm) {
  // MPI-based, confirms that all entries of contrib are true on some process.
  std::vector<int> global(contrib.size());
  comm.MaxAll(contrib.data(), global.data(), contrib.size());
  CHECK(std::all_of(global.begin(), global.end(), [](bool cond){ return cond; }));
}

//
// Tests geometry given expected values
//
template<class Mesh_type>
void
testMeshGeometry(const Teuchos::RCP<Mesh_type>& mesh,
                 const Point_List& exp_cell_centroids,
                 const Double_List& exp_cell_volumes,
                 const Point_List& exp_face_centroids,
                 const Double_List& exp_face_areas,
                 const Point_List& exp_face_normals,
                 const Point_List& exp_node_coordinates)
{
  // test cell-based quantities
  {
    std::vector<int> found(exp_cell_centroids.size(), false);
    int ncells = mesh->getNumEntities(Entity_kind::CELL,Parallel_kind::OWNED);
    for (int i = 0; i < ncells; i++) {
      auto centroid = mesh->getCellCentroid(i);

      // Search for a cell with the same centroid in the
      // expected list of centroid, check volume
      int j = 0;
      for (; j!=exp_cell_centroids.size(); ++j) {
        auto diff = exp_cell_centroids[j] - centroid;
        if (AmanziGeometry::norm(diff) < 1.0e-10) {
          CHECK_CLOSE(exp_cell_volumes[j], mesh->getCellVolume(i), 1e-10);
          break;
        }
      }
      bool lfound = (j < exp_cell_volumes.size());
      CHECK(lfound);
      if (lfound) found[j] = 1;

      // check that the outward normals sum to 0
      cEntity_ID_View cfaces;
      mesh->getCellFaces(i, cfaces);
      AmanziGeometry::Point normal_sum(mesh->getSpaceDimension());
      for (int j = 0; j < cfaces.size(); j++) {
        auto normal = mesh->getFaceNormal(cfaces[j], i);
        normal_sum = normal_sum + normal;
      }
      double val = AmanziGeometry::norm(normal_sum);
      if (val > 1.e-10) {
        for (int j = 0; j < cfaces.size(); j++) {
          auto normal = mesh->getFaceNormal(cfaces[j], i);
          std::cout << "fail cell " << i << " normal (" << j << ") " << cfaces[j] << " = " << normal << std::endl;
        }
      }
      CHECK_CLOSE(0., val, 1.0e-10);
    }
    CHECK_MPI_ALL(found, *mesh->getComm());
  }

  // test face-based quantities
  {
    std::vector<int> found(exp_face_normals.size(), false);

    int nfaces = mesh->getNumEntities(Entity_kind::FACE,Parallel_kind::OWNED);
    for (int i = 0; i < nfaces; i++) {
      AmanziGeometry::Point centroid = mesh->getFaceCentroid(i);

      int j = 0;
      for (; j < found.size(); ++j) {
        auto diff = exp_face_centroids[j] - centroid;
        if (AmanziGeometry::norm(diff) < 1.0e-10) {
          CHECK_CLOSE(exp_face_areas[j], mesh->getFaceArea(i), 1.e-10);

          // Natural normal is well-posed
          AmanziGeometry::Point natural_normal = mesh->getFaceNormal(i);

          // Check the normal with respect to each connected cell is given as the
          // natural times the orientation.
          cEntity_ID_View cellids;
          mesh->getFaceCells(i,Parallel_kind::ALL,cellids);

          for (int k = 0; k < cellids.size(); k++) {
            int orientation = 0;
            auto normal_wrt_cell = mesh->getFaceNormal(i, cellids[k], &orientation);
            if (natural_normal * orientation != normal_wrt_cell) {
              std::cout << "Fail face: " << i << " wrt cell " << cellids[k] << std::endl
                        << "  nat_normal = " << natural_normal << std::endl
                        << "  normal_wrt = " << normal_wrt_cell << std::endl
                        << "  orientation = " << orientation << std::endl;
            }
            CHECK(natural_normal * orientation == normal_wrt_cell);

            // check the cell's outward normal is indeed outward (assumes star-convex)
            AmanziGeometry::Point cellcentroid = mesh->getCellCentroid(cellids[k]);
            AmanziGeometry::Point facecentroid = mesh->getFaceCentroid(i);
            AmanziGeometry::Point outvec = facecentroid-cellcentroid;

            double dp = outvec * normal_wrt_cell;
            dp /= (AmanziGeometry::norm(outvec) * AmanziGeometry::norm(normal_wrt_cell));
            CHECK_CLOSE(1., dp, 1e-10);
          }
          break;
        }
      }
      bool lfound = (j < exp_face_areas.size());
      CHECK(lfound);
      if (lfound) found[j] = 1;
    }
    CHECK_MPI_ALL(found, *mesh->getComm());
  }

  // test the node-based quantities
  {
    std::vector<int> found(exp_node_coordinates.size(), false);
    int nnodes = mesh->getNumEntities(Entity_kind::NODE,Parallel_kind::OWNED);
    for (int i = 0; i < nnodes; i++) {
      AmanziGeometry::Point centroid;
      centroid = mesh->getNodeCoordinate(i);
      int j = 0;
      for (; j < found.size(); ++j) {
        auto diff = exp_node_coordinates[j] - centroid;
        if (AmanziGeometry::norm(diff) < 1.0e-10) break;
      }

      bool lfound = (j < exp_node_coordinates.size());
      CHECK(lfound);
      if (lfound) found[j] = 1;
    }
    CHECK_MPI_ALL(found, *mesh->getComm());
  }
}

//
// Form the expected values and call testGeometry for a 2D box
//
template<class Mesh_type>
void
testGeometryQuad(const Teuchos::RCP<Mesh_type>& mesh, int nx, int ny)
{
  // test the basic dimensionality
  CHECK_EQUAL(2, mesh->getSpaceDimension());
  int ncells = mesh->getNumEntities(Entity_kind::CELL,Parallel_kind::OWNED);
  int ncells_test = nx * ny;
  CHECK_CLOSE_SUMALL(ncells_test, ncells, *mesh->getComm());

  int nfaces = mesh->getNumEntities(Entity_kind::FACE,Parallel_kind::OWNED);
  int nfaces_test = ny * (nx+1) + nx * (ny+1);
  CHECK_CLOSE_SUMALL(nfaces_test, nfaces, *mesh->getComm());

  int nnodes = mesh->getNumEntities(Entity_kind::NODE,Parallel_kind::OWNED);
  int nnodes_test = (ny+1) * (nx+1);
  CHECK_CLOSE_SUMALL(nnodes_test, nnodes, *mesh->getComm());

  // construct expected cell volumes, centroids
  Double_List exp_cell_volumes(ncells_test, 1./nx * 1./ny);
  Point_List exp_cell_centroids;
  for (int i=0; i!=nx; ++i) {
    for (int j=0; j!=ny; ++j) {
      exp_cell_centroids.emplace_back(AmanziGeometry::Point((i+.5)/nx, (j+.5)/ny));
    }
  }

  // construct expected face centroids, areas, and normals
  Point_List exp_face_centroids;
  Double_List exp_face_areas;
  Point_List exp_face_normals;
  for (int i=0; i!=nx+1; ++i) {
    for (int j=0; j!=ny; ++j) {
      exp_face_centroids.emplace_back(AmanziGeometry::Point( ((double) i) / nx, (j+.5)/ny));
      exp_face_areas.emplace_back(1.0/ny);
      exp_face_normals.emplace_back(AmanziGeometry::Point(1.0, 0.0));
    }
  }
  for (int i=0; i!=nx; ++i) {
    for (int j=0; j!=ny+1; ++j) {
      exp_face_centroids.emplace_back(AmanziGeometry::Point((i+.5)/nx, ((double) j) / ny));
      exp_face_areas.emplace_back(1.0/nx);
      exp_face_normals.emplace_back(AmanziGeometry::Point(0.0, 1.0));
    }
  }

  // construct expected nodal locations
  Point_List exp_node_coordinates;
  for (int i=0; i!=nx+1; ++i) {
    for (int j=0; j!=ny+1; ++j) {
      exp_node_coordinates.emplace_back(AmanziGeometry::Point(((double) i)/nx, ((double) j)/ny));
    }
  }

  // run the testing
  testMeshGeometry(mesh, exp_cell_centroids, exp_cell_volumes,
                   exp_face_centroids, exp_face_areas, exp_face_normals,
                   exp_node_coordinates);
}


//
// Form the expected values and call testGeometry for a 3D cube
//
template<class Mesh_type>
void
testGeometryCube(const Teuchos::RCP<Mesh_type>& mesh, int nx, int ny, int nz)
{
  // test the basic dimensionality
  CHECK_EQUAL(3, mesh->getSpaceDimension());
  int ncells = mesh->getNumEntities(Entity_kind::CELL,Parallel_kind::OWNED);
  int ncells_test = nx * ny * nz;
  CHECK_CLOSE_SUMALL(ncells_test, ncells, *mesh->getComm());
  int nfaces = mesh->getNumEntities(Entity_kind::FACE,Parallel_kind::OWNED);
  int nfaces_test = nx * ny * (nz+1) + nx * (ny+1) * nz + (nx+1) * ny * nz;
  CHECK_CLOSE_SUMALL(nfaces_test, nfaces, *mesh->getComm());
  int nnodes = mesh->getNumEntities(Entity_kind::NODE,Parallel_kind::OWNED);
  int nnodes_test = (nx+1) * (ny+1) * (nz+1);
  CHECK_CLOSE_SUMALL(nnodes_test, nnodes, *mesh->getComm());

  // construct expected cell volumes, centroids
  Double_List exp_cell_volumes(ncells_test, 1./nx * 1./ny * 1./nz);
  Point_List exp_cell_centroids;
  for (int i=0; i!=nx; ++i) {
    for (int j=0; j!=ny; ++j) {
      for (int k=0; k!=nz; ++k) {
        exp_cell_centroids.emplace_back(AmanziGeometry::Point((i+.5)/nx, (j+.5)/ny, (k+.5)/nz));
      }
    }
  }

  // construct expected face centroids, areas, and normals
  Point_List exp_face_centroids;
  Double_List exp_face_areas;
  Point_List exp_face_normals;
  for (int i=0; i!=(nx+1); ++i) {
    for (int j=0; j!=ny; ++j) {
      for (int k=0; k!=nz; ++k) {
        exp_face_centroids.emplace_back(AmanziGeometry::Point( ((double) i) / nx, (j+.5)/ny, (k+.5)/nz));
        exp_face_areas.emplace_back(1.0/ny * 1.0/nz);
        exp_face_normals.emplace_back(AmanziGeometry::Point(1.0, 0.0, 0.0));
      }
    }
  }
  for (int i=0; i!=nx; ++i) {
    for (int j=0; j!=(ny+1); ++j) {
      for (int k=0; k!=nz; ++k) {
        exp_face_centroids.emplace_back(AmanziGeometry::Point((i+.5)/nx, ((double) j) / ny, (k+.5)/nz));
        exp_face_areas.emplace_back(1.0/nx * 1.0/nz);
        exp_face_normals.emplace_back(AmanziGeometry::Point(0.0, 1.0, 0.0));
      }
    }
  }
  for (int i=0; i!=nx; ++i) {
    for (int j=0; j!=ny; ++j) {
      for (int k=0; k!=(nz+1); ++k) {
        exp_face_centroids.emplace_back(AmanziGeometry::Point((i+.5)/nx, (j+.5)/ny, ((double) k) / nz));
        exp_face_areas.emplace_back(1.0/nx * 1.0/ny);
        exp_face_normals.emplace_back(AmanziGeometry::Point(0.0, 0.0, 1.0));
      }
    }
  }

  // construct expected nodal locations
  Point_List exp_node_coordinates;
  for (int i=0; i!=nx+1; ++i) {
    for (int j=0; j!=ny+1; ++j) {
      for (int k=0; k!=nz+1; ++k) {
        exp_node_coordinates.emplace_back(AmanziGeometry::Point(((double) i)/nx, ((double) j)/ny, ((double) k)/nz));
      }
    }
  }

  // run the testing
  testMeshGeometry(mesh, exp_cell_centroids, exp_cell_volumes,
                   exp_face_centroids, exp_face_areas, exp_face_normals,
                   exp_node_coordinates);
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
template<class Mesh_type>
void
testExteriorMapsUnitBox(const Teuchos::RCP<Mesh_type>& mesh, int nx, int ny, int nz=-1)
{
  // check faces are on the boundary
  int nbfaces = mesh->getMap(AmanziMesh::Entity_kind::BOUNDARY_FACE, false).NumGlobalElements();
  int nbfaces_test;
  if (nz < 0) {
    nbfaces_test = 2*nx + 2*ny;
  } else {
    nbfaces_test = nx*ny*2 + nx*nz*2 + ny*nz*2;
  }
  CHECK_EQUAL(nbfaces_test, nbfaces);

  auto& bfaces = mesh->getMap(AmanziMesh::Entity_kind::BOUNDARY_FACE, true);
  auto& faces = mesh->getMap(AmanziMesh::Entity_kind::FACE, true);
  auto bface_ids = mesh->getBoundaryFaces();

  for (int j=0; j!=bfaces.NumMyElements(); ++j) {
    auto bf = faces.LID(bfaces.GID(j));
    CHECK_EQUAL(bface_ids[j], bf);
    auto f_centroid = mesh->getFaceCentroid(bf);
    bool found = false;
    for (int i=0; i!=mesh->getManifoldDimension(); ++i) {
      if (std::abs(f_centroid[i]) < 1e-10 ||
          std::abs(f_centroid[i] - 1) < 1e-10) {
        found = true;
        break;
      }
    }
    if (!found) {
      std::cout << "not found: " << bf << " at " << f_centroid << std::endl;
    }
    CHECK(found);
  }

  // check nodes are on the boundary
  //
  // NOTE: this appears broken in current master, see #583
  int nbnodes = mesh->getMap(AmanziMesh::Entity_kind::BOUNDARY_NODE, false).NumGlobalElements();
  int nbnodes_test;
  if (nz < 0) {
    nbnodes_test = 2*(nx-1) + 2*(ny-1) + 4; // don't double count the corners
  } else {
    nbnodes_test = 2*(nx-1)*(ny-1) +
      2*(nx-1)*(nz-1) +
      2*(ny-1)*(nz-1) +
      4*(nx-1) + 4*(ny-1) + 4*(nz-1)
      + 8;
  }
  CHECK_EQUAL(nbnodes_test, nbnodes);

  auto& bnodes = mesh->getMap(AmanziMesh::Entity_kind::BOUNDARY_NODE, true);
  auto& nodes = mesh->getMap(AmanziMesh::Entity_kind::NODE, true);
  for (int j=0; j!=bnodes.NumMyElements(); ++j) {
    std::cout << " bnode " << j << " GID " << bnodes.GID(j) << " LID " << nodes.LID(bnodes.GID(j)) << std::endl;

    auto bn = nodes.LID(bnodes.GID(j));
    AmanziGeometry::Point nc;
    nc = mesh->getNodeCoordinate(bn);
    bool found = false;
    if (std::abs(nc[0]) < 1e-10 ||
        std::abs(nc[0] - 1) < 1e-10 ||
        std::abs(nc[1]) < 1e-10 ||
        std::abs(nc[1] - 1) < 1e-10 ||
        std::abs(nc[2]) < 1e-10 ||
        std::abs(nc[2] - 1) < 1e-10) {
      found = true;
    }
    CHECK(found);
  }
}


//
// Test a columnar system
//
template<MemSpace_kind MEM>
inline
void
testColumnsUniformDz(const MeshCache<MEM>& mesh, double dz)
{
  // tests the columnar structure of cells
  int n_columns = mesh.columns.num_columns_all;

  // also tests that cols with ghost entities are listed first
  int ncells_owned = mesh.getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);
  bool owned = true;

  for (int col=0; col!=n_columns; ++col) {
    const auto& cells = mesh.columns.cells_.template getRow<MEM>(col);

    // check all owned cells first, then all ghosted
    if (owned) {
      if (cells[0] < ncells_owned) {
        for (const auto& c : cells) {
          CHECK(c < ncells_owned);
        }
      } else {
        owned = false;
      }
    }
    if (!owned) {
      for (const auto& c : cells) {
        CHECK(c >= ncells_owned);
      }
    }

    // check geometry
    const auto& faces = mesh.columns.faces_.template getRow<MEM>(col);
    CHECK(faces.size() == (cells.size() + 1));
    for (int i=0; i!=cells.size(); ++i) {
      Entity_ID c = cells[i];
      AmanziGeometry::Point cc = mesh.getCellCentroid(c);
      std::cout << "col " << col << " cell(" << i << ") is " << c << " with cent = " << cc << std::endl;
      AmanziGeometry::Point fd = mesh.getFaceCentroid(faces[i+1]);
      AmanziGeometry::Point fu = mesh.getFaceCentroid(faces[i]);
      CHECK_CLOSE(cc[0], fu[0], 1e-10);
      CHECK_CLOSE(cc[1], fu[1], 1e-10);
      CHECK_CLOSE(cc[0], fd[0], 1e-10);
      CHECK_CLOSE(cc[1], fd[1], 1e-10);
      CHECK_CLOSE(cc[2], fd[2] + dz/2.0, 1e-10);
      CHECK_CLOSE(cc[2], fu[2] - dz/2.0, 1e-10);

      if (i+1!=cells.size()) {
        AmanziGeometry::Point cc_d = mesh.getCellCentroid(cells[i+1]);
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
