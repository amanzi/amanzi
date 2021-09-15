/*
  Copyright 2010-201x held jointly by LANL, ORNL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

#pragma once

#include "Mesh.hh"
#include "MeshFactory.hh"
#include "MeshAudit.hh"


using namespace Amanzi;
using namespace Amanzi::AmanziMesh;

using Entity_Direction_List = std::vector<int>;

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
                 const std::vector<AmanziGeometry::Point>& exp_cell_centroids,
                 const std::vector<double>& exp_cell_volumes,
                 const std::vector<AmanziGeometry::Point>& exp_face_centroids,
                 const std::vector<double>& exp_face_areas,
                 const std::vector<AmanziGeometry::Point>& exp_face_normals,
                 const std::vector<AmanziGeometry::Point>& exp_node_coordinates)
{
  // test cell-based quantities
  {
    std::vector<int> found(exp_cell_centroids.size(), false);
    int ncells = mesh->num_entities(Entity_kind::CELL,Parallel_type::OWNED);
    for (int i = 0; i < ncells; i++) {
      auto centroid = mesh->cell_centroid(i);

      // Search for a cell with the same centroid in the
      // expected list of centroid, check volume
      int j = 0;
      for (; j!=exp_cell_centroids.size(); ++j) {
        auto diff = exp_cell_centroids[j] - centroid;
        if (AmanziGeometry::norm(diff) < 1.0e-10) {
          CHECK_CLOSE(exp_cell_volumes[j], mesh->cell_volume(i), 1e-10);
          break;
        }
      }
      bool lfound = (j < exp_cell_volumes.size());
      CHECK(lfound);
      if (lfound) found[j] = 1;

      // check that the outward normals sum to 0
      Entity_ID_List cfaces;
      mesh->cell_get_faces(i, &cfaces);
      AmanziGeometry::Point normal_sum(mesh->space_dimension());
      for (int j = 0; j < cfaces.size(); j++) {
        auto normal = mesh->face_normal(cfaces[j],false,i);
        normal_sum = normal_sum + normal;
      }
      double val = AmanziGeometry::norm(normal_sum);
      if (val > 1.e-10) {
        for (int j = 0; j < cfaces.size(); j++) {
          auto normal = mesh->face_normal(cfaces[j],false,i);
          std::cout << "fail cell " << i << " normal (" << j << ") " << cfaces[j] << " = " << normal << std::endl;
        }
      }
      CHECK_CLOSE(0., val, 1.0e-10);
    }
    CHECK_MPI_ALL(found, *mesh->get_comm());
  }

  // test face-based quantities
  {
    std::vector<int> found(exp_face_normals.size(), false);

    int nfaces = mesh->num_entities(Entity_kind::FACE,Parallel_type::OWNED);
    for (int i = 0; i < nfaces; i++) {
      AmanziGeometry::Point centroid = mesh->face_centroid(i);

      int j = 0;
      for (; j < found.size(); ++j) {
        auto diff = exp_face_centroids[j] - centroid;
        if (AmanziGeometry::norm(diff) < 1.0e-10) {
          CHECK_CLOSE(exp_face_areas[j], mesh->face_area(i), 1.e-10);

          // Natural normal is well-posed
          AmanziGeometry::Point natural_normal = mesh->face_normal(i);

          // Check the normal with respect to each connected cell is given as the
          // natural times the orientation.
          Entity_ID_List cellids;
          mesh->face_get_cells(i,Parallel_type::ALL,&cellids);

          for (int k = 0; k < cellids.size(); k++) {
            int orientation = 0;
            auto normal_wrt_cell = mesh->face_normal(i, false, cellids[k], &orientation);
            CHECK(natural_normal * orientation == normal_wrt_cell);

            // check the cell's outward normal is indeed outward (assumes star-convex)
            AmanziGeometry::Point cellcentroid = mesh->cell_centroid(cellids[k]);
            AmanziGeometry::Point facecentroid = mesh->face_centroid(i);
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
    CHECK_MPI_ALL(found, *mesh->get_comm());
  }

  // test the node-based quantities
  {
    std::vector<int> found(exp_node_coordinates.size(), false);
    int nnodes = mesh->num_entities(Entity_kind::NODE,Parallel_type::OWNED);
    for (int i = 0; i < nnodes; i++) {
      AmanziGeometry::Point centroid;
      mesh->node_get_coordinates(i, &centroid);
      int j = 0;
      for (; j < found.size(); ++j) {
        auto diff = exp_node_coordinates[j] - centroid;
        if (AmanziGeometry::norm(diff) < 1.0e-10) break;
      }

      bool lfound = (j < exp_node_coordinates.size());
      CHECK(lfound);
      if (lfound) found[j] = 1;
    }
    CHECK_MPI_ALL(found, *mesh->get_comm());
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
  CHECK_EQUAL(2, mesh->space_dimension());
  int ncells = mesh->num_entities(Entity_kind::CELL,Parallel_type::OWNED);
  int ncells_test = nx * ny;
  CHECK_CLOSE_SUMALL(ncells_test, ncells, *mesh->get_comm());

  int nfaces = mesh->num_entities(Entity_kind::FACE,Parallel_type::OWNED);
  int nfaces_test = ny * (nx+1) + nx * (ny+1);
  CHECK_CLOSE_SUMALL(nfaces_test, nfaces, *mesh->get_comm());

  int nnodes = mesh->num_entities(Entity_kind::NODE,Parallel_type::OWNED);
  int nnodes_test = (ny+1) * (nx+1);
  CHECK_CLOSE_SUMALL(nnodes_test, nnodes, *mesh->get_comm());

  // construct expected cell volumes, centroids
  std::vector<double> exp_cell_volumes(ncells_test, 1./nx * 1./ny);
  std::vector<AmanziGeometry::Point> exp_cell_centroids;
  for (int i=0; i!=nx; ++i) {
    for (int j=0; j!=ny; ++j) {
      exp_cell_centroids.emplace_back(AmanziGeometry::Point((i+.5)/nx, (j+.5)/ny));
    }
  }

  // construct expected face centroids, areas, and normals
  std::vector<AmanziGeometry::Point> exp_face_centroids;
  std::vector<double> exp_face_areas;
  std::vector<AmanziGeometry::Point> exp_face_normals;
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
  std::vector<AmanziGeometry::Point> exp_node_coordinates;
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
  CHECK_EQUAL(3, mesh->space_dimension());
  int ncells = mesh->num_entities(Entity_kind::CELL,Parallel_type::OWNED);
  int ncells_test = nx * ny * nz;
  CHECK_CLOSE_SUMALL(ncells_test, ncells, *mesh->get_comm());
  int nfaces = mesh->num_entities(Entity_kind::FACE,Parallel_type::OWNED);
  int nfaces_test = nx * ny * (nz+1) + nx * (ny+1) * nz + (nx+1) * ny * nz;
  CHECK_CLOSE_SUMALL(nfaces_test, nfaces, *mesh->get_comm());
  int nnodes = mesh->num_entities(Entity_kind::NODE,Parallel_type::OWNED);
  int nnodes_test = (nx+1) * (ny+1) * (nz+1);
  CHECK_CLOSE_SUMALL(nnodes_test, nnodes, *mesh->get_comm());

  // construct expected cell volumes, centroids
  std::vector<double> exp_cell_volumes(ncells_test, 1./nx * 1./ny * 1./nz);
  std::vector<AmanziGeometry::Point> exp_cell_centroids;
  for (int i=0; i!=nx; ++i) {
    for (int j=0; j!=ny; ++j) {
      for (int k=0; k!=nz; ++k) {
        exp_cell_centroids.emplace_back(AmanziGeometry::Point((i+.5)/nx, (j+.5)/ny, (k+.5)/nz));
      }
    }
  }

  // construct expected face centroids, areas, and normals
  std::vector<AmanziGeometry::Point> exp_face_centroids;
  std::vector<double> exp_face_areas;
  std::vector<AmanziGeometry::Point> exp_face_normals;
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
  std::vector<AmanziGeometry::Point> exp_node_coordinates;
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
  int nbfaces = mesh->exterior_face_map(false).NumGlobalElements();
  int nbfaces_test;
  if (nz < 0) {
    nbfaces_test = 2*nx + 2*ny;
  } else {
    nbfaces_test = nx*ny*2 + nx*nz*2 + ny*nz*2;
  }
  CHECK_EQUAL(nbfaces_test, nbfaces);

  auto& bfaces = mesh->exterior_face_map(true);
  auto& faces = mesh->map(AmanziMesh::Entity_kind::FACE, true);
  for (int j=0; j!=bfaces.NumMyElements(); ++j) {
    auto bf = faces.LID(bfaces.GID(j));
    auto f_centroid = mesh->face_centroid(bf);
    bool found = false;
    for (int i=0; i!=mesh->manifold_dimension(); ++i) {
      if (std::abs(f_centroid[i]) < 1e-10 ||
          std::abs(f_centroid[i] - 1) < 1e-10) {
        found = true;
      }
    }
    CHECK(found);
  }

  // check nodes are on the boundary
  //
  // NOTE: this appears broken in current master, see #583
  int nbnodes = mesh->exterior_node_map(false).NumGlobalElements();
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

  auto& bnodes = mesh->exterior_node_map(true);
  auto& nodes = mesh->map(AmanziMesh::Entity_kind::NODE, true);
  for (int j=0; j!=bnodes.NumMyElements(); ++j) {
    std::cout << " bnode " << j << " GID " << bnodes.GID(j) << " LID " << nodes.LID(bnodes.GID(j)) << std::endl;

    auto bn = nodes.LID(bnodes.GID(j));
    AmanziGeometry::Point nc;
    mesh->node_get_coordinates(bn, &nc);
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
template<class Mesh_type>
void
testColumnsUniformDz(const Teuchos::RCP<Mesh_type>& mesh, double dz)
{
  // tests the columnar structure of cells
  int n_columns = mesh->num_columns(true);

  // also tests that cols with ghost entities are listed first
  int ncells_owned = mesh->num_entities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_type::OWNED);
  bool owned = true;

  for (int col=0; col!=n_columns; ++col) {
    const auto& cells = mesh->cells_of_column(col);

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
    const auto& faces = mesh->faces_of_column(col);
    CHECK(faces.size() == (cells.size() + 1));
    for (int i=0; i!=cells.size(); ++i) {
      Entity_ID c = cells[i];
      AmanziGeometry::Point cc = mesh->cell_centroid(c);
      AmanziGeometry::Point fd = mesh->face_centroid(faces[i+1]);
      AmanziGeometry::Point fu = mesh->face_centroid(faces[i]);
      CHECK_CLOSE(cc[0], fu[0], 1e-10);
      CHECK_CLOSE(cc[1], fu[1], 1e-10);
      CHECK_CLOSE(cc[0], fd[0], 1e-10);
      CHECK_CLOSE(cc[1], fd[1], 1e-10);
      CHECK_CLOSE(cc[2], fd[2] + dz/2.0, 1e-10);
      CHECK_CLOSE(cc[2], fu[2] - dz/2.0, 1e-10);

      if (i != 0) CHECK_EQUAL(cells[i-1], mesh->cell_get_cell_above(c));
      if (i != cells.size()-1) CHECK_EQUAL(cells[i+1], mesh->cell_get_cell_below(c));
    }
  }

  // test the columnar structure of nodes
  int nnodes = mesh->num_entities(AmanziMesh::Entity_kind::NODE, AmanziMesh::Parallel_type::ALL);
  for (int n=0; n!=nnodes; ++n) {
    AmanziGeometry::Point nc;
    mesh->node_get_coordinates(n, &nc);

    int nu = mesh->node_get_node_above(n);
    if (nu >= 0) {
      AmanziGeometry::Point nuc;
      mesh->node_get_coordinates(nu, &nuc);
      CHECK_CLOSE(nc[0], nuc[0], 1.e-10);
      CHECK_CLOSE(nc[1], nuc[1], 1.e-10);
      CHECK_CLOSE(nc[2], nuc[2]-dz, 1.e-10);
    }
  }

}
