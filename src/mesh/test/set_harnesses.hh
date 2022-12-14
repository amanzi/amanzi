/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

#pragma once

#include "Mesh.hh"
#include "MeshFactory.hh"
#include "MeshAudit.hh"

#include "geometry_harnesses.hh"

using namespace Amanzi;
using namespace Amanzi::AmanziMesh;

using Entity_Direction_List = std::vector<int>;


//
// Helper functions for testing sets
//

//
// Tests Mesh sets on a 3x3x3 hex
//
template <class Mesh_type>
void
testHexMeshSets3x3x3(const Teuchos::RCP<Mesh_type>& mesh, bool labeled, Framework f)
{
  auto gm = mesh->geometric_model();
  CHECK(32 == gm->size() || 33 == gm->size()); // some tests add an extra set
  CHECK_EQUAL(3, mesh->space_dimension());
  CHECK_EQUAL(3, mesh->manifold_dimension());
  for (const auto& r : *gm) {
    std::string r_name = r->get_name();

    if (r_name == "Entire Mesh") {
      if (!mesh->valid_set_type(r->get_type(), AmanziMesh::Entity_kind::CELL)) continue;

      // check that the sum of n_ents is the global number of cells, 3*3*3
      int n_ents =
        mesh->get_set_size(r_name, AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_type::OWNED);
      CHECK_CLOSE_SUMALL(3 * 3 * 3, n_ents, *mesh->get_comm());

      // check that the number of entities is the number of owned cells
      CHECK_EQUAL(
        mesh->num_entities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_type::OWNED),
        n_ents);
      // check that the number of entities is the number of all cells
      CHECK_EQUAL(
        mesh->num_entities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_type::ALL),
        mesh->get_set_size(r_name, AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_type::ALL));

    } else if (r_name == "Bottom LS" || r_name == "Cell Set 1" || r_name == "Bottom ColFunc") {
      if (!mesh->valid_set_type(r->get_type(), AmanziMesh::Entity_kind::CELL)) continue;
      if (r_name == "Cell Set 1" && f == Framework::MOAB)
        continue; // MOAB cannot handle cell sets, only mat IDs
      // the bottom layer of cells
      if (labeled) {
        int n_ents = mesh->get_set_size(
          r_name, AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_type::OWNED);
        CHECK_CLOSE_SUMALL(3 * 3, n_ents, *mesh->get_comm());

        AmanziMesh::Entity_ID_List ents;
        mesh->get_set_entities(
          r_name, AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_type::ALL, &ents);
        for (const auto& e : ents) {
          auto cc = mesh->cell_centroid(e);
          CHECK_CLOSE(0.5 / 3, cc[2], 1.e-10);
        }
      }

    } else if (r_name == "Middle LS" || r_name == "Cell Set 2" || r_name == "Middle ColFunc" ||
               r_name == "Bottom+Middle Box - Bottom LS") {
      if (!mesh->valid_set_type(r->get_type(), AmanziMesh::Entity_kind::CELL)) continue;
      if (r_name == "Cell Set 2" && f == Framework::MOAB)
        continue; // MOAB cannot handle cell sets, only mat IDs

      // the middle layer of cells
      if (labeled) {
        int n_ents = mesh->get_set_size(
          r_name, AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_type::OWNED);
        CHECK_CLOSE_SUMALL(3 * 3, n_ents, *mesh->get_comm());

        AmanziMesh::Entity_ID_List ents;
        mesh->get_set_entities(
          r_name, AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_type::ALL, &ents);
        for (const auto& e : ents) {
          auto cc = mesh->cell_centroid(e);
          CHECK_CLOSE(0.5, cc[2], 1.e-10);
        }
      }

    } else if (r_name == "Top LS" || r_name == "Cell Set 3" || r_name == "Top ColFunc" ||
               r_name == "Top Box" || r_name == "NOT_Bottom+Middle Box") {
      if (!mesh->valid_set_type(r->get_type(), AmanziMesh::Entity_kind::CELL)) continue;
      if (r_name == "Cell Set 3" && f == Framework::MOAB)
        continue; // MOAB cannot handle cell sets, only mat IDs

      // the top layer of cells
      if (labeled || r_name == "Top Box") {
        int n_ents = mesh->get_set_size(
          r_name, AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_type::OWNED);
        CHECK_CLOSE_SUMALL(3 * 3, n_ents, *mesh->get_comm());

        AmanziMesh::Entity_ID_List ents;
        mesh->get_set_entities(
          r_name, AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_type::ALL, &ents);
        for (const auto& e : ents) {
          auto cc = mesh->cell_centroid(e);
          CHECK_CLOSE(2.5 / 3, cc[2], 1.e-10);
        }
      }

    } else if (r_name == "Left Box") {
      // all should support this type
      CHECK(mesh->valid_set_type(r->get_type(), AmanziMesh::Entity_kind::CELL));
      if (!mesh->valid_set_type(r->get_type(), AmanziMesh::Entity_kind::CELL)) continue;

      // left-most layer (or west-most?)
      int n_ents =
        mesh->get_set_size(r_name, AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_type::OWNED);
      CHECK_CLOSE_SUMALL(3 * 3, n_ents, *mesh->get_comm());

      AmanziMesh::Entity_ID_List ents;
      mesh->get_set_entities(
        r_name, AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_type::ALL, &ents);
      for (const auto& e : ents) {
        auto cc = mesh->cell_centroid(e);
        CHECK_CLOSE(0.5 / 3, cc[0], 1.e-10); // left
      }

    } else if (r_name == "Bottom+Middle Box") {
      if (!mesh->valid_set_type(r->get_type(), AmanziMesh::Entity_kind::CELL)) continue;

      // the middle and bottom layers of cells
      int n_ents =
        mesh->get_set_size(r_name, AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_type::OWNED);
      CHECK_CLOSE_SUMALL(3 * 3 * 2, n_ents, *mesh->get_comm());

      AmanziMesh::Entity_ID_List ents;
      mesh->get_set_entities(
        r_name, AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_type::ALL, &ents);
      for (const auto& e : ents) {
        auto cc = mesh->cell_centroid(e);
        if (cc[2] > 0.4) {
          CHECK_CLOSE(0.5, cc[2], 1.e-10);
        } else {
          CHECK_CLOSE(0.5 / 3, cc[2], 1.e-10);
        }
      }

    } else if (r_name == "Point" || r_name == "Sample Point InCell") {
      if (!mesh->valid_set_type(r->get_type(), AmanziMesh::Entity_kind::CELL)) continue;
      // the central-most cell

      int n_ents =
        mesh->get_set_size(r_name, AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_type::OWNED);
      CHECK_CLOSE_SUMALL(1, n_ents, *mesh->get_comm());

      AmanziMesh::Entity_ID_List ents;
      mesh->get_set_entities(
        r_name, AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_type::ALL, &ents);
      for (const auto& e : ents) {
        auto cc = mesh->cell_centroid(e);
        CHECK_CLOSE(0.5, cc[0], 1.e-10);
        CHECK_CLOSE(0.5, cc[1], 1.e-10);
        CHECK_CLOSE(0.5, cc[2], 1.e-10);
      }

    } else if (r_name == "Bottom LS+Point") {
      if (!mesh->valid_set_type(r->get_type(), AmanziMesh::Entity_kind::CELL)) continue;

      // logical -- union.  Needs labeled sets
      if (labeled) {
        int n_ents = mesh->get_set_size(
          r_name, AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_type::OWNED);
        CHECK_CLOSE_SUMALL(3 * 3 + 1, n_ents, *mesh->get_comm());

        AmanziMesh::Entity_ID_List ents;
        mesh->get_set_entities(
          r_name, AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_type::ALL, &ents);
        for (const auto& e : ents) {
          auto cc = mesh->cell_centroid(e);
          if (cc[2] > 0.4) {
            CHECK_CLOSE(0.5, cc[2], 1.e-10);
          } else {
            CHECK_CLOSE(0.5 / 3, cc[2], 1.e-10);
          }
        }
      }

    } else if (r_name == "Left Box INT Bottom LS") {
      if (!mesh->valid_set_type(r->get_type(), AmanziMesh::Entity_kind::CELL)) continue;

      // 3 cells in a row
      if (labeled) {
        int n_ents = mesh->get_set_size(
          r_name, AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_type::OWNED);
        CHECK_CLOSE_SUMALL(3, n_ents, *mesh->get_comm());

        AmanziMesh::Entity_ID_List ents;
        mesh->get_set_entities(
          r_name, AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_type::ALL, &ents);
        for (const auto& e : ents) {
          auto cc = mesh->cell_centroid(e);
          CHECK_CLOSE(0.5 / 3, cc[0], 1.e-10); // left
          CHECK_CLOSE(0.5 / 3, cc[2], 1.e-10); // bottom
        }
      }

    } else if (r_name == "Enumerated") {
      if (!mesh->valid_set_type(r->get_type(), AmanziMesh::Entity_kind::CELL)) continue;

      // if (framework == Framework::MSTK) { // logical only works in MSTK
      // // cells 1,2,3 -- only supported by MSTK
      // NOT WORKING!
      // int n_ents = mesh->get_set_size(r_name, AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_type::OWNED);
      // CHECK_CLOSE_SUMALL(3, n_ents, *mesh->get_comm());

      // AmanziMesh::Entity_ID_List ents;
      // mesh->get_set_entities(r_name, AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_type::ALL, &ents);
      // std::vector<int> ents_found(3, 0);
      // for (const auto& e : ents) {
      //   int gid = mesh->cell_map(false).GID(e);
      //   if (gid == 1) ents_found[0] = 1;
      //   if (gid == 2) ents_found[1] = 1;
      //   if (gid == 3) ents_found[2] = 1;
      // }
      // CHECK_MPI_ALL(ents_found, *mesh->get_comm());
      // }

    } else if (r_name == "Face 106" || r_name == "Top Face Plane") {
      if (r_name == "Top Face Plane") {
        // all should support this type
        CHECK(mesh->valid_set_type(r->get_type(), AmanziMesh::Entity_kind::FACE));
      }
      if (!mesh->valid_set_type(r->get_type(), AmanziMesh::Entity_kind::FACE)) continue;

      // the top faces
      if (labeled || r_name == "Top Face Plane") {
        int n_ents = mesh->get_set_size(
          r_name, AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_type::OWNED);
        CHECK_CLOSE_SUMALL(3 * 3, n_ents, *mesh->get_comm());

        AmanziMesh::Entity_ID_List ents;
        mesh->get_set_entities(
          r_name, AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_type::ALL, &ents);
        for (const auto& e : ents) {
          auto fc = mesh->face_centroid(e);
          std::cout << r_name << " face: " << fc << std::endl;
          CHECK_CLOSE(1, fc[2], 1.e-10);
        }
      }

    } else if (r_name == "Face 103" || r_name == "West Face Plane") {
      if (r_name == "West Face Plane") {
        // all should support this type
        CHECK(mesh->valid_set_type(r->get_type(), AmanziMesh::Entity_kind::FACE));
      }
      if (!mesh->valid_set_type(r->get_type(), AmanziMesh::Entity_kind::FACE)) continue;

      // the top faces
      if (labeled || r_name == "West Face Plane") {
        int n_ents = mesh->get_set_size(
          r_name, AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_type::OWNED);
        CHECK_CLOSE_SUMALL(3 * 3, n_ents, *mesh->get_comm());

        AmanziMesh::Entity_ID_List ents;
        mesh->get_set_entities(
          r_name, AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_type::ALL, &ents);
        for (const auto& e : ents) {
          auto fc = mesh->face_centroid(e);
          std::cout << r_name << " face: " << fc << std::endl;
          CHECK_CLOSE(0, fc[0], 1.e-10);
        }
      }

    } else if (r_name == "Central Face Box") {
      if (r_name == "Central Face Box") {
        // all should support this type
        CHECK(mesh->valid_set_type(r->get_type(), AmanziMesh::Entity_kind::FACE));
      }
      if (!mesh->valid_set_type(r->get_type(), AmanziMesh::Entity_kind::FACE)) continue;

      // a single face in the middle bottom
      int n_ents =
        mesh->get_set_size(r_name, AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_type::OWNED);
      CHECK_CLOSE_SUMALL(1, n_ents, *mesh->get_comm());

      AmanziMesh::Entity_ID_List ents;
      mesh->get_set_entities(
        r_name, AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_type::ALL, &ents);
      for (const auto& e : ents) {
        auto fc = mesh->face_centroid(e);
        CHECK_CLOSE(0.5, fc[0], 1.e-10);
        CHECK_CLOSE(0.5, fc[1], 1.e-10);
        CHECK_CLOSE(1, fc[2], 1.e-10);
      }

    } else if (r_name == "Face 106 - Central Face Box") {
      if (!mesh->valid_set_type(r->get_type(), AmanziMesh::Entity_kind::FACE)) continue;

      if (labeled) {
        // all but the above box
        int n_ents = mesh->get_set_size(
          r_name, AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_type::OWNED);
        CHECK_CLOSE_SUMALL(3 * 3 - 1, n_ents, *mesh->get_comm());

        AmanziMesh::Entity_ID_List ents;
        mesh->get_set_entities(
          r_name, AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_type::ALL, &ents);
        for (const auto& e : ents) {
          auto fc = mesh->face_centroid(e);
          std::cout << r_name << " centroid = " << fc << std::endl;
          CHECK_CLOSE(1, fc[2], 1.e-10);
          if (std::abs(fc[0] - 0.5) < 1e-10 && std::abs(fc[1] - 0.5) < 1e-10) { CHECK(false); }
        }
      }

    } else if (r_name == "Domain Boundary") {
      if (!mesh->valid_set_type(r->get_type(), AmanziMesh::Entity_kind::FACE)) continue;

      // all faces on the boundary
      int n_ents =
        mesh->get_set_size(r_name, AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_type::OWNED);
      CHECK_CLOSE_SUMALL(3 * 3 * 6, n_ents, *mesh->get_comm());

      AmanziMesh::Entity_ID_List ents;
      mesh->get_set_entities(
        r_name, AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_type::ALL, &ents);
      for (const auto& e : ents) {
        auto fc = mesh->face_centroid(e);
        CHECK(std::abs(fc[0]) < 1e-10 || std::abs(fc[1]) < 1e-10 || std::abs(fc[2]) < 1e-10 ||
              std::abs(fc[0] - 1) < 1e-10 || std::abs(fc[1] - 1) < 1e-10 ||
              std::abs(fc[2] - 1) < 1e-10);
      }

    } else if (r_name == "Sample Point OnFace") {
      if (!mesh->valid_set_type(r->get_type(), AmanziMesh::Entity_kind::CELL)) continue;

      // // two cells -- two cells share the face
      // // does not work -- implementation currently assumes exact arithmetic
      // int n_ents = mesh->get_set_size(r_name, AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_type::OWNED);
      // CHECK_CLOSE_SUMALL(2, n_ents, *mesh->get_comm());

      // AmanziMesh::Entity_ID_List ents;
      // mesh->get_set_entities(r_name, AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_type::ALL, &ents);
      // for (const auto& e : ents) {
      //   auto cc = mesh->cell_centroid(e);
      //   CHECK_CLOSE(0.5, cc[1], 1.e-10);
      //   CHECK_CLOSE(0.5, cc[2], 1.e-10);

      //   CHECK(std::abs(0.5 - cc[0]) < 1.e-10 ||
      //         std::abs(0.5/3 - cc[0]) < 1.e-10);

      // }

    } else if (r_name == "Sample Point OnEdge") {
      if (!mesh->valid_set_type(r->get_type(), AmanziMesh::Entity_kind::CELL)) continue;

      // // four cells share the edge
      // // does not work -- implementation currently assumes exact arithmetic
      // int n_ents = mesh->get_set_size(r_name, AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_type::OWNED);
      // CHECK_CLOSE_SUMALL(4, n_ents, *mesh->get_comm());

      // AmanziMesh::Entity_ID_List ents;
      // mesh->get_set_entities(r_name, AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_type::ALL, &ents);
      // for (const auto& e : ents) {
      //   auto cc = mesh->cell_centroid(e);
      //   CHECK_CLOSE(0.5, cc[2], 1.e-10);
      //   CHECK(std::abs(0.5 - cc[0]) < 1.e-10 ||
      //         std::abs(0.5/3 - cc[0]) < 1.e-10);
      //   CHECK(std::abs(0.5 - cc[1]) < 1.e-10 ||
      //         std::abs(0.5/3 - cc[1]) < 1.e-10);
      // }

    } else if (r_name == "Sample Point OnVertex") {
      if (!mesh->valid_set_type(r->get_type(), AmanziMesh::Entity_kind::CELL)) continue;

      // // eight cells share the point
      // // does not work -- implementation currently assumes exact arithmetic
      // int n_ents = mesh->get_set_size(r_name, AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_type::OWNED);
      // CHECK_CLOSE_SUMALL(8, n_ents, *mesh->get_comm());

      // AmanziMesh::Entity_ID_List ents;
      // mesh->get_set_entities(r_name, AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_type::ALL, &ents);
      // for (const auto& e : ents) {
      //   auto cc = mesh->cell_centroid(e);

      //   CHECK(std::abs(0.5 - cc[0]) < 1.e-10 ||
      //         std::abs(0.5/3 - cc[0]) < 1.e-10);
      //   CHECK(std::abs(0.5 - cc[1]) < 1.e-10 ||
      //         std::abs(0.5/3 - cc[1]) < 1.e-10);
      //   CHECK(std::abs(0.5 - cc[2]) < 1.e-10 ||
      //         std::abs(0.5/3 - cc[2]) < 1.e-10);
      // }

    } else if (r_name == "Interior XY Plane") {
      CHECK(mesh->valid_set_type(r->get_type(), AmanziMesh::Entity_kind::NODE));

      int n_ents =
        mesh->get_set_size(r_name, AmanziMesh::Entity_kind::NODE, AmanziMesh::Parallel_type::OWNED);
      CHECK_CLOSE_SUMALL(4 * 4, n_ents, *mesh->get_comm());

      AmanziMesh::Entity_ID_List ents;
      mesh->get_set_entities(
        r_name, AmanziMesh::Entity_kind::NODE, AmanziMesh::Parallel_type::ALL, &ents);
      for (const auto& e : ents) {
        AmanziGeometry::Point nc;
        mesh->node_get_coordinates(e, &nc);
        CHECK_CLOSE(1.0 / 3, nc[2], 1.e-10);
      }

    } else if (r_name == "Top Box Nodes") {
      CHECK(mesh->valid_set_type(r->get_type(), AmanziMesh::Entity_kind::NODE));

      int n_ents =
        mesh->get_set_size(r_name, AmanziMesh::Entity_kind::NODE, AmanziMesh::Parallel_type::OWNED);
      CHECK_CLOSE_SUMALL(4, n_ents, *mesh->get_comm());

      AmanziMesh::Entity_ID_List ents;
      mesh->get_set_entities(
        r_name, AmanziMesh::Entity_kind::NODE, AmanziMesh::Parallel_type::ALL, &ents);
      for (const auto& e : ents) {
        AmanziGeometry::Point nc;
        mesh->node_get_coordinates(e, &nc);
        if (nc[0] > 0.5)
          CHECK_CLOSE(2.0 / 3, nc[0], 1.e-10);
        else
          CHECK_CLOSE(1.0 / 3, nc[0], 1.e-10);
        if (nc[1] > 0.5)
          CHECK_CLOSE(2.0 / 3, nc[1], 1.e-10);
        else
          CHECK_CLOSE(1.0 / 3, nc[1], 1.e-10);
        CHECK_CLOSE(1.0, nc[2], 1.e-10);
      }

    } else if (r_name == "Unit Hex") {
      // pass -- this is used for extraction method
      continue;
    } else {
      std::cout << "Bad set name: \"" << r_name << "\"" << std::endl;
      CHECK(false);
    }
  }
}


//
// Helper functions for testing sets
//

//
// Tests Mesh sets on a 3x3x3 hex
//
template <class Mesh_type>
void
testQuadMeshSets3x3(const Teuchos::RCP<Mesh_type>& mesh,
                    bool labeled,
                    Framework f,
                    bool extracted,
                    bool broken_labeled_set_inheritance = false)
{
  auto gm = mesh->geometric_model();
  auto num_regions = gm->size();
  CHECK(num_regions == 23 || num_regions == 24);
  CHECK_EQUAL(2, mesh->space_dimension());
  CHECK_EQUAL(2, mesh->manifold_dimension());

  for (const auto& r : *gm) {
    std::string r_name = r->get_name();
    std::cout << "Checking set: " << r->get_name() << std::endl;

    if (r_name == "Entire Mesh") {
      if (!mesh->valid_set_type(r->get_type(), AmanziMesh::Entity_kind::CELL)) continue;

      // check that the sum of n_ents is the global number of cells, 3*3
      int n_ents =
        mesh->get_set_size(r_name, AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_type::OWNED);
      CHECK_CLOSE_SUMALL(3 * 3, n_ents, *mesh->get_comm());

      // check that the number of entities is the number of owned cells
      CHECK_EQUAL(
        mesh->num_entities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_type::OWNED),
        n_ents);
      // check that the number of entities is the number of all cells
      CHECK_EQUAL(
        mesh->num_entities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_type::ALL),
        mesh->get_set_size(r_name, AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_type::ALL));


    } else if (r_name == "Top LS" || r_name == "Cell Set 3" || r_name == "Top ColFunc" ||
               r_name == "Face 106" || r_name == "Top Face Plane") {
      if (!extracted) continue;                             // only valid for extracted
      if (!labeled && r_name != "Top Face Plane") continue; // only valid for exo meshes
      if (r_name == "Top Face Plane") continue;             // This only used to work due to a bug!
      if (!mesh->valid_set_type(r->get_type(), AmanziMesh::Entity_kind::CELL)) continue;
      if (r_name == "Cell Set 3" && f == Framework::MOAB)
        continue; // MOAB cannot handle cell sets, only mat IDs
      if (r_name == "Cell Set 3" && f == Framework::MSTK)
        continue; // MSTK does not correctly deal with surface cells inferred from the cell below the top surface
      if (r_name == "Top ColFunc") continue; // inferring cell-to-cell is not implemented
      if (r_name == "Top LS") continue;      // inferring cell-to-cell is not implemented

      // this is all cells
      int n_ents =
        mesh->get_set_size(r_name, AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_type::OWNED);
      CHECK_CLOSE_SUMALL(3 * 3, n_ents, *mesh->get_comm());

      AmanziMesh::Entity_ID_List ents;
      mesh->get_set_entities(
        r_name, AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_type::ALL, &ents);

    } else if (r_name == "Box") {
      // all should support this type
      CHECK(mesh->valid_set_type(r->get_type(), AmanziMesh::Entity_kind::CELL));

      // left-most layer (or west-most?)
      int n_ents =
        mesh->get_set_size(r_name, AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_type::OWNED);
      CHECK_CLOSE_SUMALL(3, n_ents, *mesh->get_comm());

      AmanziMesh::Entity_ID_List ents;
      mesh->get_set_entities(
        r_name, AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_type::ALL, &ents);
      for (const auto& e : ents) {
        auto cc = mesh->cell_centroid(e);
        CHECK_CLOSE(0.5 / 3, cc[1], 1.e-10); // left
      }

    } else if (r_name == "NOT_Box" || r_name == "Entire Mesh - Box") {
      if (!mesh->valid_set_type(r->get_type(), AmanziMesh::Entity_kind::CELL)) continue;

      // others
      int n_ents =
        mesh->get_set_size(r_name, AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_type::OWNED);
      CHECK_CLOSE_SUMALL(2 * 3, n_ents, *mesh->get_comm());

      AmanziMesh::Entity_ID_List ents;
      mesh->get_set_entities(
        r_name, AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_type::ALL, &ents);
      for (const auto& e : ents) {
        auto cc = mesh->cell_centroid(e);
        if (cc[1] > 0.6)
          CHECK_CLOSE(2.5 / 3, cc[1], 1.e-10); // left
        else
          CHECK_CLOSE(1.5 / 3, cc[1], 1.e-10); // left
      }

    } else if (r_name == "Point" || r_name == "Sample Point InCell" ||
               r_name == "Point intersects NOT_Box") {
      if (!mesh->valid_set_type(r->get_type(), AmanziMesh::Entity_kind::CELL)) continue;
      // the central-most cell

      int n_ents =
        mesh->get_set_size(r_name, AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_type::OWNED);
      CHECK_CLOSE_SUMALL(1, n_ents, *mesh->get_comm());

      AmanziMesh::Entity_ID_List ents;
      mesh->get_set_entities(
        r_name, AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_type::ALL, &ents);
      for (const auto& e : ents) {
        auto cc = mesh->cell_centroid(e);
        CHECK_CLOSE(0.5, cc[0], 1.e-10);
        CHECK_CLOSE(0.5, cc[1], 1.e-10);
      }

    } else if (r_name == "Box+Point") {
      if (!mesh->valid_set_type(r->get_type(), AmanziMesh::Entity_kind::CELL)) continue;

      int n_ents =
        mesh->get_set_size(r_name, AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_type::OWNED);
      CHECK_CLOSE_SUMALL(3 + 1, n_ents, *mesh->get_comm());

      AmanziMesh::Entity_ID_List ents;
      mesh->get_set_entities(
        r_name, AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_type::ALL, &ents);
      for (const auto& e : ents) {
        auto cc = mesh->cell_centroid(e);
        if (cc[1] > 0.33) {
          CHECK_CLOSE(0.5, cc[1], 1.e-10);
        } else {
          CHECK_CLOSE(0.5 / 3, cc[1], 1.e-10);
        }
      }

    } else if (r_name == "Enumerated") {
      if (!mesh->valid_set_type(r->get_type(), AmanziMesh::Entity_kind::CELL)) continue;

      // if (framework == Framework::MSTK) { // logical only works in MSTK
      // // cells 1,2,3 -- only supported by MSTK
      // NOT WORKING!
      // int n_ents = mesh->get_set_size(r_name, AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_type::OWNED);
      // CHECK_CLOSE_SUMALL(3, n_ents, *mesh->get_comm());

      // AmanziMesh::Entity_ID_List ents;
      // mesh->get_set_entities(r_name, AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_type::ALL, &ents);
      // std::vector<int> ents_found(3, 0);
      // for (const auto& e : ents) {
      //   int gid = mesh->cell_map(false).GID(e);
      //   if (gid == 1) ents_found[0] = 1;
      //   if (gid == 2) ents_found[1] = 1;
      //   if (gid == 3) ents_found[2] = 1;
      // }
      // CHECK_MPI_ALL(ents_found, *mesh->get_comm());
      // }

    } else if (r_name == "Face 103" || r_name == "Side Plane") {
      if (!mesh->valid_set_type(r->get_type(), AmanziMesh::Entity_kind::FACE)) continue;

      // Face 103 is not correct in MeshExtractedManifold -- inheritance from
      // the subsurface does not seem to be quite right here.
      if (broken_labeled_set_inheritance && r_name == "Face 103") continue;

      // the top faces
      if (labeled || r_name == "Side Plane") {
        int n_ents = mesh->get_set_size(
          r_name, AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_type::OWNED);
        CHECK_CLOSE_SUMALL(3, n_ents, *mesh->get_comm());

        AmanziMesh::Entity_ID_List ents;
        mesh->get_set_entities(
          r_name, AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_type::ALL, &ents);
        for (const auto& e : ents) {
          auto fc = mesh->face_centroid(e);
          std::cout << r_name << " face: " << fc << std::endl;
          CHECK_CLOSE(0, fc[0], 1.e-10);
        }
      }

    } else if (r_name == "Central Face Box") {
      CHECK(mesh->valid_set_type(r->get_type(), AmanziMesh::Entity_kind::FACE));

      // a single face in the middle bottom
      int n_ents =
        mesh->get_set_size(r_name, AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_type::OWNED);
      CHECK_CLOSE_SUMALL(1, n_ents, *mesh->get_comm());

      AmanziMesh::Entity_ID_List ents;
      mesh->get_set_entities(
        r_name, AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_type::ALL, &ents);
      for (const auto& e : ents) {
        auto fc = mesh->face_centroid(e);
        CHECK_CLOSE(0.5, fc[1], 1.e-10);
        CHECK_CLOSE(0., fc[0], 1.e-10);
      }

    } else if (r_name == "Side Plane - Central Face Box" ||
               r_name == "Face 103 - Central Face Box") {
      if (!mesh->valid_set_type(r->get_type(), AmanziMesh::Entity_kind::FACE)) continue;

      // Face 103 is not correct in MeshExtractedManifold -- inheritance from
      // the subsurface does not seem to be quite right here.
      if (broken_labeled_set_inheritance && r_name == "Face 103 - Central Face Box") continue;

      if (labeled || r_name == "Side Plane - Central Face Box") {
        // all but the above box
        int n_ents = mesh->get_set_size(
          r_name, AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_type::OWNED);
        CHECK_CLOSE_SUMALL(3 - 1, n_ents, *mesh->get_comm());

        AmanziMesh::Entity_ID_List ents;
        mesh->get_set_entities(
          r_name, AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_type::ALL, &ents);
        for (const auto& e : ents) {
          auto fc = mesh->face_centroid(e);
          std::cout << r_name << " centroid = " << fc << std::endl;
          CHECK_CLOSE(0, fc[0], 1.e-10);
          if (fc[1] > 0.5)
            CHECK_CLOSE(2.5 / 3, fc[1], 1.e-10);
          else
            CHECK_CLOSE(0.5 / 3, fc[1], 1.e-1);
        }
      }

    } else if (r_name == "Domain Boundary") {
      if (!mesh->valid_set_type(r->get_type(), AmanziMesh::Entity_kind::FACE)) continue;

      // all faces on the boundary
      int n_ents =
        mesh->get_set_size(r_name, AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_type::OWNED);
      CHECK_CLOSE_SUMALL(3 * 4, n_ents, *mesh->get_comm());

      AmanziMesh::Entity_ID_List ents;
      mesh->get_set_entities(
        r_name, AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_type::ALL, &ents);
      for (const auto& e : ents) {
        auto fc = mesh->face_centroid(e);
        CHECK(std::abs(fc[0]) < 1e-10 || std::abs(fc[1]) < 1e-10 || std::abs(fc[0] - 1) < 1e-10 ||
              std::abs(fc[1] - 1) < 1e-10);
      }

    } else if (r_name == "Sample Point OnFace") {
      if (!mesh->valid_set_type(r->get_type(), AmanziMesh::Entity_kind::CELL)) continue;

      // // two cells -- two cells share the face
      // // does not work -- implementation currently assumes exact arithmetic
      // int n_ents = mesh->get_set_size(r_name, AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_type::OWNED);
      // CHECK_CLOSE_SUMALL(2, n_ents, *mesh->get_comm());

      // AmanziMesh::Entity_ID_List ents;
      // mesh->get_set_entities(r_name, AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_type::ALL, &ents);
      // for (const auto& e : ents) {
      //   auto cc = mesh->cell_centroid(e);
      //   CHECK_CLOSE(0.5, cc[1], 1.e-10);
      //   CHECK_CLOSE(0.5, cc[2], 1.e-10);

      //   CHECK(std::abs(0.5 - cc[0]) < 1.e-10 ||
      //         std::abs(0.5/3 - cc[0]) < 1.e-10);

      // }

    } else if (r_name == "Sample Point OnVertex") {
      if (!mesh->valid_set_type(r->get_type(), AmanziMesh::Entity_kind::CELL)) continue;

      // // four cells share the edge
      // // does not work -- implementation currently assumes exact arithmetic
      // int n_ents = mesh->get_set_size(r_name, AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_type::OWNED);
      // CHECK_CLOSE_SUMALL(4, n_ents, *mesh->get_comm());

      // AmanziMesh::Entity_ID_List ents;
      // mesh->get_set_entities(r_name, AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_type::ALL, &ents);
      // for (const auto& e : ents) {
      //   auto cc = mesh->cell_centroid(e);
      //   CHECK_CLOSE(0.5, cc[2], 1.e-10);
      //   CHECK(std::abs(0.5 - cc[0]) < 1.e-10 ||
      //         std::abs(0.5/3 - cc[0]) < 1.e-10);
      //   CHECK(std::abs(0.5 - cc[1]) < 1.e-10 ||
      //         std::abs(0.5/3 - cc[1]) < 1.e-10);
      // }

    } else if (r_name == "Interior XY Plane") {
      CHECK(mesh->valid_set_type(r->get_type(), AmanziMesh::Entity_kind::NODE));

      int n_ents =
        mesh->get_set_size(r_name, AmanziMesh::Entity_kind::NODE, AmanziMesh::Parallel_type::OWNED);
      CHECK_CLOSE_SUMALL(4, n_ents, *mesh->get_comm());

      AmanziMesh::Entity_ID_List ents;
      mesh->get_set_entities(
        r_name, AmanziMesh::Entity_kind::NODE, AmanziMesh::Parallel_type::ALL, &ents);
      for (const auto& e : ents) {
        AmanziGeometry::Point nc;
        mesh->node_get_coordinates(e, &nc);
        CHECK_CLOSE(1.0 / 3, nc[1], 1.e-10);
      }

    } else if (r_name == "Top Box Nodes") {
      CHECK(mesh->valid_set_type(r->get_type(), AmanziMesh::Entity_kind::NODE));

      int n_ents =
        mesh->get_set_size(r_name, AmanziMesh::Entity_kind::NODE, AmanziMesh::Parallel_type::OWNED);
      CHECK_CLOSE_SUMALL(4, n_ents, *mesh->get_comm());

      AmanziMesh::Entity_ID_List ents;
      mesh->get_set_entities(
        r_name, AmanziMesh::Entity_kind::NODE, AmanziMesh::Parallel_type::ALL, &ents);
      for (const auto& e : ents) {
        AmanziGeometry::Point nc;
        mesh->node_get_coordinates(e, &nc);
        if (nc[0] > 0.5)
          CHECK_CLOSE(2.0 / 3, nc[0], 1.e-10);
        else
          CHECK_CLOSE(1.0 / 3, nc[0], 1.e-10);
        if (nc[1] > 0.5)
          CHECK_CLOSE(2.0 / 3, nc[1], 1.e-10);
        else
          CHECK_CLOSE(1.0 / 3, nc[1], 1.e-10);
      }
    } else {
      std::cout << "Bad set name: \"" << r_name << "\"" << std::endl;
      CHECK(false);
    }
  }
}
