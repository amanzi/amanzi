/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:

*/

//!

#include <iostream>
#include "stdlib.h"
#include "math.h"
#include "UnitTest++.h"
#include "../Mesh_simple.hh"
#include "AmanziComm.hh"
#include "GenerationSpec.hh"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"

SUITE(MeshSimple)
{
  TEST(SETS)
  {
    using namespace std;

    auto comm = Amanzi::getDefaultComm();

    std::string expcsetnames[6] = { "Bottom Box",     "Bottom+Middle Box",
                                    "Vertical Box",   "Bottom ColFunc",
                                    "Middle ColFunc", "Top ColFunc" };
    unsigned int csetsize, expcsetsizes[6] = { 9, 18, 9, 9, 9, 9 };

    unsigned int expcsetcells[6][18] = {
      { 0, 3, 6, 1, 4, 7, 2, 5, 8, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
      { 0, 9, 3, 12, 6, 15, 1, 10, 4, 13, 7, 16, 2, 11, 5, 14, 8, 17 },
      { 1, 10, 19, 4, 13, 22, 7, 16, 25, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
      { 0, 3, 6, 1, 4, 7, 2, 5, 8, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
      { 9, 12, 15, 10, 13, 16, 11, 14, 17, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
      { 18, 21, 24, 19, 22, 25, 20, 23, 26, 0, 0, 0, 0, 0, 0, 0, 0, 0 }
    };

    std::string expfsetnames[2] = { "ZLO FACE Plane", "YLO FACE Box" };

    unsigned int fsetsize, expfsetsizes[2] = { 9, 3 };


    int expfsetfaces[2][9] = { { 0, 3, 6, 1, 4, 7, 2, 5, 8 },
                               { 60, 61, 62, -1, -1, -1, -1, -1, -1 } };

    std::string expnsetnames[3] = { "Sample Point 1",
                                    "INTERIOR XY PLANE",
                                    "TOP BOX" };

    unsigned int nsetsize, expnsetsizes[3] = { 1, 16, 4 };
    int expnsetnodes[3][16] = {
      { 21, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
      { 16, 20, 24, 28, 17, 21, 25, 29, 18, 22, 26, 30, 19, 23, 27, 31 },
      { 53, 57, 54, 58, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 }
    };

    std::string infilename = "test/hex_3x3x3.xml";
    Teuchos::ParameterXMLFileReader xmlreader(infilename);

    Teuchos::ParameterList reg_spec(xmlreader.getParameters());

    Teuchos::RCP<Amanzi::AmanziGeometry::GeometricModel> gm = Teuchos::rcp(
      new Amanzi::AmanziGeometry::GeometricModel(3, reg_spec, *comm));

    // Create a mesh consisting of 3x3x3 elements (4x4x4 nodes)

    Amanzi::AmanziMesh::Mesh_simple mesh(
      0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 3, 3, 3, comm, gm);

    Teuchos::ParameterList::ConstIterator i;
    for (i = reg_spec.begin(); i != reg_spec.end(); i++) {
      const std::string reg_name = reg_spec.name(i);

      Teuchos::ParameterList reg_params = reg_spec.sublist(reg_name);

      // See if the geometric model has a region by this name
      Teuchos::RCP<const Amanzi::AmanziGeometry::Region> reg =
        gm->FindRegion(reg_name);

      CHECK(reg.get());

      // Do their names match ?

      CHECK_EQUAL(reg->name(), reg_name);


      // Get the region info directly from the XML and compare

      Teuchos::ParameterList::ConstIterator j = reg_params.begin();

      std::string shape = reg_params.name(j);

      if (shape == "region: plane") {
        if (reg_name == "ZLO FACE Plane") {
          // Do we have a valid sideset by this name

          CHECK(mesh.valid_set_name(reg_name, Amanzi::AmanziMesh::FACE));

          int j;
          for (j = 0; j < 2; j++) {
            if (expfsetnames[j] == reg_name) break;
          }

          CHECK(j < 2);


          // Verify that we can get the right number of entities in the set

          int set_size =
            mesh.get_set_size(reg_name,
                              Amanzi::AmanziMesh::FACE,
                              Amanzi::AmanziMesh::Parallel_type::ALL);

          CHECK_EQUAL(expfsetsizes[j], set_size);


          // Verify that we can get the correct set entities

          Amanzi::AmanziMesh::Entity_ID_List setents;
          mesh.get_set_entities(reg_name,
                                Amanzi::AmanziMesh::FACE,
                                Amanzi::AmanziMesh::Parallel_type::ALL,
                                setents);

          CHECK_ARRAY_EQUAL(expfsetfaces[j], setents, set_size);

        } else if (reg_name == "INTERIOR XY PLANE") {
          // Do we have a valid nodeset by this name

          CHECK(mesh.valid_set_name(reg_name, Amanzi::AmanziMesh::NODE));

          int j;
          for (j = 0; j < 3; j++) {
            if (expnsetnames[j] == reg_name) break;
          }

          CHECK(j < 3);


          // Verify that we can get the right number of entities in the set

          int set_size =
            mesh.get_set_size(reg_name,
                              Amanzi::AmanziMesh::NODE,
                              Amanzi::AmanziMesh::Parallel_type::ALL);

          CHECK_EQUAL(expnsetsizes[j], set_size);


          // Verify that we can get the correct set entities

          Amanzi::AmanziMesh::Entity_ID_List setents;
          mesh.get_set_entities(reg_name,
                                Amanzi::AmanziMesh::NODE,
                                Amanzi::AmanziMesh::Parallel_type::ALL,
                                setents);

          CHECK_ARRAY_EQUAL(expnsetnodes[j], setents, set_size);
        }

      } else if (shape == "region: box") {
        Teuchos::ParameterList box_params = reg_params.sublist(shape);
        Teuchos::Array<double> pmin =
          box_params.get<Teuchos::Array<double>>("low coordinate");
        Teuchos::Array<double> pmax =
          box_params.get<Teuchos::Array<double>>("high coordinate");

        if (pmin[0] == pmax[0] || pmin[1] == pmax[1] || pmin[2] == pmax[2]) {
          // This is a reduced dimensionality box - request a faceset or nodeset

          if (reg_name == "YLO FACE BOX") {
            // Do we have a valid sideset by this name

            CHECK(mesh.valid_set_name(reg_name, Amanzi::AmanziMesh::FACE));

            int j;
            for (j = 0; j < 2; j++) {
              if (expfsetnames[j] == reg_name) break;
            }

            CHECK(j < 2);


            // Verify that we can get the right number of entities in the set

            int set_size =
              mesh.get_set_size(reg_name,
                                Amanzi::AmanziMesh::FACE,
                                Amanzi::AmanziMesh::Parallel_type::ALL);

            CHECK_EQUAL(expfsetsizes[j], set_size);


            // Verify that we can get the correct set entities

            Amanzi::AmanziMesh::Entity_ID_List setents;
            mesh.get_set_entities(reg_name,
                                  Amanzi::AmanziMesh::FACE,
                                  Amanzi::AmanziMesh::Parallel_type::ALL,
                                  setents);

            CHECK_ARRAY_EQUAL(expfsetfaces[j], setents, set_size);

          } else if (reg_name == "TOP BOX") {
            // Do we have a valid set by this name

            CHECK(mesh.valid_set_name(reg_name, Amanzi::AmanziMesh::NODE));

            int j;
            for (j = 0; j < 3; j++) {
              if (expnsetnames[j] == reg_name) break;
            }

            CHECK(j < 3);


            // Verify that we can get the right number of entities in the set

            int set_size =
              mesh.get_set_size(reg_name,
                                Amanzi::AmanziMesh::NODE,
                                Amanzi::AmanziMesh::Parallel_type::ALL);

            CHECK_EQUAL(expnsetsizes[j], set_size);


            // Verify that we can get the correct set entities

            Amanzi::AmanziMesh::Entity_ID_List setents;
            mesh.get_set_entities(reg_name,
                                  Amanzi::AmanziMesh::NODE,
                                  Amanzi::AmanziMesh::Parallel_type::ALL,
                                  setents);

            CHECK_ARRAY_EQUAL(expnsetnodes[j], setents, set_size);
          }
        } else {
          // Do we have a valid cellset by this name

          CHECK(mesh.valid_set_name(reg_name, Amanzi::AmanziMesh::CELL));

          // Find the expected cell set info corresponding to this name

          int j;
          for (j = 0; j < 3; j++)
            if (reg_name == expcsetnames[j]) break;

          CHECK(j < 3);

          // Verify that we can get the right number of entities in the set

          int set_size =
            mesh.get_set_size(reg_name,
                              Amanzi::AmanziMesh::CELL,
                              Amanzi::AmanziMesh::Parallel_type::ALL);

          CHECK_EQUAL(expcsetsizes[j], set_size);

          // Verify that we can get the correct set entities

          Amanzi::AmanziMesh::Entity_ID_List setents;
          mesh.get_set_entities(reg_name,
                                Amanzi::AmanziMesh::CELL,
                                Amanzi::AmanziMesh::Parallel_type::ALL,
                                setents);

          CHECK_ARRAY_EQUAL(expcsetcells[j], setents, set_size);
        }
      } else if (shape == "region: point") {
        // Do we have a valid set by this name

        CHECK(mesh.valid_set_name(reg_name, Amanzi::AmanziMesh::NODE));

        int j;
        for (j = 0; j < 3; j++) {
          if (expnsetnames[j] == reg_name) break;
        }

        CHECK(j < 3);


        // Verify that we can get the right number of entities in the set

        int set_size =
          mesh.get_set_size(reg_name,
                            Amanzi::AmanziMesh::NODE,
                            Amanzi::AmanziMesh::Parallel_type::ALL);

        CHECK_EQUAL(expnsetsizes[j], set_size);


        // Verify that we can get the correct set entities

        Amanzi::AmanziMesh::Entity_ID_List setents;
        mesh.get_set_entities(reg_name,
                              Amanzi::AmanziMesh::NODE,
                              Amanzi::AmanziMesh::Parallel_type::ALL,
                              setents);

        CHECK_ARRAY_EQUAL(expnsetnodes[j], setents, set_size);

      } else if (shape == "region: labeled set") {
        std::cerr << "Mesh framework cannot do labeled sets" << std::endl;

        CHECK(false);

      } else if (shape == "region: color function") {
        // Do we have a valid set by this name

        CHECK(mesh.valid_set_name(reg_name, Amanzi::AmanziMesh::CELL));

        int j;
        for (j = 0; j < 6; j++) {
          if (expcsetnames[j] == reg_name) break;
        }

        CHECK(j < 6);


        // Verify that we can get the right number of entities in the set

        int set_size =
          mesh.get_set_size(reg_name,
                            Amanzi::AmanziMesh::CELL,
                            Amanzi::AmanziMesh::Parallel_type::ALL);

        CHECK_EQUAL(expcsetsizes[j], set_size);


        // Verify that we can get the correct set entities

        Amanzi::AmanziMesh::Entity_ID_List setents;
        mesh.get_set_entities(reg_name,
                              Amanzi::AmanziMesh::CELL,
                              Amanzi::AmanziMesh::Parallel_type::ALL,
                              setents);

        CHECK_ARRAY_EQUAL(expcsetcells[j], setents, set_size);
      }
    }
  }
}
