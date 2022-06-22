/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:

*/

//!

#include <fstream>

#include "AmanziMap.hh"
#include "AmanziComm.hh"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "UnitTest++.h"

#include "../Mesh_MSTK.hh"


TEST(MSTK_READ_NONMANIFOLD_SURFACES)
{
  std::string expcsetnames[2] = { "FRACTURE 1", "FRACTURE 2" };

  auto comm = Amanzi::getDefaultComm();

  Teuchos::ParameterList parameterlist;

  // create a sublist name Regions and put a reference to it in
  // reg_spec and other sublists as references. Turns out it is
  // important to define reg_spec and other lists below as references
  // - otherwise, a new copy is made of the sublist that is retrieved

  Teuchos::ParameterList& reg_spec = parameterlist.sublist("regions");

  Teuchos::ParameterList& frac1 = reg_spec.sublist("FRACTURE 1");
  Teuchos::ParameterList& frac1_def = frac1.sublist("region: plane");
  Teuchos::Array<double> loc1 = Teuchos::tuple(0.0, 0.0, 0.5);
  Teuchos::Array<double> dir1 = Teuchos::tuple(0.0, 0.0, 1.0);
  frac1_def.set<Teuchos::Array<double>>("point", loc1);
  frac1_def.set<Teuchos::Array<double>>("normal", dir1);
  frac1_def.set<std::string>("entity", "cell");

  Teuchos::ParameterList& frac2 = reg_spec.sublist("FRACTURE 2");
  Teuchos::ParameterList& frac2_def = frac2.sublist("region: plane");
  Teuchos::Array<double> loc2 = Teuchos::tuple(0.0, 0.5, 0.0);
  Teuchos::Array<double> dir2 = Teuchos::tuple(0.0, -1.0, 0.0);
  frac2_def.set<Teuchos::Array<double>>("point", loc2);
  frac2_def.set<Teuchos::Array<double>>("normal", dir2);
  frac2_def.set<std::string>("entity", "cell");


  Teuchos::RCP<Amanzi::AmanziGeometry::GeometricModel> gm = Teuchos::rcp(
    new Amanzi::AmanziGeometry::GeometricModel(3, reg_spec, *comm));


  // Load a mesh consisting of quads representing two perpendicular
  // intersecting surfaces in 3 space. Each surface is meshed with 100
  // regular quads

  Teuchos::RCP<Amanzi::AmanziMesh::Mesh> mesh(
    new Amanzi::AmanziMesh::Mesh_MSTK("test/fractures.exo", comm, 3, gm));


  CHECK_EQUAL(2, mesh->manifold_dimension());
  CHECK_EQUAL(3, mesh->space_dimension());

  int nc = mesh->num_entities(Amanzi::AmanziMesh::CELL,
                              Amanzi::AmanziMesh::Parallel_type::OWNED);

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
      if (reg_name == "FRACTURE 1" || reg_name == "FRACTURE 2") {
        // Do we have a valid cellset by this name

        CHECK(mesh->valid_set_name(reg_name, Amanzi::AmanziMesh::CELL));

        if (expcsetnames[0] != reg_name && expcsetnames[1] != reg_name) break;

        // Verify that we can get the number of entities in the set

        int set_size =
          mesh->get_set_size(reg_name,
                             Amanzi::AmanziMesh::CELL,
                             Amanzi::AmanziMesh::Parallel_type::OWNED);

        CHECK_EQUAL(nc / 2, set_size);

        // Verify that we can retrieve the set entities

        Amanzi::AmanziMesh::Entity_ID_List setents;
        mesh->get_set_entities(reg_name,
                               Amanzi::AmanziMesh::CELL,
                               Amanzi::AmanziMesh::Parallel_type::OWNED,
                               &setents);

        if (reg_name == "FRACTURE 1") {
          // VERIFY THAT EACH 2D CELL HAS A NORMAL THAT POINTS IN [0,0,1] DIR
          std::vector<Amanzi::AmanziGeometry::Point> ccoords;
          for (auto const& cellid : setents) {
            mesh->cell_get_coordinates(cellid, &ccoords);

            // Assumes that for 2D cells the coordinates are returned in
            // ccw direction around the cell
            Amanzi::AmanziGeometry::Point vec0 = ccoords[2] - ccoords[1];
            Amanzi::AmanziGeometry::Point vec1 = ccoords[0] - ccoords[1];
            Amanzi::AmanziGeometry::Point cpvec = vec0 ^ vec1;
            double len = norm(cpvec);
            cpvec /= len;

            CHECK_CLOSE(0.0, cpvec[0], 1.0e-08);
            CHECK_CLOSE(0.0, cpvec[1], 1.0e-08);
            CHECK_CLOSE(1.0, cpvec[2], 1.0e-08);
          }
        } else {
          // VERIFY THAT EACH 2D CELL HAS A NORMAL THAT POINTS IN [0,-1,0] DIR
          std::vector<Amanzi::AmanziGeometry::Point> ccoords;
          for (auto const& cellid : setents) {
            mesh->cell_get_coordinates(cellid, &ccoords);

            // Assumes that for 2D cells the coordinates are returned in
            // ccw direction around the cell
            Amanzi::AmanziGeometry::Point vec0 = ccoords[2] - ccoords[1];
            Amanzi::AmanziGeometry::Point vec1 = ccoords[0] - ccoords[1];
            Amanzi::AmanziGeometry::Point cpvec = vec0 ^ vec1;
            double len = norm(cpvec);
            cpvec /= len;

            CHECK_CLOSE(0.0, cpvec[0], 1.0e-08);
            CHECK_CLOSE(-1.0, cpvec[1], 1.0e-08);
            CHECK_CLOSE(0.0, cpvec[2], 1.0e-08);
          }
        }
      }
    }
  }


  // Check some things about the non-manifold edges where the two
  // surfaces come together

  int nf = mesh->num_entities(Amanzi::AmanziMesh::FACE,
                              Amanzi::AmanziMesh::Parallel_type::OWNED);
  for (int f = 0; f < nf; f++) {
    Amanzi::AmanziMesh::Entity_ID_List fcells;
    mesh->face_get_cells_host(f, Amanzi::AmanziMesh::Parallel_type::ALL, &fcells);

    std::vector<Amanzi::AmanziGeometry::Point> fnormals(fcells.size());

    Amanzi::AmanziGeometry::Point fcen = mesh->face_centroid(f);
    for (int j = 0; j < fcells.size(); j++) {
      int dir;
      fnormals[j] = mesh->face_normal(f, false, fcells[j], &dir);
      fnormals[j] /= norm(fnormals[j]);

      // Check that the face normal with respect to the cell is
      // actually pointing out of the face

      Amanzi::AmanziGeometry::Point ccen = mesh->cell_centroid(fcells[j]);
      Amanzi::AmanziGeometry::Point nvec = fcen - ccen;
      nvec /= norm(nvec);

      double dp = nvec * fnormals[j];
      CHECK_CLOSE(1.0, dp, 1.0e-12);
    }

    if (fcells.size() == 2) { // Is fnormals[0] = -fnormals[1]?
      for (int d = 0; d < 3; d++)
        CHECK_CLOSE(fnormals[0][d], -fnormals[1][d], 1.0e-12);
    }
  }
}
