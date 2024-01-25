/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Rao Garimella, others
*/

#include <iostream>
#include <UnitTest++.h>

#include <AmanziComm.hh>

#include "dbc.hh"
#include "GeometricModel.hh"
#include "Mesh.hh"
#include "../MeshFactory.hh"
#include "../MeshException.hh"

// -------------------------------------------------------------
// check_preference
// -------------------------------------------------------------

SUITE(MeshFramework)
{
  TEST(EXTRACT_COMM_SPLIT)
  {
    auto comm = Amanzi::getDefaultComm();

    double x0(0), y0(0), z0(0);
    double x1(10), y1(10), z1(10);
    int nx(10), ny(10), nz(10);

    // make regions
    Teuchos::ParameterList parameterlist;
    Teuchos::ParameterList& reg_spec = parameterlist.sublist("regions");
    Teuchos::ParameterList& top_surface = reg_spec.sublist("Top Surface");
    Teuchos::ParameterList& top_surface_def = top_surface.sublist("region: box");
    Teuchos::Array<double> low = Teuchos::tuple(0.0, 0.0, 10.0);
    Teuchos::Array<double> high = Teuchos::tuple(10.0, 10.0, 10.0);
    top_surface_def.set<Teuchos::Array<double>>("low coordinate", low);
    top_surface_def.set<Teuchos::Array<double>>("high coordinate", high);

    Teuchos::RCP<Amanzi::AmanziGeometry::GeometricModel> gm =
      Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(3, reg_spec, *comm));

    std::vector<std::string> setnames{ "Top Surface" };

    Amanzi::AmanziMesh::Preference pref;

    auto plist = Teuchos::rcp(new Teuchos::ParameterList());
    plist->set<std::string>("partitioner", "metis");
    plist->set<bool>("create subcommunicator", true);
    Amanzi::AmanziMesh::MeshFactory meshfactory(comm, gm, plist);

    if (Amanzi::AmanziMesh::framework_enabled(Amanzi::AmanziMesh::Framework::MSTK)) {
      pref.clear();
      pref.push_back(Amanzi::AmanziMesh::Framework::MSTK);
      meshfactory.set_preference(pref);

      auto parent_mesh = meshfactory.create(x0, y0, z0, x1, y1, z1, nx, ny, nz);
      CHECK(!parent_mesh.is_null());

      auto setents = parent_mesh->getSetEntities("Top Surface",
              Amanzi::AmanziMesh::Entity_kind::FACE,
              Amanzi::AmanziMesh::Parallel_kind::OWNED);
      int local_count = setents.size();
      int global_min_count;
      std::cout << "Count on rank " << comm->getRank() << " is " << local_count << std::endl;
      Teuchos::reduceAll(*comm, Teuchos::REDUCE_MIN, 1, &local_count, &global_min_count);
      AMANZI_ASSERT(global_min_count == 0); // this test only useful for min = 0

      int local_cells = parent_mesh->getNumEntities(Amanzi::AmanziMesh::Entity_kind::CELL,
                                                    Amanzi::AmanziMesh::Parallel_kind::OWNED);
      int global_cells;
      Teuchos::reduceAll(*comm, Teuchos::REDUCE_SUM, 1, &local_cells, &global_cells);
      CHECK(global_cells == 1000);

      // BELOW SEGMENT DEMONSTRATES MSTK ISSUE #112
      int max_cell = -1;
      for (int c = 0; c != local_cells; ++c) {
        max_cell = std::max(
          max_cell, parent_mesh->getMap(Amanzi::AmanziMesh::Entity_kind::CELL, false)->getGlobalElement(c));
      }
      int global_max_cell;
      Teuchos::reduceAll(*comm, Teuchos::REDUCE_MAX, 1, &max_cell, &global_max_cell);
      std::cout << "cell counts: total = " << global_cells << ", max = " << global_max_cell
                << std::endl;
      // REMOVE COMMENT ONCE MSTK IS UPDATED
      CHECK(global_max_cell == 999);
      // END SEGMENT

      // Which constructor?
      //auto newmesh = meshfactory.create(parent_mesh, setnames, Amanzi::AmanziMesh::Entity_kind::FACE, flatten);
      //if (local_count == 0) {
      //  CHECK(newmesh == Teuchos::null);
      //} else {
      //  CHECK(newmesh != Teuchos::null);
      //}
    }
  }
}
