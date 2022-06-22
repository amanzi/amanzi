/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:

*/

//!

#include <UnitTest++.h>

#include <fstream>

#include "../Mesh_MSTK.hh"

#include "MeshAudit.hh"

#include "AmanziMap.hh"
#include "AmanziComm.hh"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

// NOTE: WHEN THESE TESTS ARE RUN IN PARALLEL, THE CELLS ABOVE THE
// DEFORMED CELLS WILL NOT SHIFT DOWN CORRECTLY BECAUSE THE COLUMNS
// ARE NOT ALL ON ONE PROCESSOR. IN A REAL SUBSURFACE DEFORMATION
// SIMULATION, THE PARALLEL PARTITIONING WILL BE IN THE LATERAL DIRECTION
// ONLY AND EACH COLUMN WILL BE ON A SINGLE PROCESSOR, SO THE DEFORMATION
// OF CELLS WILL BE PROPERLY COMMUNICATED TO THE LAYERS ABOVE IT

TEST(MSTK_DEFORM_VOLS_2D)
{
  auto comm = Amanzi::getDefaultComm();

  // Define a box region to capture bottom boundary
  Teuchos::ParameterList param_list;
  Teuchos::ParameterList& regions_list = param_list.sublist("regions");
  Teuchos::ParameterList& botreg_list = regions_list.sublist("Bottom Region");
  Teuchos::ParameterList& botreg_def = botreg_list.sublist("region: box");
  Teuchos::Array<double> lo_coord = Teuchos::tuple(-5.1, -0.01);
  Teuchos::Array<double> hi_coord = Teuchos::tuple(5.1, 0.01);
  botreg_def.set<Teuchos::Array<double>>("low coordinate", lo_coord);
  botreg_def.set<Teuchos::Array<double>>("high coordinate", hi_coord);

  //  Teuchos::writeParameterListToXmlOStream(param_list,std::cout);

  if (comm->getSize() > 1) return;

  // Create a geometric model from region spec
  Teuchos::RCP<Amanzi::AmanziGeometry::GeometricModel> gm = Teuchos::rcp(
    new Amanzi::AmanziGeometry::GeometricModel(3, regions_list, *comm));

  // Generate a mesh consisting of 3x3 elements
  auto plist = Teuchos::rcp(new Teuchos::ParameterList());
  plist->sublist("unstructured")
    .sublist("expert")
    .set<std::string>("partitioner", "zoltan_rcb");
  bool request_faces = true;
  bool request_edges = true;
  Teuchos::RCP<Amanzi::AmanziMesh::Mesh> mesh(
    new Amanzi::AmanziMesh::Mesh_MSTK(-5.0,
                                      0.0,
                                      5.0,
                                      10.0,
                                      10,
                                      10,
                                      comm,
                                      gm,
                                      plist,
                                      request_faces,
                                      request_edges));

  CHECK_EQUAL(mesh->build_columns(), 1);

  int nc = mesh->num_entities(Amanzi::AmanziMesh::CELL,
                              Amanzi::AmanziMesh::Parallel_type::ALL);

  // Request target volume of 50% for some cells at the bottom of the center
  // column The others are unconstrained except for a barrier of minimum volume

  std::vector<double> orig_volumes, target_volumes, target_weights, min_volumes;
  orig_volumes.reserve(nc);
  target_volumes.reserve(nc);
  target_weights.reserve(nc);
  min_volumes.reserve(nc);

  for (int i = 0; i < nc; i++) {
    orig_volumes[i] = mesh->cell_volume_host(i);
    // target_volumes[i] = orig_volumes[i];
    target_volumes[i] = 0.0;
    min_volumes[i] = 0.25 * orig_volumes[i];

    Amanzi::AmanziGeometry::Point ccen = mesh->cell_centroid_host(i);
    if (fabs(ccen[0]) < 3.0)
      if (ccen[1] > 3.1 && ccen[1] < 4.1)
        target_volumes[i] = 0.90 * orig_volumes[i];
  }

  Amanzi::AmanziMesh::Entity_ID_List fixed_nodes;
  mesh->get_set_entities("Bottom Region",
                         Amanzi::AmanziMesh::NODE,
                         Amanzi::AmanziMesh::Parallel_type::ALL,
                         fixed_nodes);

  bool move_vertical = true;
  int status =
    mesh->deform(target_volumes, min_volumes, fixed_nodes, move_vertical);
  CHECK(status);

  for (int i = 0; i < nc; i++) {
    if (target_volumes[i] > 0.0 && target_volumes[i] < orig_volumes[i]) {
      double voldiff =
        (mesh->cell_volume_host(i) - target_volumes[i]) / target_volumes[i];

      // Check if volume difference is with 5% of target volume
      CHECK_CLOSE(0, voldiff, 5.0e-02);
    }

    // Check that we didn't fall below the minimum prescribed volume
    if (mesh->cell_volume_host(i) < min_volumes[i])
      std::cerr << "Cell volume = " << mesh->cell_volume_host(i)
                << " is less than min volume = " << min_volumes[i] << std::endl;
    //    CHECK(mesh->cell_volume(i,false) >= min_volumes[i]);
  }

  mesh->write_to_exodus_file("deformed2.exo");
}


TEST(MSTK_DEFORM_VOLS_3D)
{
  auto comm = Amanzi::getDefaultComm();

  // Define a box region to capture bottom boundary
  Teuchos::ParameterList param_list;
  Teuchos::ParameterList& regions_list = param_list.sublist("regions");
  Teuchos::ParameterList& botreg_list = regions_list.sublist("Bottom Region");
  Teuchos::ParameterList& botreg_def = botreg_list.sublist("region: box");
  Teuchos::Array<double> lo_coord = Teuchos::tuple(-0.1, -0.1, -0.01);
  Teuchos::Array<double> hi_coord = Teuchos::tuple(10.1, 10.1, 0.01);
  botreg_def.set<Teuchos::Array<double>>("low coordinate", lo_coord);
  botreg_def.set<Teuchos::Array<double>>("high coordinate", hi_coord);

  //  Teuchos::writeParameterListToXmlOStream(param_list,std::cout);

  // Create a geometric model from region spec
  auto gm = Teuchos::rcp(
    new Amanzi::AmanziGeometry::GeometricModel(3, regions_list, *comm));

  // Generate a mesh consisting of 10x10 elements
  auto plist = Teuchos::rcp(new Teuchos::ParameterList());
  plist->sublist("unstructured")
    .sublist("expert")
    .set<std::string>("partitioner", "zoltan_rcb");
  bool request_faces = true;
  bool request_edges = true;
  Teuchos::RCP<Amanzi::AmanziMesh::Mesh> mesh(
    new Amanzi::AmanziMesh::Mesh_MSTK(0.0,
                                      0.0,
                                      0.0,
                                      10.0,
                                      1.0,
                                      10.0,
                                      10,
                                      1,
                                      10,
                                      comm,
                                      gm,
                                      plist,
                                      request_faces,
                                      request_edges));

  CHECK_EQUAL(mesh->build_columns(), 1);

  int ncused = mesh->num_entities(Amanzi::AmanziMesh::CELL,
                                  Amanzi::AmanziMesh::Parallel_type::ALL);
  int ncowned = mesh->num_entities(Amanzi::AmanziMesh::CELL,
                                   Amanzi::AmanziMesh::Parallel_type::OWNED);

  // Request target volume of 50% for some cells at the bottom of the center
  // column The others are unconstrained except for a barrier of minimum volume

  std::vector<double> orig_volumes, target_volumes, target_weights, min_volumes;
  orig_volumes.reserve(ncused);
  target_volumes.reserve(ncused);
  target_weights.reserve(ncused);
  min_volumes.reserve(ncused);

  for (int i = 0; i < ncused; i++) {
    orig_volumes[i] = mesh->cell_volume_host(i);
    target_volumes[i] = orig_volumes[i];
    min_volumes[i] = 0.90 * orig_volumes[i];

    Amanzi::AmanziGeometry::Point ccen = mesh->cell_centroid_host(i);
    if (ccen[0] > 2.0 && ccen[0] < 8.0) { // row of cells along x axis
      if (ccen[2] > 3.1 && ccen[2] < 4.1) {
        target_volumes[i] = 0.85 * orig_volumes[i];
        min_volumes[i] = 0.770 * orig_volumes[i];
      } else if (ccen[2] >= 4.1) {
        min_volumes[i] = 0.7 * orig_volumes[i];
      }
    }
  }

  Amanzi::AmanziMesh::Entity_ID_List fixed_nodes;
  mesh->get_set_entities_and_vofs("Bottom Region",
                                  Amanzi::AmanziMesh::NODE,
                                  Amanzi::AmanziMesh::Parallel_type::ALL,
                                  fixed_nodes,
                                  NULL);

  bool move_vertical = true;
  int status =
    mesh->deform(target_volumes, min_volumes, fixed_nodes, move_vertical);
  CHECK(status);

  for (int i = 0; i < ncowned; i++) {
    if (target_volumes[i] > 0.0 && target_volumes[i] < orig_volumes[i]) {
      double voldiff =
        (mesh->cell_volume_host(i) - target_volumes[i]) / target_volumes[i];

      // Check if volume difference is with 5% of target volume
      CHECK_CLOSE(0, voldiff, 5e-02);
    }

    // Check that we didn't fall below the minimum prescribed volume
    // but because we asked for many volumes to stay exactly the same
    // as the original (min_volume = target_volume = original_volume)
    // give some margin of numerical error

    double eps = 1.0e-6 * orig_volumes[i];
    CHECK(mesh->cell_volume_host(i) + eps >= min_volumes[i]);
    if (!(mesh->cell_volume_host(i) + eps >= min_volumes[i])) {
      double diff = mesh->cell_volume_host(i) - min_volumes[i];
      std::cerr << "Cell Global ID "
                << mesh->getGlobalElement(i, Amanzi::AmanziMesh::CELL)
                << " Cell Local ID " << i << " Rank " << comm->getRank()
                << ": Min volume = " << min_volumes[i] << "    "
                << "Cell volume = " << mesh->cell_volume_host(i)
                << "  Diff = " << diff << std::endl;
    }
  }

  mesh->write_to_exodus_file("deformed3.exo");
}
