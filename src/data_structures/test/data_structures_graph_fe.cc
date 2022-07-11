/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon (coonet@ornl.gov)
*/

//! <MISSING_ONELINE_DOCSTRING>

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "UnitTest++.h"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"

#include "AmanziComm.hh"
#include "AmanziMap.hh"
#include "VerboseObject.hh"
#include "Mesh.hh"
#include "MeshFactory.hh"
#include "GraphFE.hh"

/* *****************************************************************
 * this test is a null test -- all entries are local
 ***************************************************************** */
TEST(FE_GRAPH_NEAREST_NEIGHBOR_TPFA)
{
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;

  auto comm = getDefaultComm();
  int getRank = comm->getRank();

  if (getRank == 0)
    std::cout << "Test: FD like graph, null off-proc assembly" << std::endl;

  // read parameter list
  std::string xmlFileName = "test/operator_convergence.xml";
  Teuchos::ParameterXMLFileReader xmlreader(xmlFileName);
  Teuchos::ParameterList plist = xmlreader.getParameters();

  Amanzi::VerboseObject::global_hide_line_prefix = true;

  // create a mesh
  Teuchos::ParameterList region_list =
    plist.get<Teuchos::ParameterList>("regions");
  Teuchos::RCP<GeometricModel> gm =
    Teuchos::rcp(new GeometricModel(2, region_list, *comm));

  Preference pref;
  pref.clear();
  pref.push_back(Framework::MSTK);

  MeshFactory meshfactory(comm, gm);
  meshfactory.set_preference(pref);
  Teuchos::RCP<Mesh> mesh = meshfactory.create(0.0, 0.0, 1.0, 1.0, 10, 10);
  //  Teuchos::RCP<const Mesh> mesh =
  //  meshfactory.create("test/median32x33.exo");

  // grab the maps
  int ncells = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  auto cell_map = mesh->cell_map(false);
  auto cell_map_ghosted = mesh->cell_map(true);

  // create the graphs, one to test local, the other to test global insertion
  GraphFE graph_local(cell_map, cell_map_ghosted, cell_map_ghosted, 5);
  GraphFE graph_global(cell_map, cell_map_ghosted, cell_map_ghosted, 5);

  Kokkos::View<Entity_ID*,Kokkos::HostSpace> faces;
  Kokkos::View<Entity_ID*,Kokkos::HostSpace> face_cells;
  std::vector<int> neighbor_cells;
  for (int c = 0; c != ncells; ++c) {
    neighbor_cells.resize(0);
    neighbor_cells.push_back(c);

    mesh->cell_get_faces(c, faces);
    for (int n = 0; n != faces.size(); ++n) {
      mesh->face_get_cells_host(
        faces[n], AmanziMesh::Parallel_type::ALL, face_cells);
      if (face_cells.size() > 1) {
        neighbor_cells.push_back(c == face_cells[0] ? face_cells[1] :
                                                      face_cells[0]);
      }
    }

    GO global_c = cell_map->getGlobalElement(c);
    std::vector<GO> global_neighbors(neighbor_cells.size());
    for (int n=0; n!=neighbor_cells.size(); ++n)
      global_neighbors[n] = cell_map_ghosted->getGlobalElement(neighbor_cells[n]);

    graph_local.insertLocalIndices(c, neighbor_cells.size(), &neighbor_cells[0]);
    graph_global.insertGlobalIndices(
      global_c, neighbor_cells.size(), &global_neighbors[0]);
  }

  graph_local.fillComplete(cell_map, cell_map);
  graph_global.fillComplete(cell_map, cell_map);
}


/* *****************************************************************
 * this test is a null test -- all entries are local
 ***************************************************************** */
TEST(FE_GRAPH_FACE_FACE)
{
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;

  auto comm = getDefaultComm();
  int getRank = comm->getRank();

  if (getRank == 0)
    std::cout << "Test: FE like graph, off-proc assembly" << std::endl;

  // read parameter list
  std::string xmlFileName = "test/operator_convergence.xml";
  Teuchos::ParameterXMLFileReader xmlreader(xmlFileName);
  Teuchos::ParameterList plist = xmlreader.getParameters();

  Amanzi::VerboseObject::global_hide_line_prefix = true;

  // create a mesh
  Teuchos::ParameterList region_list =
    plist.get<Teuchos::ParameterList>("regions");
  Teuchos::RCP<GeometricModel> gm =
    Teuchos::rcp(new GeometricModel(2, region_list, *comm));

  Preference pref;
  pref.clear();
  pref.push_back(Framework::MSTK);

  MeshFactory meshfactory(comm, gm);
  meshfactory.set_preference(pref);
  Teuchos::RCP<Mesh> mesh = meshfactory.create(0.0, 0.0, 1.0, 1.0, 10, 10);
  //  Teuchos::RCP<const Mesh> mesh =
  //  meshfactory.create("test/median32x33.exo");

  // grab the maps
  int ncells = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  auto face_map = mesh->face_map(false);
  auto face_map_ghosted = mesh->face_map(true);

  // create the graph
  GraphFE graph_local(face_map, face_map_ghosted, face_map_ghosted, 7);
  GraphFE graph_global(face_map, face_map_ghosted, face_map_ghosted, 7);

  Kokkos::View<Entity_ID*,Kokkos::HostSpace> faces;
  Kokkos::View<Entity_ID*,Kokkos::HostSpace> face_cells;
  for (int c = 0; c != ncells; ++c) {
    mesh->cell_get_faces(c, faces);

    std::vector<GO> global_faces(faces.size());
    for (int n=0; n!=faces.size(); ++n)
      global_faces[n] = face_map_ghosted->getGlobalElement(faces[n]);
    
    for (int n=0; n!=faces.size(); ++n) {
      graph_local.insertLocalIndices(faces[n], faces.size(), &faces[0]);
      AMANZI_ASSERT(global_faces[n] >= 0);
      graph_global.insertGlobalIndices(
        global_faces[n], global_faces.size(), &global_faces[0]);
    }
  }

  graph_local.fillComplete(face_map, face_map);
  graph_global.fillComplete(face_map, face_map);
}
