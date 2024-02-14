/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

//! <MISSING_ONELINE_DOCSTRING>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>
#include <numeric>

#include "UnitTest++.h"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"

#include "MeshFactory.hh"

#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_SerialDenseVector.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"

#include "AmanziComm.hh"
#include "AmanziMatrix.hh"
#include "VerboseObject.hh"
#include "Mesh.hh"
#include "MeshFactory.hh"
#include "GraphFE.hh"
#include "MatrixFE.hh"

/* *****************************************************************
 * this test is a null test -- all entries are local
 * **************************************************************** */
TEST(FE_MATRIX_NEAREST_NEIGHBOR_TPFA)
{
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;

  auto comm = getDefaultComm();
  int getRank = comm->getRank();

  if (getRank == 0) std::cout << "Test: FD like matrix, null off-proc assembly" << std::endl;

  // read parameter list
  std::string xmlFileName = "test/operator_convergence.xml";
  Teuchos::ParameterXMLFileReader xmlreader(xmlFileName);
  Teuchos::ParameterList plist = xmlreader.getParameters();

  Amanzi::VerboseObject::global_hide_line_prefix = true;

  // create a mesh
  Teuchos::ParameterList region_list = plist.get<Teuchos::ParameterList>("regions");
  Teuchos::RCP<GeometricModel> gm = Teuchos::rcp(new GeometricModel(2, region_list, *comm));

  Preference pref;
  pref.clear();
  pref.push_back(Framework::MSTK);

  MeshFactory meshfactory(comm, gm);
  meshfactory.set_preference(pref);
  Teuchos::RCP<Mesh> mesh = meshfactory.create(0.0, 0.0, 1.0, 1.0, 10, 10);
  //  Teuchos::RCP<const Mesh> mesh =
  //  meshfactory.create("test/median32x33.exo");

  // grab the maps
  int ncells = mesh->getNumEntities(AmanziMesh::CELL, AmanziMesh::Parallel_kind::OWNED);
  auto cell_map = mesh->getMap(Entity_kind::CELL, false);
  auto cell_map_ghosted = mesh->getMap(Entity_kind::CELL, true);

  // create the graph
  auto graph = Teuchos::rcp(new GraphFE(cell_map, cell_map_ghosted, cell_map_ghosted, 5));

  std::vector<int> neighbor_cells;
  for (int c = 0; c != ncells; ++c) {
    neighbor_cells.resize(0);
    neighbor_cells.push_back(c);

    auto faces = mesh->getCellFaces(c);
    for (int n = 0; n != faces.size(); ++n) {
      auto face_cells = mesh->getFaceCells(faces[n]);
      if (face_cells.size() > 1) {
        neighbor_cells.push_back(c == face_cells[0] ? face_cells[1] : face_cells[0]);
      }
    }

    graph->insertLocalIndices(c, neighbor_cells.size(), &neighbor_cells[0]);
  }

  graph->fillComplete(cell_map, cell_map);

  // create the test matrix
  MatrixFE matrix(graph);

  // and the control matrix
  Matrix_type control(graph->getGraph());

  for (int c = 0; c != ncells; ++c) {
    neighbor_cells.resize(0);
    neighbor_cells.push_back(c);

    auto faces = mesh->getCellFaces(c);
    for (int n = 0; n != faces.size(); ++n) {
      auto face_cells = mesh->getFaceCells(faces[n]);
      if (face_cells.size() > 1) {
        neighbor_cells.push_back(c == face_cells[0] ? face_cells[1] : face_cells[0]);
      }
    }

    Teuchos::SerialDenseMatrix<int, double> vals(neighbor_cells.size(), 1);
    vals.random();

    matrix.sumIntoLocalValues(c, neighbor_cells.size(), &vals(0, 0), neighbor_cells.data());

    std::vector<GO> neighbor_cell_gids(neighbor_cells.size());
    for (int n = 0; n != neighbor_cells.size(); ++n) {
      neighbor_cell_gids[n] = cell_map_ghosted->getGlobalElement(neighbor_cells[n]);
      AMANZI_ASSERT(neighbor_cell_gids[n] >= 0);
    }
    AMANZI_ASSERT(cell_map_ghosted->getGlobalElement(c) >= 0);
    control.sumIntoGlobalValues(cell_map_ghosted->getGlobalElement(c),
                                neighbor_cells.size(),
                                &vals(0, 0),
                                neighbor_cell_gids.data());
  }

  matrix.fillComplete();

  control.fillComplete(graph->getDomainMap(), graph->getRangeMap());

  auto os = Teuchos::getFancyOStream(Teuchos::rcpFromRef(std::cout));
  matrix.getMatrix()->describe(*os, Teuchos::VERB_EXTREME);
  control.describe(*os, Teuchos::VERB_EXTREME);


  // check matrix equality
  for (int c = 0; c != ncells; ++c) {
    MatrixFE::cValues_host_view_type mat_vals, ctrl_vals;
    MatrixFE::cLocalIndices_host_view_type mat_inds, ctrl_inds;
    int ctrl_nentries, mat_nentries;

    matrix.getLocalRowView(c, mat_inds, mat_vals);
    mat_nentries = mat_vals.size();

    control.getLocalRowView(c, ctrl_inds, ctrl_vals);
    ctrl_nentries = ctrl_vals.size();

    CHECK_EQUAL(ctrl_nentries, mat_nentries);
    for (int i = 0; i != std::min(ctrl_nentries, mat_nentries); ++i) {
      CHECK_EQUAL(ctrl_inds[i], mat_inds[i]);
      CHECK_CLOSE(ctrl_vals[i], mat_vals[i], 1.e-10);
    }
  }
}


/* *****************************************************************
 * this test is a real test with FE-like assembly of face-face system
 * **************************************************************** */
TEST(FE_MATRIX_FACE_FACE)
{
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;

  auto comm = getDefaultComm();
  int getRank = comm->getRank();

  if (getRank == 0) std::cout << "Test: FE like matrix, off-proc assembly" << std::endl;

  // read parameter list
  std::string xmlFileName = "test/operator_convergence.xml";
  Teuchos::ParameterXMLFileReader xmlreader(xmlFileName);
  Teuchos::ParameterList plist = xmlreader.getParameters();

  Amanzi::VerboseObject::global_hide_line_prefix = true;

  // create a mesh
  Teuchos::ParameterList region_list = plist.get<Teuchos::ParameterList>("regions");
  Teuchos::RCP<GeometricModel> gm = Teuchos::rcp(new GeometricModel(2, region_list, *comm));

  MeshFactory meshfactory(comm, gm);
  meshfactory.set_preference({ Framework::MSTK });
  Teuchos::RCP<Mesh> mesh = meshfactory.create(0.0, 0.0, 1.0, 1.0, 10, 10);

  // grab the maps
  int ncells = mesh->getNumEntities(AmanziMesh::CELL, AmanziMesh::Parallel_kind::OWNED);
  auto face_map = mesh->getMap(Entity_kind::FACE, false);
  auto face_map_ghosted = mesh->getMap(Entity_kind::FACE, true);

  // create the graph
  Teuchos::RCP<GraphFE> graph =
    Teuchos::rcp(new GraphFE(face_map, face_map_ghosted, face_map_ghosted, 7));

  for (int c = 0; c != ncells; ++c) {
    auto faces = mesh->getCellFaces(c);
    for (int n = 0; n != faces.size(); ++n) {
      graph->insertLocalIndices(faces[n], faces.size(), &faces(0));
    }
  }

  graph->fillComplete(face_map, face_map);

  // create the test matrix
  MatrixFE matrix(graph);

  // and the control matrix
  Matrix_type control(graph->getGraph());

  for (int c = 0; c != ncells; ++c) {
    auto faces = mesh->getCellFaces(c);

    Teuchos::SerialDenseVector<int, GO> face_gids(faces.size());
    for (int n = 0; n != faces.size(); ++n) {
      face_gids[n] = face_map_ghosted->getGlobalElement(faces[n]);
      AMANZI_ASSERT(face_gids[n] > -1);
    }

    Teuchos::SerialDenseMatrix<int, double> vals(faces.size(), faces.size());
    vals.random();
    matrix.sumIntoLocalValuesTransposed(faces.data(), faces.data(), vals);

    for (int i = 0; i != faces.size(); ++i) {
      control.sumIntoGlobalValues(face_gids[i], faces.size(), &vals(0, i), &face_gids[0]);
    }
  }

  matrix.fillComplete();
  control.fillComplete(graph->getDomainMap(), graph->getRangeMap());

  // check matrix equality
  int nfaces = mesh->getNumEntities(AmanziMesh::FACE, AmanziMesh::Parallel_kind::OWNED);
  for (int f = 0; f != nfaces; ++f) {
    MatrixFE::cValues_host_view_type mat_vals, ctrl_vals;
    MatrixFE::cLocalIndices_host_view_type mat_inds, ctrl_inds;
    int ctrl_nentries, mat_nentries;

    matrix.getLocalRowView(f, mat_inds, mat_vals);
    mat_nentries = mat_vals.size();

    control.getLocalRowView(f, ctrl_inds, ctrl_vals);
    ctrl_nentries = ctrl_vals.size();

    CHECK_EQUAL(ctrl_nentries, mat_nentries);
    for (int i = 0; i != std::min(ctrl_nentries, mat_nentries); ++i) {
      CHECK_EQUAL(ctrl_inds[i], mat_inds[i]);
      CHECK_CLOSE(ctrl_vals[i], mat_vals[i], 1.e-10);
    }
  }
}
