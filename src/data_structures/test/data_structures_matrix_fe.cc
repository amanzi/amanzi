/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  Data Structures

*/

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "UnitTest++.h"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"

#include "MeshFactory.hh"

#include "Teuchos_TimeMonitor.hpp"
#include "Epetra_IntSerialDenseVector.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_FECrsGraph.h"
#include "Epetra_FECrsMatrix.h"

#include "AmanziComm.hh"
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
  using namespace Amanzi::Operators;

  auto comm = getDefaultComm();
  int MyPID = comm->MyPID();

  if (MyPID == 0) std::cout << "Test: FD like matrix, null off-proc assembly" << std::endl;

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
  //  Teuchos::RCP<const Mesh> mesh = meshfactory.create("test/median32x33.exo");

  // grab the maps
  int ncells =
    mesh->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);
  Teuchos::RCP<Epetra_Map> cell_map =
    Teuchos::rcp(new Epetra_Map(mesh->getMap(AmanziMesh::Entity_kind::CELL, false)));
  Teuchos::RCP<Epetra_Map> cell_map_ghosted =
    Teuchos::rcp(new Epetra_Map(mesh->getMap(AmanziMesh::Entity_kind::CELL, true)));

  // create the graph
  int ierr(0);
  Teuchos::RCP<GraphFE> graph =
    Teuchos::rcp(new GraphFE(cell_map, cell_map_ghosted, cell_map_ghosted, 5));

  std::vector<int> neighbor_cells;
  for (int c = 0; c != ncells; ++c) {
    neighbor_cells.resize(0);
    neighbor_cells.push_back(c);

    auto faces = mesh->getCellFaces(c);
    for (int n = 0; n != faces.size(); ++n) {
      auto face_cells = mesh->getFaceCells(faces[n], AmanziMesh::Parallel_kind::ALL);
      if (face_cells.size() > 1) {
        neighbor_cells.push_back(c == face_cells[0] ? face_cells[1] : face_cells[0]);
      }
    }

    ierr |= graph->InsertMyIndices(c, neighbor_cells.size(), &neighbor_cells[0]);
    CHECK(!ierr);
  }

  ierr |= graph->FillComplete(cell_map, cell_map);
  CHECK(!ierr);

  // create the test matrix
  MatrixFE matrix(graph);

  // and the control matrix
  Epetra_FECrsMatrix control(Copy, graph->Graph());

  for (int c = 0; c != ncells; ++c) {
    neighbor_cells.resize(0);
    neighbor_cells.push_back(c);

    auto faces = mesh->getCellFaces(c);
    for (int n = 0; n != faces.size(); ++n) {
      auto face_cells = mesh->getFaceCells(faces[n], AmanziMesh::Parallel_kind::ALL);
      if (face_cells.size() > 1) {
        neighbor_cells.push_back(c == face_cells[0] ? face_cells[1] : face_cells[0]);
      }
    }

    Epetra_SerialDenseVector vals(neighbor_cells.size());
    vals.Random();

    ierr |= matrix.SumIntoMyValues(c, neighbor_cells.size(), &vals[0], &neighbor_cells[0]);
    AMANZI_ASSERT(!ierr);
    CHECK(!ierr);

    std::vector<int> neighbor_cell_gids(neighbor_cells.size());
    for (int n = 0; n != neighbor_cells.size(); ++n) {
      neighbor_cell_gids[n] = cell_map_ghosted->GID(neighbor_cells[n]);
      AMANZI_ASSERT(neighbor_cell_gids[n] >= 0);
    }
    AMANZI_ASSERT(cell_map_ghosted->GID(c) >= 0);
    ierr |= control.SumIntoGlobalValues(
      cell_map_ghosted->GID(c), neighbor_cells.size(), &vals[0], &neighbor_cell_gids[0]);
    AMANZI_ASSERT(!ierr);
    CHECK(!ierr);
  }

  ierr |= matrix.FillComplete();
  CHECK(!ierr);

  ierr |= control.GlobalAssemble();
  CHECK(!ierr);

  // check matrix equality
  for (int c = 0; c != ncells; ++c) {
    int nentries(0);
    std::vector<double> mat_vals(5);
    std::vector<double> ctrl_vals(5);
    std::vector<int> mat_inds(5);
    std::vector<int> ctrl_inds(5);

    ierr |= matrix.Matrix()->ExtractMyRowCopy(c, 5, nentries, &mat_vals[0], &mat_inds[0]);
    CHECK(!ierr);
    mat_vals.resize(nentries);
    mat_inds.resize(nentries);

    ierr |= control.ExtractMyRowCopy(c, 5, nentries, &ctrl_vals[0], &ctrl_inds[0]);
    CHECK(!ierr);
    ctrl_vals.resize(nentries);
    ctrl_inds.resize(nentries);

    CHECK(mat_vals == ctrl_vals);
    CHECK(mat_inds == ctrl_inds);
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
  using namespace Amanzi::Operators;

  auto comm = getDefaultComm();
  int MyPID = comm->MyPID();

  if (MyPID == 0) std::cout << "Test: FE like matrix, off-proc assembly" << std::endl;

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
  //  Teuchos::RCP<const Mesh> mesh = meshfactory.create("test/median32x33.exo");

  // grab the maps
  int ncells =
    mesh->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);
  Teuchos::RCP<Epetra_Map> face_map =
    Teuchos::rcp(new Epetra_Map(mesh->getMap(AmanziMesh::Entity_kind::FACE, false)));
  Teuchos::RCP<Epetra_Map> face_map_ghosted =
    Teuchos::rcp(new Epetra_Map(mesh->getMap(AmanziMesh::Entity_kind::FACE, true)));

  // create the graph
  int ierr(0);
  Teuchos::RCP<GraphFE> graph =
    Teuchos::rcp(new GraphFE(face_map, face_map_ghosted, face_map_ghosted, 5));

  for (int c = 0; c != ncells; ++c) {
    auto cfaces = mesh->getCellFaces(c);
    AmanziMesh::Entity_ID_View faces;
    faces.fromConst(cfaces);

    for (int n = 0; n != faces.size(); ++n) {
      ierr |= graph->InsertMyIndices(faces[n], faces.size(), &faces[0]);
      CHECK(!ierr);
    }
  }

  ierr |= graph->FillComplete(face_map, face_map);
  CHECK(!ierr);

  // create the test matrix
  MatrixFE matrix(graph);

  // and the control matrix
  Epetra_FECrsMatrix control(Copy, graph->Graph());

  for (int c = 0; c != ncells; ++c) {
    auto faces = mesh->getCellFaces(c);

    Epetra_IntSerialDenseVector face_gids(faces.size());
    for (int n = 0; n != faces.size(); ++n) {
      face_gids[n] = face_map_ghosted->GID(faces[n]);
      AMANZI_ASSERT(face_gids[n] > -1);
    }

    Epetra_SerialDenseMatrix vals(faces.size(), faces.size());
    vals.Random();

    ierr |= matrix.SumIntoMyValues(&faces[0], vals);
    CHECK(!ierr);

    ierr |= control.SumIntoGlobalValues(face_gids, vals);
    CHECK(!ierr);
  }

  ierr |= matrix.FillComplete();
  CHECK(!ierr);

  ierr |= control.GlobalAssemble();
  CHECK(!ierr);

  // check matrix equality
  int nfaces =
    mesh->getNumEntities(AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_kind::OWNED);
  for (int f = 0; f != nfaces; ++f) {
    int nentries(0);
    std::vector<double> mat_vals(7);
    std::vector<double> ctrl_vals(7);
    std::vector<int> mat_inds(7);
    std::vector<int> ctrl_inds(7);

    ierr |= matrix.Matrix()->ExtractMyRowCopy(f, 7, nentries, &mat_vals[0], &mat_inds[0]);
    CHECK(!ierr);
    mat_vals.resize(nentries);
    mat_inds.resize(nentries);

    ierr |= control.ExtractMyRowCopy(f, 7, nentries, &ctrl_vals[0], &ctrl_inds[0]);
    CHECK(!ierr);
    ctrl_vals.resize(nentries);
    ctrl_inds.resize(nentries);

    CHECK(mat_inds == ctrl_inds);
    CHECK(mat_vals == ctrl_vals);
    if (!(mat_vals == ctrl_vals)) {
      std::cout << "Bad mat: ";
      for (std::vector<double>::const_iterator it = mat_vals.begin(); it != mat_vals.end();
           ++it)
        std::cout << " " << *it;
      std::cout << std::endl << "   ctrl: ";
      for (std::vector<double>::const_iterator it = ctrl_vals.begin(); it != ctrl_vals.end();
           ++it)
        std::cout << " " << *it;
      std::cout << std::endl;
    }
  }
}
