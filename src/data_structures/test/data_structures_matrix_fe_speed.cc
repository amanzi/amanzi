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
#include "GraphFE.hh"
#include "MatrixFE.hh"

/* *****************************************************************
 * this test is a null test -- all entries are local
* **************************************************************** */

TEST(FE_MATRIX_NEAREST_NEIGHBOR_TPFA_Epetra_FECrs)
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
  Teuchos::RCP<Teuchos::Time> timer = Teuchos::TimeMonitor::getNewTimer("Cell FECrsMatrix");

  Amanzi::VerboseObject::global_hide_line_prefix = true;

  // create a mesh
  Teuchos::ParameterList region_list = plist.get<Teuchos::ParameterList>("regions");
  GeometricModelPtr gm = new GeometricModel(2, region_list, *comm);

  Preference pref;
  pref.clear();
  pref.push_back(Framework::MSTK);

  MeshFactory meshfactory(comm, gm);
  meshfactory.set_preference(pref);
  Teuchos::RCP<Mesh> mesh = meshfactory.create(0.0, 0.0, 1.0, 1.0, 100, 1000);
  //Teuchos::RCP<const Mesh> mesh = meshfactory.create("test/median32x33.exo");

  // grab the maps
  int ncells = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_kind::OWNED);
  Teuchos::RCP<Epetra_Map> cell_map = Teuchos::rcp(new Epetra_Map(mesh->cell_map(false)));
  Teuchos::RCP<Epetra_Map> cell_map_ghosted = Teuchos::rcp(new Epetra_Map(mesh->cell_map(true)));

  // create the graph
  int ierr(0);
  Teuchos::RCP<GraphFE> graph =
    Teuchos::rcp(new GraphFE(cell_map, cell_map_ghosted, cell_map_ghosted, 5));

  Entity_ID_View faces;
  Entity_ID_View face_cells;
  std::vector<int> neighbor_cells;
  for (int c = 0; c != ncells; ++c) {
    neighbor_cells.resize(0);
    neighbor_cells.push_back(c);

    mesh->cell_get_faces(c, &faces);
    for (int n = 0; n != faces.size(); ++n) {
      mesh->face_get_cells(faces[n], AmanziMesh::Parallel_kind::ALL, &face_cells);
      if (face_cells.size() > 1) {
        neighbor_cells.push_back(c == face_cells[0] ? face_cells[1] : face_cells[0]);
      }
    }

    ierr |= graph->InsertMyIndices(c, neighbor_cells.size(), &neighbor_cells[0]);
    CHECK(!ierr);
  }

  ierr |= graph->FillComplete(cell_map, cell_map);
  CHECK(!ierr);

  // and the control matrix
  Epetra_FECrsMatrix control(Copy, graph->Graph());
  {
    Teuchos::TimeMonitor monitor(*timer);
    for (int c = 0; c != ncells; ++c) {
      neighbor_cells.resize(0);
      neighbor_cells.push_back(c);

      mesh->cell_get_faces(c, &faces);
      for (int n = 0; n != faces.size(); ++n) {
        mesh->face_get_cells(faces[n], AmanziMesh::Parallel_kind::ALL, &face_cells);
        if (face_cells.size() > 1) {
          neighbor_cells.push_back(c == face_cells[0] ? face_cells[1] : face_cells[0]);
        }
      }

      std::vector<double> vals(neighbor_cells.size(), 1);
      std::vector<int> neighbor_cell_gids(neighbor_cells.size());
      for (int n = 0; n != neighbor_cells.size(); ++n) {
        neighbor_cell_gids[n] = cell_map_ghosted->GID(neighbor_cells[n]);
        AMANZI_ASSERT(neighbor_cell_gids[n] >= 0);
      }
      AMANZI_ASSERT(cell_map_ghosted->GID(c) >= 0);
      ierr |= control.SumIntoGlobalValues(
        cell_map_ghosted->GID(c), neighbor_cells.size(), &vals[0], &neighbor_cell_gids[0]);
      AMANZI_ASSERT(!ierr);
    }

    ierr |= control.GlobalAssemble();
  }

  CHECK(!ierr);
}

TEST(FE_MATRIX_NEAREST_NEIGHBOR_TPFA_Epetra_FECrs_Nonlocal)
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
  Teuchos::RCP<Teuchos::Time> timer = Teuchos::TimeMonitor::getNewTimer("Cell FECrsMatrix offproc");

  Amanzi::VerboseObject::global_hide_line_prefix = true;

  // create a mesh
  Teuchos::ParameterList region_list = plist.get<Teuchos::ParameterList>("regions");
  GeometricModelPtr gm = new GeometricModel(2, region_list, *comm);

  Preference pref;
  pref.clear();
  pref.push_back(Framework::MSTK);

  MeshFactory meshfactory(comm, gm);
  meshfactory.set_preference(pref);
  Teuchos::RCP<Mesh> mesh = meshfactory.create(0.0, 0.0, 1.0, 1.0, 100, 1000);
  //  Teuchos::RCP<const Mesh> mesh = meshfactory.create("test/median32x33.exo");

  // grab the maps
  int ncells = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_kind::OWNED);
  Teuchos::RCP<Epetra_Map> cell_map = Teuchos::rcp(new Epetra_Map(mesh->cell_map(false)));
  Teuchos::RCP<Epetra_Map> cell_map_ghosted = Teuchos::rcp(new Epetra_Map(mesh->cell_map(true)));

  // create the graph
  int ierr(0);
  Teuchos::RCP<Epetra_FECrsGraph> graph =
    Teuchos::rcp(new Epetra_FECrsGraph(Copy, *cell_map, *cell_map_ghosted, 5, false, true));

  Entity_ID_View faces;
  Entity_ID_View face_cells;
  std::vector<int> neighbor_cells;
  for (int c = 0; c != ncells; ++c) {
    neighbor_cells.resize(0);
    neighbor_cells.push_back(c);

    mesh->cell_get_faces(c, &faces);
    for (int n = 0; n != faces.size(); ++n) {
      mesh->face_get_cells(faces[n], AmanziMesh::Parallel_kind::ALL, &face_cells);
      if (face_cells.size() > 1) {
        neighbor_cells.push_back(c == face_cells[0] ? face_cells[1] : face_cells[0]);
      }
    }

    std::vector<int> neighbor_cell_gids(neighbor_cells.size());
    for (int n = 0; n != neighbor_cells.size(); ++n)
      neighbor_cell_gids[n] = cell_map_ghosted->GID(neighbor_cells[n]);

    ierr |= graph->InsertGlobalIndices(
      cell_map_ghosted->GID(c), neighbor_cells.size(), &neighbor_cell_gids[0]);
    CHECK(!ierr);
  }

  ierr |= graph->FillComplete(*cell_map, *cell_map);
  CHECK(!ierr);

  // and the control matrix
  Epetra_FECrsMatrix control(Copy, *graph);

  {
    Teuchos::TimeMonitor monitor(*timer);

    for (int c = 0; c != ncells; ++c) {
      neighbor_cells.resize(0);
      neighbor_cells.push_back(c);

      mesh->cell_get_faces(c, &faces);
      for (int n = 0; n != faces.size(); ++n) {
        mesh->face_get_cells(faces[n], AmanziMesh::Parallel_kind::ALL, &face_cells);
        if (face_cells.size() > 1) {
          neighbor_cells.push_back(c == face_cells[0] ? face_cells[1] : face_cells[0]);
        }
      }

      std::vector<int> neighbor_cell_gids(neighbor_cells.size());
      for (int n = 0; n != neighbor_cells.size(); ++n) {
        neighbor_cell_gids[n] = cell_map_ghosted->GID(neighbor_cells[n]);
        AMANZI_ASSERT(neighbor_cell_gids[n] >= 0);
      }

      std::vector<double> vals(neighbor_cells.size(), 1);
      AMANZI_ASSERT(cell_map_ghosted->GID(c) >= 0);
      ierr |= control.SumIntoGlobalValues(
        cell_map_ghosted->GID(c), neighbor_cells.size(), &vals[0], &neighbor_cell_gids[0]);
      AMANZI_ASSERT(!ierr);
    }

    ierr |= control.GlobalAssemble();
  }
  CHECK(!ierr);
}

TEST(FE_MATRIX_NEAREST_NEIGHBOR_TPFA_MatrixFE)
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
  Teuchos::RCP<Teuchos::Time> timer = Teuchos::TimeMonitor::getNewTimer("Cell MatrixFE");

  Amanzi::VerboseObject::global_hide_line_prefix = true;

  // create a mesh
  Teuchos::ParameterList region_list = plist.get<Teuchos::ParameterList>("regions");
  GeometricModelPtr gm = new GeometricModel(2, region_list, *comm);

  Preference pref;
  pref.clear();
  pref.push_back(Framework::MSTK);

  MeshFactory meshfactory(comm, gm);
  meshfactory.set_preference(pref);
  Teuchos::RCP<Mesh> mesh = meshfactory.create(0.0, 0.0, 1.0, 1.0, 100, 1000);
  //  Teuchos::RCP<const Mesh> mesh = meshfactory.create("test/median32x33.exo");

  // grab the maps
  int ncells = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_kind::OWNED);
  Teuchos::RCP<Epetra_Map> cell_map = Teuchos::rcp(new Epetra_Map(mesh->cell_map(false)));
  Teuchos::RCP<Epetra_Map> cell_map_ghosted = Teuchos::rcp(new Epetra_Map(mesh->cell_map(true)));

  // create the graph
  int ierr(0);
  Teuchos::RCP<GraphFE> graph =
    Teuchos::rcp(new GraphFE(cell_map, cell_map_ghosted, cell_map_ghosted, 5));

  Entity_ID_View faces;
  Entity_ID_View face_cells;
  std::vector<int> neighbor_cells;
  for (int c = 0; c != ncells; ++c) {
    neighbor_cells.resize(0);
    neighbor_cells.push_back(c);

    mesh->cell_get_faces(c, &faces);
    for (int n = 0; n != faces.size(); ++n) {
      mesh->face_get_cells(faces[n], AmanziMesh::Parallel_kind::ALL, &face_cells);
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
  {
    Teuchos::TimeMonitor monitor(*timer);

    for (int c = 0; c != ncells; ++c) {
      neighbor_cells.resize(0);
      neighbor_cells.push_back(c);

      mesh->cell_get_faces(c, &faces);
      for (int n = 0; n != faces.size(); ++n) {
        mesh->face_get_cells(faces[n], AmanziMesh::Parallel_kind::ALL, &face_cells);
        if (face_cells.size() > 1) {
          neighbor_cells.push_back(c == face_cells[0] ? face_cells[1] : face_cells[0]);
        }
      }

      std::vector<double> vals(neighbor_cells.size(), 1);
      ierr |= matrix.SumIntoMyValues(c, neighbor_cells.size(), &vals[0], &neighbor_cells[0]);
      AMANZI_ASSERT(!ierr);
    }

    ierr |= matrix.FillComplete();
  }
  CHECK(!ierr);
}


TEST(FE_MATRIX_FACE_FACE_Epetra_FECrsMatrix2)
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
  Teuchos::RCP<Teuchos::Time> timer = Teuchos::TimeMonitor::getNewTimer("Block FECrsMatrix");

  Amanzi::VerboseObject::global_hide_line_prefix = true;

  // create a mesh
  Teuchos::ParameterList region_list = plist.get<Teuchos::ParameterList>("regions");
  GeometricModelPtr gm = new GeometricModel(2, region_list, *comm);

  Preference pref;
  pref.clear();
  pref.push_back(Framework::MSTK);

  MeshFactory meshfactory(comm, gm);
  meshfactory.set_preference(pref);
  Teuchos::RCP<Mesh> mesh = meshfactory.create(0.0, 0.0, 1.0, 1.0, 100, 1000);
  //  Teuchos::RCP<const Mesh> mesh = meshfactory.create("test/median32x33.exo");

  // grab the maps
  int ncells = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_kind::OWNED);
  Teuchos::RCP<Epetra_Map> face_map = Teuchos::rcp(new Epetra_Map(mesh->face_map(false)));
  Teuchos::RCP<Epetra_Map> face_map_ghosted = Teuchos::rcp(new Epetra_Map(mesh->face_map(true)));

  // create the graph
  int ierr(0);
  Teuchos::RCP<GraphFE> graph =
    Teuchos::rcp(new GraphFE(face_map, face_map_ghosted, face_map_ghosted, 5));

  Entity_ID_View faces;
  Entity_ID_View face_cells;
  for (int c = 0; c != ncells; ++c) {
    mesh->cell_get_faces(c, &faces);

    for (int n = 0; n != faces.size(); ++n) {
      ierr |= graph->InsertMyIndices(faces[n], faces.size(), &faces[0]);
      CHECK(!ierr);
    }
  }

  ierr |= graph->FillComplete(face_map, face_map);
  CHECK(!ierr);

  // and the control matrix
  Epetra_FECrsMatrix control(Copy, graph->Graph());

  {
    Teuchos::TimeMonitor monitor(*timer);
    for (int c = 0; c != ncells; ++c) {
      mesh->cell_get_faces(c, &faces);
      Epetra_IntSerialDenseVector face_gids(faces.size());
      for (int n = 0; n != faces.size(); ++n) { face_gids[n] = face_map_ghosted->GID(faces[n]); }

      Epetra_SerialDenseMatrix vals(faces.size(), faces.size());
      ierr |= control.SumIntoGlobalValues(face_gids, vals);
    }
    ierr |= control.GlobalAssemble();
  }
  CHECK(!ierr);
}

TEST(FE_MATRIX_FACE_FACE_Epetra_FECrsMatrix_offproc2)
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
  Teuchos::RCP<Teuchos::Time> timer =
    Teuchos::TimeMonitor::getNewTimer("Block FECrsMatrix offproc");

  Amanzi::VerboseObject::global_hide_line_prefix = true;

  // create a mesh
  Teuchos::ParameterList region_list = plist.get<Teuchos::ParameterList>("regions");
  GeometricModelPtr gm = new GeometricModel(2, region_list, *comm);

  Preference pref;
  pref.clear();
  pref.push_back(Framework::MSTK);

  MeshFactory meshfactory(comm, gm);
  meshfactory.set_preference(pref);
  Teuchos::RCP<Mesh> mesh = meshfactory.create(0.0, 0.0, 1.0, 1.0, 100, 1000);
  //  Teuchos::RCP<const Mesh> mesh = meshfactory.create("test/median32x33.exo");

  // grab the maps
  int ncells = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_kind::OWNED);
  Teuchos::RCP<Epetra_Map> face_map = Teuchos::rcp(new Epetra_Map(mesh->face_map(false)));
  Teuchos::RCP<Epetra_Map> face_map_ghosted = Teuchos::rcp(new Epetra_Map(mesh->face_map(true)));

  // create the graph
  int ierr(0);
  Teuchos::RCP<Epetra_FECrsGraph> graph =
    Teuchos::rcp(new Epetra_FECrsGraph(Copy, *face_map, *face_map_ghosted, 8, false, true));

  Entity_ID_View faces;
  for (int c = 0; c != ncells; ++c) {
    mesh->cell_get_faces(c, &faces);

    std::vector<int> face_gids(faces.size());
    for (int n = 0; n != faces.size(); ++n) face_gids[n] = face_map_ghosted->GID(faces[n]);

    ierr |=
      graph->InsertGlobalIndices(face_gids.size(), &face_gids[0], face_gids.size(), &face_gids[0]);
    AMANZI_ASSERT(!ierr);
  }

  ierr |= graph->FillComplete(*face_map, *face_map);
  CHECK(!ierr);

  // and the control matrix
  Epetra_FECrsMatrix control(Copy, *graph);

  {
    Teuchos::TimeMonitor monitor(*timer);
    for (int c = 0; c != ncells; ++c) {
      mesh->cell_get_faces(c, &faces);
      Epetra_IntSerialDenseVector face_gids(faces.size());
      for (int n = 0; n != faces.size(); ++n) { face_gids[n] = face_map_ghosted->GID(faces[n]); }

      Epetra_SerialDenseMatrix vals(faces.size(), faces.size());
      ierr |= control.SumIntoGlobalValues(face_gids, vals);
    }

    ierr |= control.GlobalAssemble();
  }
  CHECK(!ierr);
}

TEST(FE_MATRIX_FACE_FACE_MatrixFE2)
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
  Teuchos::RCP<Teuchos::Time> timer = Teuchos::TimeMonitor::getNewTimer("Block MatrixFE");

  Amanzi::VerboseObject::global_hide_line_prefix = true;

  // create a mesh
  Teuchos::ParameterList region_list = plist.get<Teuchos::ParameterList>("regions");
  GeometricModelPtr gm = new GeometricModel(2, region_list, *comm);

  Preference pref;
  pref.clear();
  pref.push_back(Framework::MSTK);

  MeshFactory meshfactory(comm, gm);
  meshfactory.set_preference(pref);
  Teuchos::RCP<Mesh> mesh = meshfactory.create(0.0, 0.0, 1.0, 1.0, 100, 1000);
  //  Teuchos::RCP<const Mesh> mesh = meshfactory.create("test/median32x33.exo");

  // grab the maps
  int ncells = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_kind::OWNED);
  Teuchos::RCP<Epetra_Map> face_map = Teuchos::rcp(new Epetra_Map(mesh->face_map(false)));
  Teuchos::RCP<Epetra_Map> face_map_ghosted = Teuchos::rcp(new Epetra_Map(mesh->face_map(true)));

  // create the graph
  int ierr(0);
  Teuchos::RCP<GraphFE> graph =
    Teuchos::rcp(new GraphFE(face_map, face_map_ghosted, face_map_ghosted, 5));

  Entity_ID_View faces;
  Entity_ID_View face_cells;
  for (int c = 0; c != ncells; ++c) {
    mesh->cell_get_faces(c, &faces);

    for (int n = 0; n != faces.size(); ++n) {
      ierr |= graph->InsertMyIndices(faces[n], faces.size(), &faces[0]);
      CHECK(!ierr);
    }
  }

  ierr |= graph->FillComplete(face_map, face_map);
  CHECK(!ierr);

  // create the test matrix
  MatrixFE matrix(graph);

  {
    Teuchos::TimeMonitor monitor(*timer);

    for (int c = 0; c != ncells; ++c) {
      mesh->cell_get_faces(c, &faces);

      Epetra_SerialDenseMatrix vals(faces.size(), faces.size());
      ierr |= matrix.SumIntoMyValues(&faces[0], vals);
    }

    ierr |= matrix.FillComplete();
  }
  CHECK(!ierr);
}
