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
  using namespace Amanzi::Operators;

  auto comm = getDefaultComm();
  int getRank = comm->getRank();

  if (getRank == 0)
    std::cout << "Test: FD like matrix, null off-proc assembly" << std::endl;

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

  // create the graph
  int ierr(0);
  auto graph = Teuchos::rcp(new GraphFE(cell_map, cell_map_ghosted, cell_map_ghosted, 5));
  
  Entity_ID_View faces;
  Entity_ID_View face_cells;
  std::vector<int> neighbor_cells;
  for (int c = 0; c != ncells; ++c) {
    neighbor_cells.resize(0);
    neighbor_cells.push_back(c);

    mesh->cell_get_faces(c, faces);
    for (int n = 0; n != faces.size(); ++n) {
      mesh->face_get_cells(
        faces[n], AmanziMesh::Parallel_type::ALL, face_cells);
      if (face_cells.size() > 1) {
        neighbor_cells.push_back(c == face_cells[0] ? face_cells[1] :
                                                      face_cells[0]);
      }
    }

    ierr |=
      graph->InsertMyIndices(c, neighbor_cells.size(), &neighbor_cells[0]);
    CHECK(!ierr);
  }

  ierr |= graph->FillComplete(cell_map, cell_map);
  CHECK(!ierr);

  // create the test matrix
  MatrixFE matrix(graph);

  // and the control matrix
  Matrix_type control(graph->Graph());

  for (int c = 0; c != ncells; ++c) {
    neighbor_cells.resize(0);
    neighbor_cells.push_back(c);

    mesh->cell_get_faces(c, faces);
    for (int n = 0; n != faces.size(); ++n) {
      mesh->face_get_cells(
        faces[n], AmanziMesh::Parallel_type::ALL, face_cells);
      if (face_cells.size() > 1) {
        neighbor_cells.push_back(c == face_cells[0] ? face_cells[1] :
                                                      face_cells[0]);
      }
    }

    Teuchos::SerialDenseMatrix<int,double> vals(neighbor_cells.size(), 1);
    vals.random();
    
    ierr |= matrix.SumIntoMyValues(c, neighbor_cells.size(), &vals(0,0), neighbor_cells.data());
    AMANZI_ASSERT(!ierr);
    CHECK(!ierr);

    std::vector<int> neighbor_cell_gids(neighbor_cells.size());
    for (int n=0; n!=neighbor_cells.size(); ++n) {
      neighbor_cell_gids[n] = cell_map_ghosted->getGlobalElement(neighbor_cells[n]);
      AMANZI_ASSERT(neighbor_cell_gids[n] >= 0);
    }
    AMANZI_ASSERT(cell_map_ghosted->getGlobalElement(c) >= 0);
    control.sumIntoGlobalValues(cell_map_ghosted->getGlobalElement(c), neighbor_cells.size(), &vals(0,0), neighbor_cell_gids.data());
    CHECK(!ierr);
  }

  ierr |= matrix.FillComplete();
  CHECK(!ierr);

  control.fillComplete(graph->DomainMap(), graph->RangeMap());

  auto os = Teuchos::getFancyOStream(Teuchos::rcpFromRef(std::cout));
  matrix.Matrix()->describe(*os, Teuchos::VERB_EXTREME);
  control.describe(*os, Teuchos::VERB_EXTREME);
  
  CHECK(!ierr);

  // check matrix equality
  for (int c=0; c!=ncells; ++c) {
    const double *ctrl_vals, *mat_vals;
    const int *ctrl_inds, *mat_inds;
    int ctrl_nentries, mat_nentries;

    ierr |= matrix.ExtractMyRowCopy(c, mat_nentries, mat_vals, mat_inds);
    CHECK(!ierr);

    ierr |= control.getLocalRowView(c, ctrl_nentries, ctrl_vals, ctrl_inds);
    CHECK(!ierr);

    CHECK_EQUAL(ctrl_nentries, mat_nentries);
    for (int i=0; i!=std::min(ctrl_nentries, mat_nentries); ++i) {
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
  using namespace Amanzi::Operators;

  auto comm = getDefaultComm();
  int getRank = comm->getRank();

  if (getRank == 0)
    std::cout << "Test: FE like matrix, off-proc assembly" << std::endl;

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
  int ierr(0);
  Teuchos::RCP<GraphFE> graph =
    Teuchos::rcp(new GraphFE(face_map, face_map_ghosted, face_map_ghosted, 5));

  Entity_ID_View faces;
  Entity_ID_View face_cells;
  for (int c = 0; c != ncells; ++c) {
    mesh->cell_get_faces(c, faces);

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
  Matrix_type control(graph->Graph());

  for (int c = 0; c != ncells; ++c) {
    mesh->cell_get_faces(c, faces);

    Teuchos::SerialDenseVector<int,int> face_gids(faces.size());
    for (int n=0; n!=faces.size(); ++n) {
      face_gids[n] = face_map_ghosted->getGlobalElement(faces[n]);
      AMANZI_ASSERT(face_gids[n] > -1);
    }

    Teuchos::SerialDenseMatrix<int,double> vals(faces.size(), faces.size());
    vals.random();

    ierr |= matrix.SumIntoMyValues_Transposed(faces.data(), faces.data(), vals);
    CHECK(!ierr);

    for (int i=0; i!=faces.size(); ++i) {
      control.sumIntoGlobalValues(face_gids[i], faces.size(), &vals(0,i), &face_gids[0]);
    }
  }

  ierr |= matrix.FillComplete();
  control.fillComplete(graph->DomainMap(), graph->RangeMap());


  // check matrix equality
  int nfaces = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED);
  for (int f=0; f!=nfaces; ++f) {
    const double *ctrl_vals, *mat_vals;
    const int *ctrl_inds, *mat_inds;
    int ctrl_nentries, mat_nentries;

    ierr |= matrix.ExtractMyRowCopy(f, mat_nentries, mat_vals, mat_inds);
    CHECK(!ierr);

    ierr |= control.getLocalRowView(f, ctrl_nentries, ctrl_vals, ctrl_inds);
    CHECK(!ierr);

    CHECK_EQUAL(ctrl_nentries, mat_nentries);
    for (int i=0; i!=std::min(ctrl_nentries, mat_nentries); ++i) {
      CHECK_EQUAL(ctrl_inds[i], mat_inds[i]);
      CHECK_CLOSE(ctrl_vals[i], mat_vals[i], 1.e-10);
    }
  }
}
