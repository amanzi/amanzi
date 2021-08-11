/*
  Data Structures

  Copyright 2010-2012 held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@lanl.gov)
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

#include "AmanziComm.hh"
#include "VerboseObject.hh"
#include "Mesh.hh"
#include "MeshFactory.hh"
#include "GraphFE.hh"

/* *****************************************************************
* this test is a null test -- all entries are local
***************************************************************** */
TEST(FE_GRAPH_NEAREST_NEIGHBOR_TPFA) {
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;

  auto comm = getDefaultComm();
  int MyPID = comm->MyPID();

  if (MyPID == 0) std::cout << "Test: FD like graph, null off-proc assembly" << std::endl;

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

  MeshFactory meshfactory(comm,gm);
  meshfactory.set_preference(pref);
  Teuchos::RCP<Mesh> mesh = meshfactory.create(0.0, 0.0, 1.0, 1.0, 10, 10);
  //  Teuchos::RCP<const Mesh> mesh = meshfactory.create("test/median32x33.exo");

  // grab the maps
  int ncells = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  Teuchos::RCP<Epetra_Map> cell_map = Teuchos::rcp(new Epetra_Map(mesh->cell_map(false)));
  Teuchos::RCP<Epetra_Map> cell_map_ghosted = Teuchos::rcp(new Epetra_Map(mesh->cell_map(true)));

  // create the graphs, one to test local, the other to test global insertion
  int ierr(0);
  GraphFE graph_local(cell_map, cell_map_ghosted, cell_map_ghosted, 5);
  GraphFE graph_global(cell_map, cell_map_ghosted, cell_map_ghosted, 5);
  
  Entity_ID_List faces;
  Entity_ID_List face_cells;
  std::vector<int> neighbor_cells;
  for (int c=0; c!=ncells; ++c) {
    neighbor_cells.resize(0);
    neighbor_cells.push_back(c);
    
    mesh->cell_get_faces(c, &faces);
    for (int n=0; n!=faces.size(); ++n) {
      mesh->face_get_cells(faces[n], AmanziMesh::Parallel_type::ALL, &face_cells);
      if (face_cells.size() > 1) {
        neighbor_cells.push_back(c == face_cells[0] ? face_cells[1] : face_cells[0]);
      }	
    }

    int global_c = cell_map->GID(c);
    std::vector<int> global_neighbors(neighbor_cells.size());
    for (int n=0; n!=neighbor_cells.size(); ++n) global_neighbors[n] = cell_map_ghosted->GID(neighbor_cells[n]);

    ierr |= graph_local.InsertMyIndices(c, neighbor_cells.size(), &neighbor_cells[0]);
    CHECK(!ierr);
    ierr |= graph_global.InsertGlobalIndices(global_c, neighbor_cells.size(), &global_neighbors[0]);
    CHECK(!ierr);
  }

  ierr |= graph_local.FillComplete(cell_map, cell_map);
  CHECK(!ierr);
  ierr |= graph_global.FillComplete(cell_map, cell_map);
  CHECK(!ierr);
}


/* *****************************************************************
* this test is a null test -- all entries are local
***************************************************************** */
TEST(FE_GRAPH_FACE_FACE) {
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;

  auto comm = getDefaultComm();
  int MyPID = comm->MyPID();

  if (MyPID == 0) std::cout << "Test: FE like graph, off-proc assembly" << std::endl;

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

  MeshFactory meshfactory(comm,gm);
  meshfactory.set_preference(pref);
  Teuchos::RCP<Mesh> mesh = meshfactory.create(0.0, 0.0, 1.0, 1.0, 10, 10);
  //  Teuchos::RCP<const Mesh> mesh = meshfactory.create("test/median32x33.exo");

  // grab the maps
  int ncells = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  Teuchos::RCP<Epetra_Map> face_map = Teuchos::rcp(new Epetra_Map(mesh->face_map(false)));
  Teuchos::RCP<Epetra_Map> face_map_ghosted = Teuchos::rcp(new Epetra_Map(mesh->face_map(true)));

  // create the graph
  int ierr(0);
  GraphFE graph_local(face_map, face_map_ghosted, face_map_ghosted, 5);
  GraphFE graph_global(face_map, face_map_ghosted, face_map_ghosted, 5);
  
  Entity_ID_List faces;
  Entity_ID_List face_cells;
  for (int c=0; c!=ncells; ++c) {
    mesh->cell_get_faces(c, &faces);

    std::vector<int> global_faces(faces.size());
    for (int n=0; n!=faces.size(); ++n) global_faces[n] = face_map_ghosted->GID(faces[n]);
    
    for (int n=0; n!=faces.size(); ++n) {
      ierr |= graph_local.InsertMyIndices(faces[n], faces.size(), &faces[0]);
      CHECK(!ierr);
      AMANZI_ASSERT(global_faces[n] >= 0);
      ierr |= graph_global.InsertGlobalIndices(global_faces[n], global_faces.size(), &global_faces[0]);
      CHECK(!ierr);
    }
  }

  ierr |= graph_local.FillComplete(face_map, face_map);
  CHECK(!ierr);
  ierr |= graph_global.FillComplete(face_map, face_map);
  CHECK(!ierr);
}
