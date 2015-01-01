/*
  This is the operator component of the Amanzi code. 

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

#include "MeshFactory.hh"
#include "GMVMesh.hh"

#include "GraphFE.hh"

/* *****************************************************************
 * this test is a null test -- all entries are local
* **************************************************************** */
TEST(FE_GRAPH_NEAREST_NEIGHBOR_TPFA) {
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;

  Epetra_MpiComm comm(MPI_COMM_WORLD);
  int MyPID = comm.MyPID();

  if (MyPID == 0) std::cout << "Test: 2D steady-state elliptic solver, mixed discretization" << std::endl;

  // read parameter list
  std::string xmlFileName = "test/operator_convergence.xml";
  Teuchos::ParameterXMLFileReader xmlreader(xmlFileName);
  Teuchos::ParameterList plist = xmlreader.getParameters();

  Amanzi::VerboseObject::hide_line_prefix = true;

  // create a mesh 
  Teuchos::ParameterList region_list = plist.get<Teuchos::ParameterList>("Regions");
  GeometricModelPtr gm = new GeometricModel(2, region_list, &comm);

  FrameworkPreference pref;
  pref.clear();
  pref.push_back(MSTK);

  MeshFactory meshfactory(&comm);
  meshfactory.preference(pref);
  Teuchos::RCP<Mesh> mesh = meshfactory(0.0, 0.0, 1.0, 1.0, 10, 10, gm);
  //  Teuchos::RCP<const Mesh> mesh = meshfactory("test/median32x33.exo", gm);

  // grab the maps
  int ncells = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  Teuchos::RCP<Epetra_Map> cell_map = Teuchos::rcp(new Epetra_Map(mesh->cell_map(false)));
  Teuchos::RCP<Epetra_Map> cell_map_ghosted = Teuchos::rcp(new Epetra_Map(mesh->cell_map(true)));

  // create the graph
  int ierr(0);
  GraphFE graph(cell_map, cell_map_ghosted, cell_map_ghosted, 5);
  
  Entity_ID_List faces;
  Entity_ID_List face_cells;
  std::vector<int> neighbor_cells;
  for (int c=0; c!=ncells; ++c) {
    neighbor_cells.resize(0);
    neighbor_cells.push_back(c);
    
    mesh->cell_get_faces(c, AmanziMesh::USED, &faces);
    for (int f=0; f!=faces.size(); ++f) {
      mesh->face_get_cells(f, AmanziMesh::USED, &face_cells);
      if (face_cells.size() > 1) {
	neighbor_cells.push_back(c == face_cells[0] ? face_cells[1] : face_cells[0]);
      }	
    }

    ierr |= graph.InsertMyIndices(c, neighbor_cells.size(), &neighbor_cells[0]);
    CHECK(!ierr);
  }

  ierr |= graph.FillCompete(cell_map, cell_map);
  CHECK(!ierr);
}


/* *****************************************************************
 * this test is a null test -- all entries are local
* **************************************************************** */
TEST(FE_GRAPH_FACE_FACE) {
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;

  Epetra_MpiComm comm(MPI_COMM_WORLD);
  int MyPID = comm.MyPID();

  if (MyPID == 0) std::cout << "Test: 2D steady-state elliptic solver, mixed discretization" << std::endl;

  // read parameter list
  std::string xmlFileName = "test/operator_convergence.xml";
  Teuchos::ParameterXMLFileReader xmlreader(xmlFileName);
  Teuchos::ParameterList plist = xmlreader.getParameters();

  Amanzi::VerboseObject::hide_line_prefix = true;

  // create a mesh 
  Teuchos::ParameterList region_list = plist.get<Teuchos::ParameterList>("Regions");
  GeometricModelPtr gm = new GeometricModel(2, region_list, &comm);

  FrameworkPreference pref;
  pref.clear();
  pref.push_back(MSTK);

  MeshFactory meshfactory(&comm);
  meshfactory.preference(pref);
  Teuchos::RCP<Mesh> mesh = meshfactory(0.0, 0.0, 1.0, 1.0, 10, 10, gm);
  //  Teuchos::RCP<const Mesh> mesh = meshfactory("test/median32x33.exo", gm);

  // grab the maps
  int ncells = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  Teuchos::RCP<Epetra_Map> face_map = Teuchos::rcp(new Epetra_Map(mesh->face_map(false)));
  Teuchos::RCP<Epetra_Map> face_map_ghosted = Teuchos::rcp(new Epetra_Map(mesh->face_map(true)));

  // create the graph
  int ierr(0);
  GraphFE graph(face_map, face_map_ghosted, face_map_ghosted, 5);
  
  Entity_ID_List faces;
  Entity_ID_List face_cells;
  for (int c=0; c!=ncells; ++c) {
    mesh->cell_get_faces(c, AmanziMesh::USED, &faces);

    for (int n=0; n!=faces.size(); ++n) {
      ierr |= graph.InsertMyIndices(faces[n], faces.size(), &faces[0]);
      CHECK(!ierr);
    }
  }

  ierr |= graph.FillCompete(face_map, face_map);
  CHECK(!ierr);

}
