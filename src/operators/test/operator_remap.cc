/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <vector>

#include "UnitTest++.h"

#include "Epetra_MultiVector.h"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"

#include "GMVMesh.hh"
#include "MeshFactory.hh"

#include "AdvectionRiemann.hh"
#include "OperatorDefs.hh"


/* *****************************************************************
* Remap of polynomilas in two dimensions
***************************************************************** */
TEST(REMAP_CONSTANT_2D) {
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;

  Epetra_MpiComm comm(MPI_COMM_WORLD);
  int MyPID = comm.MyPID();

  if (MyPID == 0) std::cout << "\nTest: remap of constant functions in 2D." << std::endl;

  // create initial mesh
  FrameworkPreference pref;
  pref.clear();
  pref.push_back(MSTK);

  MeshFactory meshfactory(&comm);
  meshfactory.preference(pref);

  Teuchos::RCP<const Mesh> mesh1 = meshfactory(0.0, 0.0, 1.0, 1.0, 7, 7);

  int ncells_owned = mesh1->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  int ncells_wghost = mesh1->num_entities(AmanziMesh::CELL, AmanziMesh::USED);

  // create deformed mesh
  Teuchos::RCP<Mesh> mesh2 = meshfactory(0.0, 0.0, 1.0, 1.0, 7, 7);

  int nnodes_owned = mesh2->num_entities(AmanziMesh::NODE, AmanziMesh::OWNED);

  AmanziGeometry::Point xv(2);
  Entity_ID_List nodeids;
  AmanziGeometry::Point_List new_positions, final_positions;

  for (int v = 0; v < nnodes_owned; ++v) {
    mesh2->node_get_coordinates(v, &xv);
    nodeids.push_back(v);
    new_positions.push_back(xv);
  }
    
  mesh2->deform(nodeids, new_positions, false, &final_positions);

  // create and initialize cell-based field 
  Teuchos::RCP<Epetra_MultiVector> p1 = Teuchos::rcp(new Epetra_MultiVector(mesh1->cell_map(true), 1));

  for (int c = 0; c < ncells_wghost; c++) {
    const AmanziGeometry::Point& xc = mesh1->cell_centroid(c);
    (*p1)[0][c] = xc[0] + 2 * xc[1];
  }

  // remap cell-based field
  Teuchos::RCP<Epetra_MultiVector> p2 = Teuchos::rcp(new Epetra_MultiVector(mesh2->cell_map(true), 1));

  Teuchos::ParameterList plist;
  plist.set<std::string>("discretization", "DG order 0: face");

  plist.sublist("schema domain")
      .set<std::string>("base", "face")
      .set<Teuchos::Array<std::string> >("location", std::vector<std::string>({"cell"}))
      .set<Teuchos::Array<std::string> >("type", std::vector<std::string>({"scalar"}))
      .set<Teuchos::Array<int> >("number", std::vector<int>({1}));

  plist.sublist("schema range") = plist.sublist("schema domain");

  Teuchos::RCP<AdvectionRiemann> op = Teuchos::rcp(new AdvectionRiemann(plist, mesh1));
  

  // calculate error
}



