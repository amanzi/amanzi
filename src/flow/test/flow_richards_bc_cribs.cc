/*
  This is the flow component of the Amanzi code. 

  Copyright 2010-2012 held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
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

#include "Mesh.hh"
#include "MeshFactory.hh"
#include "GMVMesh.hh"

#include "State.hh"
#include "Richards_PK.hh"


/* **************************************************************** */
TEST(FLOW_3D_RICHARDS) {
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::AmanziFlow;

  Epetra_MpiComm comm(MPI_COMM_WORLD);
  int MyPID = comm.MyPID();
  if (MyPID == 0) cout << "Test: 3D Richards, crib model" << endl;

  /* read parameter list */
  string xmlFileName = "test/flow_richards_bc_cribs.xml";
  ParameterXMLFileReader xmlreader(xmlFileName);
  ParameterList plist = xmlreader.getParameters();

  /* create a mesh framework */
  ParameterList region_list = plist.get<Teuchos::ParameterList>("Regions");
  GeometricModelPtr gm = new GeometricModel(3, region_list, &comm);

  FrameworkPreference pref;
  pref.clear();
  pref.push_back(MSTK);

  MeshFactory factory(&comm);
  factory.preference(pref);  
  ParameterList mesh_list = plist.get<ParameterList>("Mesh").get<ParameterList>("Unstructured");
  ParameterList factory_list = mesh_list.get<ParameterList>("Generate Mesh");
  Teuchos::RCP<Mesh> mesh(factory(factory_list, gm));

  /* create a simple state and populate it */
  Amanzi::VerboseObject::hide_line_prefix = true;

  ParameterList state_list = plist.get<ParameterList>("State");
  RCP<State> S = rcp(new State(state_list));
  S->RegisterDomainMesh(rcp_const_cast<Mesh>(mesh));

  Richards_PK* RPK = new Richards_PK(plist, S);
  S->Setup();
  S->InitializeFields();
  RPK->InitializeFields();
  S->CheckAllFieldsInitialized();

  /* initialize the Richards process kernel */
  RPK->InitPK();
  RPK->InitializeAuxiliaryData();
  RPK->InitSteadyState(0.0, 1e-7);

  /* solve the problem */
  RPK->AdvanceToSteadyState(0.0, 1e-7);
  RPK->CommitState(S);

  /* derive dependent variable */
  const Epetra_MultiVector& p = *S->GetFieldData("pressure")->ViewComponent("cell");
  const Epetra_MultiVector& ws = *S->GetFieldData("water_saturation")->ViewComponent("cell");
  const Epetra_MultiVector& K = *S->GetFieldData("permeability")->ViewComponent("cell");

  GMV::open_data_file(*mesh, (std::string)"flow.gmv");
  GMV::start_data();
  GMV::write_cell_data(p, 0, "pressure");
  GMV::write_cell_data(ws, 0, "saturation");
  GMV::write_cell_data(*K(0), 0, "permeability_x");
  GMV::write_cell_data(*K(1), 0, "permeability_y");
  GMV::write_cell_data(*K(2), 0, "permeability_z");
  GMV::close_data_file();

  /* check the pressure profile */
  int ncells = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  //for (int c = 0; c < ncells; c++) cout << (mesh->cell_centroid(c))[2] << " " << pressure[c] << endl;
  for (int c = 0; c < ncells; c++) CHECK(p[0][c] > 4500.0 && p[0][c] < 101325.0);

  delete RPK;
}
