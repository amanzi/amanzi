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

#include "MeshFactory.hh"
#include "MeshAudit.hh"
#include "GMVMesh.hh"

#include "State.hh"
#include "Richards_PK.hh"

/* **************************************************************** */
TEST(FLOW_2D_RICHARDS) {
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Flow;

  Epetra_MpiComm comm(MPI_COMM_WORLD);
  int MyPID = comm.MyPID();

  if (MyPID == 0) std::cout << "Test: 2D Richards, 2-layer model" << std::endl;

  /* read parameter list */
  std::string xmlFileName = "test/flow_richards_2D.xml";
  ParameterXMLFileReader xmlreader(xmlFileName);
  ParameterList plist = xmlreader.getParameters();

  /* create an MSTK mesh framework */
  ParameterList region_list = plist.get<Teuchos::ParameterList>("Regions");
  GeometricModelPtr gm = new GeometricModel(2, region_list, &comm);

  FrameworkPreference pref;
  pref.clear();
  pref.push_back(MSTK);

  MeshFactory meshfactory(&comm);
  meshfactory.preference(pref);
  RCP<const Mesh> mesh = meshfactory(0.0, -2.0, 1.0, 0.0, 18, 18, gm);

  /* create a simple state and populate it */
  Amanzi::VerboseObject::hide_line_prefix = false;

  ParameterList state_list;
  RCP<State> S = rcp(new State(state_list));
  S->RegisterDomainMesh(rcp_const_cast<Mesh>(mesh));

  Richards_PK* RPK = new Richards_PK(plist, S);
  S->Setup();
  S->InitializeFields();
  S->InitializeEvaluators();
  RPK->InitializeFields();
  S->CheckAllFieldsInitialized();

  /* modify the default state for the problem at hand */
  std::string passwd("state"); 
  Epetra_MultiVector& K = *S->GetFieldData("permeability", passwd)->ViewComponent("cell");
  
  AmanziMesh::Entity_ID_List block;
  mesh->get_set_entities("Material 1", AmanziMesh::CELL, AmanziMesh::OWNED, &block);
  for (int i = 0; i != block.size(); ++i) {
    int c = block[i];
    K[0][c] = 0.1;
    K[1][c] = 2.0;
  }

  mesh->get_set_entities("Material 2", AmanziMesh::CELL, AmanziMesh::OWNED, &block);
  for (int i = 0; i != block.size(); ++i) {
    int c = block[i];
    K[0][c] = 0.5;
    K[1][c] = 0.5;
  }

  *S->GetScalarData("fluid_density", passwd) = 1.0;
  *S->GetScalarData("fluid_viscosity", passwd) = 1.0;
  Epetra_Vector& gravity = *S->GetConstantVectorData("gravity", passwd);
  gravity[1] = -1.0;

  /* create the initial pressure function */
  Epetra_MultiVector& p = *S->GetFieldData("pressure", passwd)->ViewComponent("cell");

  for (int c = 0; c < p.MyLength(); c++) {
    const Point& xc = mesh->cell_centroid(c);
    p[0][c] = xc[1] * (xc[1] + 2.0);
  }

  /* initialize the Richards process kernel */
  RPK->Initialize(S.ptr());
  RPK->ti_specs_sss().T1 = 100.0;
  RPK->ti_specs_sss().max_itrs = 400;

  RPK->InitializeAuxiliaryData();
  RPK->InitSteadyState(0.0, 1e-8);

  /* solve the problem */
  RPK->AdvanceToSteadyState(0.0, 0.1);
  RPK->CommitState(0.0, S.ptr());

  if (MyPID == 0) {
    GMV::open_data_file(*mesh, (std::string)"flow.gmv");
    GMV::start_data();
    GMV::write_cell_data(p, 0, "pressure");
    GMV::close_data_file();
  }

  /* check the pressure */
  int ncells = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  for (int c = 0; c < ncells; c++) CHECK(p[0][c] > 0.0 && p[0][c] < 2.1);

  delete RPK;
}
