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

// TPLs
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "UnitTest++.h"

// Amanzi
#include "GMVMesh.hh"
#include "MeshAudit.hh"
#include "MeshFactory.hh"
#include "State.hh"

// Flow
#include "Richards_PK.hh"
#include "Richards_SteadyState.hh"

/* **************************************************************** */
TEST(FLOW_BOUNDARY_SOLVER) {
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Flow;

  Epetra_MpiComm comm(MPI_COMM_WORLD);
  int MyPID = comm.MyPID();
  if (MyPID == 0) std::cout << "Test: Richards, boundary solver" << std::endl;

  // read parameter list
  std::string xmlFileName = "test/flow_richards_boundary_solver.xml";
  Teuchos::RCP<Teuchos::ParameterList> plist = Teuchos::getParametersFromXmlFile(xmlFileName);

  // create a mesh framework
  Teuchos::ParameterList region_list = plist->get<Teuchos::ParameterList>("Regions");
  GeometricModelPtr gm = new GeometricModel(2, region_list, &comm);

  FrameworkPreference pref;
  pref.clear();
  pref.push_back(MSTK);
  pref.push_back(STKMESH);

  double bottom = -10.;
  MeshFactory meshfactory(&comm);
  meshfactory.preference(pref);
  Teuchos::RCP<const Mesh> mesh = meshfactory(0.0, bottom, 1.0, 0.0, 1, 10, gm);

  // create a simple state and populate it
  Teuchos::ParameterList state_list = plist->sublist("State");
  Teuchos::RCP<State> S = Teuchos::rcp(new State(state_list));
  S->RegisterDomainMesh(Teuchos::rcp_const_cast<Mesh>(mesh));

  Teuchos::RCP<TreeVector> soln = Teuchos::rcp(new TreeVector());
  Teuchos::RCP<Richards_PK> RPK = Teuchos::rcp(new Richards_PK(plist, "Flow", S, soln));

  RPK->Setup();
  S->Setup();
  S->InitializeFields();
  S->InitializeEvaluators();

  // modify the default state for the problem at hand
  std::string passwd("flow"); 
  Epetra_MultiVector& K = *S->GetFieldData("permeability", passwd)->ViewComponent("cell");
  
  AmanziMesh::Entity_ID_List block;
  mesh->get_set_entities("Material 1", AmanziMesh::CELL, AmanziMesh::OWNED, &block);
  for (int i = 0; i != block.size(); ++i) {
    int c = block[i];
    K[0][c] = 1e-13;
    K[1][c] = 1e-13;
  }

  S->GetField("permeability", "flow")->set_initialized();

  double atm_pressure = 101325.;
  double rho = *S->GetScalarData("fluid_density", passwd);

  Epetra_Vector& gravity = *S->GetConstantVectorData("gravity", "state");
  gravity[1] = -9.8;
  S->GetField("gravity", "state")->set_initialized();

  // create the initial pressure function
  Epetra_MultiVector& p = *S->GetFieldData("pressure", passwd)->ViewComponent("cell");

  for (int c = 0; c < p.MyLength(); c++) {
    const Point& xc = mesh->cell_centroid(c);
    p[0][c] = atm_pressure + rho*gravity[1]*(xc[1] - bottom);
  }

   
  // initialize the Richards process kernel
  RPK->Initialize();
  S->CheckAllFieldsInitialized();
  
  double boundary_val;

  int nfaces = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::OWNED);
  // for (int f = 0; f < nfaces; f++) {
  //   AmanziMesh::Entity_ID_List cells;
  //   mesh->face_get_cells(f, AmanziMesh::USED, &cells);
  //   if ((cells.size() == 1)&&(cells[0] == 9)){
  int f = 29;
  boundary_val = RPK->BoundaryFaceValue(f, *S->GetFieldData("pressure", passwd));
  std::cout<<": "<<f<<" "<<" "<< boundary_val<<"\n";

  //   }
  // }

  //exit(0);


  // std::cout<<p<<"\n";
  // exit(0);

  // // solve the problem 
  // TI_Specs ti_specs;
  // ti_specs.T0 = 0.0;
  // ti_specs.dT0 = 1.0;
  // ti_specs.T1 = 1.0;
  // ti_specs.max_itrs = 400;

  // AdvanceToSteadyState(S, *RPK, ti_specs, soln);
  // RPK->CommitStep(0.0, 1.0);  // dummy times



  if (MyPID == 0) {
    GMV::open_data_file(*mesh, (std::string)"flow.gmv");
    GMV::start_data();
    GMV::write_cell_data(p, 0, "pressure");
    GMV::close_data_file();
  }

  // check the pressure
  // int ncells = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  // for (int c = 0; c < ncells; c++) CHECK(p[0][c] > -4.0 && p[0][c] < 0.01);

  CHECK (fabs(78341.9 - boundary_val) < 1e-1);

}

