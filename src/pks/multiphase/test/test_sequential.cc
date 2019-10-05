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
#include "Mesh.hh"
#include "GMVMesh.hh"

#include "State.hh"
#include "Saturation_PK.hh"
#include "Pressure_PK.hh"
#include "OperatorDefs.hh"


/* **************************************************************** */
TEST(MULTIPHASE_SATURATION) {
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;

  Epetra_MpiComm comm(MPI_COMM_WORLD);
  int MyPID = comm.MyPID();

  if (MyPID == 0) std::cout << "Test: steady state 2-phase immiscible flow, uniform rectangular mesh" << std::endl;

  /* read parameter list */
  std::string xmlFileName = "test/test_sequential.xml";
  ParameterXMLFileReader xmlreader(xmlFileName);
  ParameterList plist = xmlreader.getParameters();
  Teuchos::RCP<Teuchos::ParameterList> parameter_list = 
     Teuchos::rcp(new Teuchos::ParameterList(plist));
  Teuchos::RCP<Teuchos::ParameterList> pk_tree_list_S = Teuchos::sublist(parameter_list, "PK Tree Saturation");
  Teuchos::RCP<Teuchos::ParameterList> pk_tree_list_P = Teuchos::sublist(parameter_list, "PK Tree Pressure");
  Teuchos::RCP<Teuchos::ParameterList> state_list = Teuchos::sublist(parameter_list, "State");
  Teuchos::RCP<Teuchos::ParameterList> time_list = Teuchos::sublist(parameter_list, "Cycle Driver");
  ParameterList mesh_list = plist.sublist("Mesh").sublist("Unstructured").sublist("Generate Mesh");
  std::vector<int> mesh_size = mesh_list.get<Teuchos::Array<int> > ("Number of Cells").toVector(); 

  /* create a MSTK mesh framework */
  ParameterList region_list = plist.get<Teuchos::ParameterList>("Regions");
  GeometricModelPtr gm = new GeometricModel(2, region_list, &comm);

  FrameworkPreference pref;
  pref.clear();
  pref.push_back(MSTK);
  pref.push_back(STKMESH);

  MeshFactory meshfactory(&comm);
  meshfactory.preference(pref);
  RCP<const Mesh> mesh = meshfactory(0.0, 0.0, 1.0, 0.1, mesh_size[0], 1, gm);
  int ncells = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  int nfaces_wghost = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::ALL);

  /* create a simple state populate it */
  RCP<State> S = rcp(new State(*state_list));
  //RCP<State> S = rcp(new State());
  S->RegisterDomainMesh(rcp_const_cast<Mesh>(mesh));
  std::string passwd("state");

  // create the global solution vector
  Teuchos::RCP<TreeVector> pk_sol_s = Teuchos::rcp(new TreeVector());
  Teuchos::RCP<TreeVector> pk_sol_p = Teuchos::rcp(new TreeVector());

  // Create pressure and saturation PK
  Multiphase::Saturation_PK* SPK = new Multiphase::Saturation_PK(*pk_tree_list_S, parameter_list, S, pk_sol_s);
  Multiphase::Pressure_PK* PPK = new Multiphase::Pressure_PK(*pk_tree_list_P, parameter_list, S, pk_sol_p);
  PPK->Setup();
  SPK->Setup();

  // Set up and initialize state
  S->Setup();
  S->InitializeFields();

  // Init the PKs
  double dT, t_final, t_sim;
  //t_final = 0.6;
  //dT = 0.1;
  //t_final = 0.1;
  dT = time_list->get<double>("time step");
  t_final = time_list->get<double>("end time");
  t_sim = 0.0;
  //dT = 0.3/mesh_size[0];
  PPK->Initialize();
  SPK->Initialize();
  S->CheckAllFieldsInitialized();

  bool failed = true;
  while (!(t_sim > t_final || abs(t_sim - t_final) < 1e-8)) {
    PPK->AdvanceStep(t_sim, t_sim + dT, false);
    PPK->CommitState(S.ptr());
    failed = SPK->AdvanceStep(t_sim, t_sim + dT, false);
    if (!failed) {
      t_sim += dT;
      SPK->CommitState(S.ptr());
    } else {
      dT = dT / 2.0;
    }
  }
  std::cout << "t_sim: " << t_sim << "\n";
  
  // get the data for viz
  Epetra_MultiVector& ws_c = *S->GetFieldData("water_saturation", passwd)->ViewComponent("cell", false);
  Epetra_MultiVector& p = *S->GetFieldData("pressure", passwd)->ViewComponent("cell", false);
  if (MyPID == 0) {
    GMV::open_data_file(*mesh, (std::string)"test_sequential.gmv");
    GMV::start_data();
    GMV::write_cell_data(ws_c, 0, "water_saturation");
    GMV::write_cell_data(p, 0, "pressure");
    GMV::close_data_file();
  }

  delete PPK;
  delete SPK;
}
