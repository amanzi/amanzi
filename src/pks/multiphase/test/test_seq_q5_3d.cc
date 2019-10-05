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
  //std::string xmlFileName = "test/test_seq_q5_3d.xml";
  std::string xmlFileName = "test/test_seq_spe10.xml";
  ParameterXMLFileReader xmlreader(xmlFileName);
  ParameterList plist = xmlreader.getParameters();
  Teuchos::RCP<Teuchos::ParameterList> parameter_list = 
     Teuchos::rcp(new Teuchos::ParameterList(plist));
  Teuchos::RCP<Teuchos::ParameterList> pk_tree_list_S = Teuchos::sublist(parameter_list, "PK Tree Saturation");
  Teuchos::RCP<Teuchos::ParameterList> pk_tree_list_P = Teuchos::sublist(parameter_list, "PK Tree Pressure");
  Teuchos::RCP<Teuchos::ParameterList> state_list = Teuchos::sublist(parameter_list, "State");
  Teuchos::RCP<Teuchos::ParameterList> time_list = Teuchos::sublist(parameter_list, "Cycle Driver");
  //ParameterList mesh_list = plist.sublist("Mesh").sublist("Unstructured").sublist("Generate Mesh");
  ParameterList domain_list = plist.sublist("Domain");
  //std::vector<int> mesh_size = mesh_list.get<Teuchos::Array<int> > ("Number of Cells").toVector();
  //std::vector<double> high_coor = mesh_list.get<Teuchos::Array<double> > ("Domain High Coordinate").toVector();
  //std::vector<double> low_coor = mesh_list.get<Teuchos::Array<double> > ("Domain Low Coordinate").toVector();
  std::string file_name = plist.sublist("Mesh").sublist("Unstructured").sublist("Read Mesh File").get<std::string>("File");
  int geom_dim = domain_list.get<int> ("Spatial Dimension");

  /* create a MSTK mesh framework */
  ParameterList region_list = plist.get<Teuchos::ParameterList>("Regions");
  GeometricModelPtr gm = new GeometricModel(geom_dim, region_list, &comm);

  FrameworkPreference pref;
  pref.clear();
  pref.push_back(MSTK);
  pref.push_back(STKMESH);

  MeshFactory meshfactory(&comm);
  meshfactory.preference(pref);
  //RCP<const Mesh> mesh = meshfactory(low_coor[0], low_coor[1], low_coor[2], 
  //                                   high_coor[0], high_coor[1], high_coor[2], 
  //                                   mesh_size[0], mesh_size[1], mesh_size[2], gm);
  RCP<const Mesh> mesh = meshfactory(file_name, gm);

  /* create a simple state populate it */
  RCP<State> S = rcp(new State(*state_list));
  S->RegisterDomainMesh(rcp_const_cast<Mesh>(mesh));
  std::string passwd("state");

  std::vector<Teuchos::RCP<Visualization> > visualization_;
  for (State::mesh_iterator mesh=S->mesh_begin();
       mesh!=S->mesh_end(); ++mesh) {
    if (parameter_list->isSublist("Visualization Data")) {
      Teuchos::ParameterList& vis_plist = parameter_list->sublist("Visualization Data");
      Teuchos::RCP<Visualization> vis =
        Teuchos::rcp(new Visualization(vis_plist, &comm));
      vis->set_mesh(mesh->second.first);
      vis->CreateFiles();
      visualization_.push_back(vis);
    }
  }

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
  
  // write initial state
  for (std::vector<Teuchos::RCP<Visualization> >::iterator vis=visualization_.begin();
       vis!=visualization_.end(); ++vis) {
    WriteVis((*vis).ptr(), S.ptr());
  }

  // Init the PKs
  //double dT, dT_actual, t_final, t_sim;
  //t_final = 0.1;
  //t_sim = 0.0;
  //dT_actual = 0.0;
  //dT = 0.05;
  PPK->Initialize();
  SPK->Initialize();
  S->CheckAllFieldsInitialized();

  bool failed = true;
  double t_sim = 0.0;
  double dT = time_list->get<double>("time step");
  double t_final = time_list->get<double>("end time");
  while (!((t_sim > t_final) || abs(t_sim - t_final) < 1.0e-8)) {
    PPK->AdvanceStep(t_sim, t_sim + dT, false);
    PPK->CommitState(S.ptr());
    failed = SPK->AdvanceStep(t_sim, t_sim + dT, false);
    if (!failed) {
      t_sim += dT;
      SPK->CommitState(S.ptr());
      S->advance_time(dT);
      S->advance_cycle();
      for (std::vector<Teuchos::RCP<Visualization> >::iterator vis=visualization_.begin();
           vis!=visualization_.end(); ++vis) {
        WriteVis((*vis).ptr(), S.ptr());
      } 
    } else {
      dT = dT / 2.0;
    }
    std::cout << "t_sim: " << t_sim << "\n";
  }

  /*
  // get the data for viz
  Epetra_MultiVector& ws_c = *S->GetFieldData("water_saturation", passwd)->ViewComponent("cell", false);
  Epetra_MultiVector& p = *S->GetFieldData("pressure", passwd)->ViewComponent("cell", false);
  if (MyPID == 0) {
    GMV::open_data_file(*mesh, (std::string)"test_seq_q5_3d.gmv");
    GMV::start_data();
    GMV::write_cell_data(ws_c, 0, "water_saturation");
    GMV::write_cell_data(p, 0, "pressure");
    GMV::close_data_file();
  }
  */

  delete PPK;
  delete SPK;
}
