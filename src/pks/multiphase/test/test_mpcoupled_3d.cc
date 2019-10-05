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
#include "TimerManager.hh"

#include "State.hh"
#include "MPCoupled_PK.hh"
#include "OperatorDefs.hh"
#include "Visualization.hh"


/* **************************************************************** */
TEST(MULTIPHASE_MPCOUPLED_3D) {
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;

  Epetra_MpiComm comm(MPI_COMM_WORLD);
  int MyPID = comm.MyPID();

  if (MyPID == 0) std::cout << "Test: multiphase coupled, uniform rectangular mesh" << std::endl;

  /* read parameter list */
  std::string xmlFileName = "test/weak_scaling_copper.xml";
  //std::string xmlFileName = "test/spe10_3d.xml";
  ParameterXMLFileReader xmlreader(xmlFileName);
  ParameterList plist = xmlreader.getParameters();
  Teuchos::RCP<Teuchos::ParameterList> parameter_list = 
     Teuchos::rcp(new Teuchos::ParameterList(plist));
  Teuchos::RCP<Teuchos::ParameterList> pk_tree_list = Teuchos::sublist(parameter_list, "PK Tree");
  Teuchos::RCP<Teuchos::ParameterList> state_list = Teuchos::sublist(parameter_list, "State");
  Teuchos::RCP<Teuchos::ParameterList> time_list = Teuchos::sublist(parameter_list, "Cycle Driver");
  ParameterList mesh_list = plist.sublist("Mesh").sublist("Unstructured").sublist("Generate Mesh");
  std::vector<int> mesh_size = mesh_list.get<Teuchos::Array<int> > ("Number of Cells").toVector();
  std::vector<double> high_coor = mesh_list.get<Teuchos::Array<double> > ("Domain High Coordinate").toVector();
  std::vector<double> low_coor = mesh_list.get<Teuchos::Array<double> > ("Domain Low Coordinate").toVector();
  //std::string file_name = plist.sublist("Mesh").sublist("Unstructured").sublist("Read Mesh File").get<std::string>("File");
  //std::cout << file_name << endl;

  /* create a MSTK mesh framework */
  ParameterList region_list = plist.get<Teuchos::ParameterList>("Regions");
  GeometricModelPtr gm = new GeometricModel(3, region_list, &comm);

  FrameworkPreference pref;
  pref.clear();
  pref.push_back(MSTK);
  pref.push_back(STKMESH);

  MeshFactory meshfactory(&comm);
  meshfactory.preference(pref);
  RCP<const Mesh> mesh = meshfactory(low_coor[0], low_coor[1], low_coor[2],
                                     high_coor[0], high_coor[1], high_coor[2],
                                     mesh_size[0], mesh_size[1], mesh_size[2], gm);
  //RCP<const Mesh> mesh = meshfactory(file_name, gm);
  //int ncells = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  //int nfaces_wghost = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::USED);

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
  
  // Create Multiphase PK
  Teuchos::RCP<TreeVector> soln = Teuchos::rcp(new TreeVector());
  Multiphase::MPCoupled_PK* MPK = new Multiphase::MPCoupled_PK(*pk_tree_list, parameter_list, S, soln);
  S->Setup();
  S->InitializeFields();

  // Init the PK and solve for pressure
  MPK->Initialize();
  
  // write initial state
  for (std::vector<Teuchos::RCP<Visualization> >::iterator vis=visualization_.begin();
       vis!=visualization_.end(); ++vis) {
    WriteVis((*vis).ptr(), S.ptr());
  }
  
  bool failed = true;
  //double dT = 0.1;
  double t_sim = 0.0;
  //double t_final = 1.0;
  double dT = time_list->get<double>("time step");
  double t_final = time_list->get<double>("end time");
  while (!(t_sim > t_final || abs(t_sim - t_final) < 1.0e-10)) {
    failed = MPK->AdvanceStep(t_sim, t_sim + dT, false);
    if (failed) {
      dT /= 2.0;
    } else {
      t_sim += dT;
      MPK->CommitStep(t_sim, t_sim + dT);
      S->advance_time(dT);
      S->advance_cycle();
      std::cout << "State time: " << S->time() << "; State cycle: " << S->cycle() << endl;
      for (std::vector<Teuchos::RCP<Visualization> >::iterator vis=visualization_.begin();
           vis!=visualization_.end(); ++vis) {
        WriteVis((*vis).ptr(), S.ptr());
        std::cout << "writing visualization file" << std::endl;
      } 
    } 
    std::cout << "t_sim: " << t_sim << "\n";
  }
  //MPI_Comm mpi_comm(MPI_COMM_WORLD);
  //Amanzi::timer_manager.parSync(mpi_comm);
  //Amanzi::timer_manager.print();
  
  delete MPK;
}