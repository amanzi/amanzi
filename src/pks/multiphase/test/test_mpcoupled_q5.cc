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

double rss_usage() { // return ru_maxrss in MBytes
#if (defined(__unix__) || defined(__unix) || defined(unix) || defined(__APPLE__) || defined(__MACH__))
  struct rusage usage;
  getrusage(RUSAGE_SELF, &usage);
#if (defined(__APPLE__) || defined(__MACH__))
  return static_cast<double>(usage.ru_maxrss)/1024.0/1024.0;
#else
  return static_cast<double>(usage.ru_maxrss)/1024.0;
#endif
#else
  return 0.0;
#endif
}


/* **************************************************************** */
TEST(MULTIPHASE_MPCOUPLED) {
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;

  Epetra_MpiComm comm(MPI_COMM_WORLD);
  int MyPID = comm.MyPID();

  if (MyPID == 0) std::cout << "Test: multiphase coupled, uniform rectangular mesh" << std::endl;

  /* read parameter list */
  //std::string xmlFileName = "test/test_eigenvalues.xml";
  std::string xmlFileName = "test/spe10_2d.xml";
  //std::string xmlFileName = "test/test_mpcoupled_grav.xml";
  //std::string xmlFileName = "test/test_mpcoupled_q5.xml";
  ParameterXMLFileReader xmlreader(xmlFileName);
  ParameterList plist = xmlreader.getParameters();
  Teuchos::RCP<Teuchos::ParameterList> parameter_list = 
     Teuchos::rcp(new Teuchos::ParameterList(plist));
  Teuchos::RCP<Teuchos::ParameterList> pk_tree_list = Teuchos::sublist(parameter_list, "PK Tree");
  Teuchos::RCP<Teuchos::ParameterList> state_list = Teuchos::sublist(parameter_list, "State");
  Teuchos::RCP<Teuchos::ParameterList> time_list = Teuchos::sublist(parameter_list, "Cycle Driver");
  Teuchos::RCP<Teuchos::ParameterList> flow_list = Teuchos::sublist(parameter_list, "Flow");
  ParameterList mesh_list = plist.sublist("Mesh").sublist("Unstructured").sublist("Generate Mesh");
  std::vector<int> mesh_size = mesh_list.get<Teuchos::Array<int> > ("Number of Cells").toVector();
  std::vector<double> high_coor = mesh_list.get<Teuchos::Array<double> > ("Domain High Coordinate").toVector();
  std::vector<double> low_coor = mesh_list.get<Teuchos::Array<double> > ("Domain Low Coordinate").toVector();
  std::string file_name = plist.sublist("Mesh").sublist("Unstructured").sublist("Read Mesh File").get<std::string>("File");
  std::cout << file_name << endl;

  /* create a MSTK mesh framework */
  ParameterList region_list = plist.get<Teuchos::ParameterList>("Regions");
  GeometricModelPtr gm = new GeometricModel(2, region_list, &comm);

  FrameworkPreference pref;
  pref.clear();
  pref.push_back(MSTK);
  pref.push_back(STKMESH);

  MeshFactory meshfactory(&comm);
  meshfactory.preference(pref);
  //RCP<const Mesh> mesh = meshfactory(low_coor[0], low_coor[1],
  //                                   high_coor[0], high_coor[1],
  //                                   mesh_size[0], mesh_size[1], gm);
  RCP<const Mesh> mesh = meshfactory(file_name, gm);
  int ncells = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  int nfaces_wghost = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::ALL);

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

  // write initial state
  for (std::vector<Teuchos::RCP<Visualization> >::iterator vis=visualization_.begin();
       vis!=visualization_.end(); ++vis) {
    WriteVis((*vis).ptr(), S.ptr());
  } 

  // Init the PK and solve for pressure
  MPK->Initialize();
  

  bool failed = true;
  double t_sim = 0.0;
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

  
  // get the data for viz
  Epetra_MultiVector& perm = *S->GetFieldData("permeability",passwd)->ViewComponent("cell");
  Epetra_MultiVector& p1 = *S->GetFieldData("phase1_pressure",passwd)->ViewComponent("cell");
  Epetra_MultiVector& p2 = *S->GetFieldData("phase2_pressure",passwd)->ViewComponent("cell");
  Epetra_MultiVector& s1 = *S->GetFieldData("phase1_saturation",passwd)->ViewComponent("cell");
  Epetra_MultiVector& s2 = *S->GetFieldData("phase2_saturation",passwd)->ViewComponent("cell");
  if (MyPID == 0) {
    GMV::open_data_file(*mesh, (std::string)"test_mpcoupled_q5.gmv");
    GMV::start_data();
    GMV::write_cell_data(perm, 0, "permeability");
    GMV::write_cell_data(p1, 0, "p1");
    GMV::write_cell_data(s1, 0, "s1");
    GMV::write_cell_data(p2, 0, "p2");
    GMV::write_cell_data(s2, 0, "s2");
    GMV::close_data_file();
  }

  Amanzi::VerboseObject *vo = new Amanzi::VerboseObject("MPCoupled_PK::", *flow_list); 

  if (vo->os_OK(Teuchos::VERB_MEDIUM)) {
    double global_ncells(0.0);
    double local_ncells(0.0);
    for (State::mesh_iterator mesh = S->mesh_begin(); mesh != S->mesh_end(); ++mesh) {
      Epetra_Map cell_map = (mesh->second.first)->cell_map(false);
      global_ncells += cell_map.NumGlobalElements();
      local_ncells += cell_map.NumMyElements();
    }    

    double mem = rss_usage();
    
    double percell(mem);
    if (local_ncells > 0) {
      percell = mem/local_ncells;
    }

    double max_percell(0.0);
    double min_percell(0.0);
    comm.MinAll(&percell,&min_percell,1);
    comm.MaxAll(&percell,&max_percell,1);

    double total_mem(0.0);
    double max_mem(0.0);
    double min_mem(0.0);
    comm.SumAll(&mem,&total_mem,1);
    comm.MinAll(&mem,&min_mem,1);
    comm.MaxAll(&mem,&max_mem,1);

    Teuchos::OSTab tab = vo->getOSTab();
    *vo->os() << "======================================================================" << std::endl;
    *vo->os() << "Simulation made " << S->cycle() << " cycles.\n";
    *vo->os() << "All meshes combined have " << global_ncells << " cells.\n";
    *vo->os() << "Memory usage (high water mark):\n";
    *vo->os() << std::fixed << std::setprecision(1);
    *vo->os() << "  Maximum per core:   " << std::setw(7) << max_mem 
               << " MBytes,  maximum per cell: " << std::setw(7) << max_percell*1024*1024 
               << " Bytes" << std::endl;
    *vo->os() << "  Minumum per core:   " << std::setw(7) << min_mem 
               << " MBytes,  minimum per cell: " << std::setw(7) << min_percell*1024*1024 
               << " Bytes" << std::endl;
    *vo->os() << "  Total:              " << std::setw(7) << total_mem 
               << " MBytes,  total per cell:   " << std::setw(7) << total_mem/global_ncells*1024*1024 
               << " Bytes" << std::endl;
  }

  double doubles_count(0.0);
  for (State::field_iterator field = S->field_begin(); field != S->field_end(); ++field) {
    doubles_count += static_cast<double>(field->second->GetLocalElementCount());
  }
  double global_doubles_count(0.0);
  double min_doubles_count(0.0);
  double max_doubles_count(0.0);
  comm.SumAll(&doubles_count,&global_doubles_count,1);
  comm.MinAll(&doubles_count,&min_doubles_count,1);
  comm.MaxAll(&doubles_count,&max_doubles_count,1);

  Teuchos::OSTab tab = vo->getOSTab();
  *vo->os() << "Doubles allocated in state fields " << std::endl;
  *vo->os() << "  Maximum per core:   " << std::setw(7)
             << max_doubles_count*8/1024/1024 << " MBytes" << std::endl;
  *vo->os() << "  Minimum per core:   " << std::setw(7)
             << min_doubles_count*8/1024/1024 << " MBytes" << std::endl; 
  *vo->os() << "  Total:              " << std::setw(7)
             << global_doubles_count*8/1024/1024 << " MBytes" << std::endl;

  delete MPK;
}
