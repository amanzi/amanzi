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
#include "Reduced2p2c_PK.hh"
#include "OperatorDefs.hh"
#include "TimeStepManager.hh"

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
TEST(MULTIPHASE_REDUCED_2P2C) {
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Flow;

  Comm_ptr_type comm = Amanzi::getDefaultComm();
  int MyPID = comm->MyPID();

  if (MyPID == 0) std::cout << "Test: multiphase component pk, uniform rectangular mesh" << std::endl;

  /* read parameter list */
  std::string xmlFileName = "test/test_reduced_2p2c.xml";
  ParameterXMLFileReader xmlreader(xmlFileName);
  ParameterList plist = xmlreader.getParameters();
  Teuchos::RCP<Teuchos::ParameterList> parameter_list = 
     Teuchos::rcp(new Teuchos::ParameterList(plist));
  Teuchos::ParameterList pk_tree_list = plist.sublist("PK Tree");
  Teuchos::RCP<Teuchos::ParameterList> state_list = Teuchos::sublist(parameter_list, "state");
  Teuchos::RCP<Teuchos::ParameterList> time_list = Teuchos::sublist(parameter_list, "cycle driver");
  Teuchos::RCP<Teuchos::ParameterList> MPMC_specs = Teuchos::sublist(parameter_list, "MPMC Specs");

  /* create a MSTK mesh framework */
  ParameterList region_list = plist.get<Teuchos::ParameterList>("regions");
  auto gm = Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(2, region_list, *comm));

  MeshFactory meshfactory(comm, gm);
  meshfactory.set_preference(Preference({Framework::MSTK, Framework::STK}));

  RCP<const Mesh> mesh = meshfactory.create(0.0, 0.0, 200.0, 20.0, 200, 10);
  int ncells = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  int nfaces_wghost = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::ALL);

  /* create a simple state populate it */
  RCP<State> S = rcp(new State(*state_list));
  //RCP<State> S = rcp(new State());
  S->RegisterDomainMesh(rcp_const_cast<Mesh>(mesh));
  std::string passwd("state");

  std::vector<Teuchos::RCP<Visualization> > visualization_;
  for (State::mesh_iterator mesh=S->mesh_begin(); mesh!=S->mesh_end(); ++mesh) {
    if (parameter_list->isSublist("visualization data")) {
      Teuchos::ParameterList& vis_plist = parameter_list->sublist("visualization data");
      Teuchos::RCP<Visualization> vis = Teuchos::rcp(new Visualization(vis_plist));
      vis->set_mesh(mesh->second.first);
      vis->CreateFiles();
      visualization_.push_back(vis);
    }
  }

  //std::vector<double> t_, tp_start_, tp_end_, tp_dt_, tp_max_cycle_, tp_max_dt_;  
  double max_dt_, min_dt_;

  max_dt_ = time_list->get<double>("max time step size", 864000);
  min_dt_ = time_list->get<double>("min time step size", 5400);
  double dT = time_list->get<double>("time step");
  double t_final = time_list->get<double>("end time");

  const Teuchos::RCP<CompositeVectorSpace> cvs = Teuchos::rcp(new CompositeVectorSpace());
  cvs->SetMesh(mesh);
  cvs->SetGhosted(true);
  cvs->SetComponent("cell", AmanziMesh::CELL, 1);
  cvs->SetOwned(false);

  Teuchos::RCP<CompositeVector> saturation = Teuchos::rcp(new CompositeVector(*cvs));
  Teuchos::RCP<CompositeVector> pressure = Teuchos::rcp(new CompositeVector(*cvs));
  Teuchos::RCP<CompositeVector> rhl = Teuchos::rcp(new CompositeVector(*cvs));
  saturation->PutScalar(1.0);
  pressure->PutScalar(1.0e6);
  rhl->PutScalar(0.0);

  Teuchos::RCP<CompositeVector> saturation_old = Teuchos::rcp(new CompositeVector(*cvs));
  Teuchos::RCP<CompositeVector> pressure_old = Teuchos::rcp(new CompositeVector(*cvs));
  Teuchos::RCP<CompositeVector> rhl_old = Teuchos::rcp(new CompositeVector(*cvs));

  Teuchos::RCP<TreeVector> sat_tree = Teuchos::rcp(new TreeVector());
  Teuchos::RCP<TreeVector> pres_tree = Teuchos::rcp(new TreeVector());
  Teuchos::RCP<TreeVector> rhl_tree = Teuchos::rcp(new TreeVector());
  sat_tree->SetData(saturation);
  pres_tree->SetData(pressure);
  rhl_tree->SetData(rhl);

  Teuchos::RCP<TreeVector> sat_tree_old = Teuchos::rcp(new TreeVector());
  Teuchos::RCP<TreeVector> pres_tree_old = Teuchos::rcp(new TreeVector());
  Teuchos::RCP<TreeVector> rhl_tree_old = Teuchos::rcp(new TreeVector());
  sat_tree_old->SetData(saturation_old);
  pres_tree_old->SetData(pressure_old);
  rhl_tree_old->SetData(rhl_old);

  Teuchos::RCP<CompositeVector> res1 = Teuchos::rcp(new CompositeVector(*cvs));
  Teuchos::RCP<CompositeVector> res2 = Teuchos::rcp(new CompositeVector(*cvs));
  Teuchos::RCP<CompositeVector> res3 = Teuchos::rcp(new CompositeVector(*cvs));
  Teuchos::RCP<TreeVector> f = Teuchos::rcp(new TreeVector());
  Teuchos::RCP<TreeVector> f1 = Teuchos::rcp(new TreeVector());
  Teuchos::RCP<TreeVector> f2 = Teuchos::rcp(new TreeVector());
  Teuchos::RCP<TreeVector> f3 = Teuchos::rcp(new TreeVector());
  f1->SetData(res1);
  f2->SetData(res2);
  f3->SetData(res3);
  f->PushBack(f1);
  f->PushBack(f2);
  f->PushBack(f3); 

  // create the global solution vector
  Teuchos::RCP<TreeVector> soln = Teuchos::rcp(new TreeVector());
  Teuchos::RCP<TreeVector> u_new = Teuchos::rcp(new TreeVector());
  u_new->PushBack(pres_tree);
  u_new->PushBack(sat_tree);
  u_new->PushBack(rhl_tree);
  Teuchos::RCP<TreeVector> u_old = Teuchos::rcp(new TreeVector());
  u_old->PushBack(pres_tree_old);
  u_old->PushBack(sat_tree_old);
  u_old->PushBack(rhl_tree_old);

  // Create pressure PK
  Multiphase::Reduced2p2c_PK* MPMC = new Multiphase::Reduced2p2c_PK(pk_tree_list, parameter_list, S, soln);
  S->Setup();
  S->InitializeFields();

  // Init the PK and solve for pressure
  MPMC->Initialize(S.ptr());

  // write initial state
  for (std::vector<Teuchos::RCP<Visualization> >::iterator vis=visualization_.begin();
       vis!=visualization_.end(); ++vis) {
    WriteVis((*vis).ptr(), S.ptr());
  }

  bool failed = true;
  double t_sim = 0.0;
  //double dT = 5000.0*365.0*24.0*3600.0;
  //double t_final = 10*dT;
  while (!(t_sim > t_final || abs(t_sim - t_final) < 1e-8)) {
    if (t_sim + dT - t_final > 1e-12) dT = t_final - t_sim;
    failed = MPMC->AdvanceStep(t_sim, t_sim + dT, false);
    if (failed) {
      dT = MPMC->get_dt();
    } else {
      t_sim += dT;
      MPMC->CommitStep(t_sim, t_sim + dT, S); 
      S->advance_time(dT);
      S->advance_cycle();
      if (MyPID == 0) {
        std::cout << "State time: " << S->time() << "; State cycle: " << S->cycle() << endl;
      }
      for (std::vector<Teuchos::RCP<Visualization> >::iterator vis=visualization_.begin();
           vis!=visualization_.end(); ++vis) {
        WriteVis((*vis).ptr(), S.ptr());
        //std::cout << "writing visualization file" << std::endl;
      } 
      dT = MPMC->get_dt();
    } 
  //std::cout << "Time step: " << std::round(t_sim/dT) << "\n"; 
  }

  //MPMC->Functional(0.0, 0.1, u_new, u_new, f);
  //MPMC->UpdatePreconditioner(0.0, u_new, 0.1);

  // get the data for viz
  Epetra_MultiVector& pw = *S->GetFieldData("pressure_w", passwd)->ViewComponent("cell", false);
  Epetra_MultiVector& pn = *S->GetFieldData("pressure_n", passwd)->ViewComponent("cell", false);
  Epetra_MultiVector& sw = *S->GetFieldData("saturation_w", passwd)->ViewComponent("cell", false);
  Epetra_MultiVector& sn = *S->GetFieldData("saturation_n", passwd)->ViewComponent("cell", false);
  Epetra_MultiVector& hydrogen_density = *soln->SubVector(2)->Data()->ViewComponent("cell", false);
  if (MyPID == 0) {
    GMV::open_data_file(*mesh, (std::string)"test_reduced_2p2c.gmv");
    GMV::start_data();
    GMV::write_cell_data(pw, 0, "pressure_w");
    GMV::write_cell_data(pn, 0, "pressure_n");
    GMV::write_cell_data(sw, 0, "saturation_w");
    GMV::write_cell_data(sn, 0, "saturation_n");
    GMV::write_cell_data(hydrogen_density, 0, "hydrogen_density");
    GMV::close_data_file();
  }
  MPI_Comm mpi_comm(MPI_COMM_WORLD);
  Amanzi::timer_manager.parSync(mpi_comm);
  Amanzi::timer_manager.print();

  Amanzi::VerboseObject *vo = new Amanzi::VerboseObject("Reduced2p2c::", *MPMC_specs);

  Teuchos::OSTab tab = vo->getOSTab();
  *vo->os() << "======================================================================" << std::endl;
  *vo->os() << "Accumulate setup time: " << MPMC->accumulateSetupTime_ << std::endl;
  *vo->os() << "Accumulate solve time: " << MPMC->accumulateSolveTime_ << std::endl;

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
    comm->MinAll(&percell,&min_percell,1);
    comm->MaxAll(&percell,&max_percell,1);

    double total_mem(0.0);
    double max_mem(0.0);
    double min_mem(0.0);
    comm->SumAll(&mem,&total_mem,1);
    comm->MinAll(&mem,&min_mem,1);
    comm->MaxAll(&mem,&max_mem,1);

    //Teuchos::OSTab tab = vo->getOSTab();
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
  comm->SumAll(&doubles_count,&global_doubles_count,1);
  comm->MinAll(&doubles_count,&min_doubles_count,1);
  comm->MaxAll(&doubles_count,&max_doubles_count,1);

  //Teuchos::OSTab tab = vo->getOSTab();
  *vo->os() << "Doubles allocated in state fields " << std::endl;
  *vo->os() << "  Maximum per core:   " << std::setw(7)
             << max_doubles_count*8/1024/1024 << " MBytes" << std::endl;
  *vo->os() << "  Minimum per core:   " << std::setw(7)
             << min_doubles_count*8/1024/1024 << " MBytes" << std::endl; 
  *vo->os() << "  Total:              " << std::setw(7)
             << global_doubles_count*8/1024/1024 << " MBytes" << std::endl;

  delete MPMC;
}
