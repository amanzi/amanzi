/*
  Multi-Process Coordinator

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon
          Daniil Svyatskiy

  Implementation for the CycleDriver.  CycleDriver is basically just a class to hold
  the cycle driver, which runs the overall, top level timestep loop.  It
  instantiates states, ensures they are initialized, and runs the timestep loop
  including Vis and restart/checkpoint dumps.  It contains one and only one PK
  -- most likely this PK is an MPC of some type -- to do the actual work.
*/

#include <ios>
#include <iostream>
#include <unistd.h>
#include <sys/resource.h>

// TPLs
#define BOOST_FILESYSTEM_NO_DEPRECATED
#include "boost/filesystem/operations.hpp"

#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

// Amanzi
#include "Checkpoint.hh"
#include "errors.hh"
#include "IO.hh"
#include "PK.hh"
#include "PK_Factory.hh"
#include "ObservationData.hh"
#include "TimeStepManager.hh"
#include "TreeVector.hh"
#include "State.hh"
#include "Visualization.hh"
#include "MeshInfo.hh"

// MPC
#include "CycleDriver.hh"
#include "FlexibleObservations.hh"

#define DEBUG_MODE 1

namespace Amanzi {

bool reset_info_compfunc(const std::pair<double,double>& x, const std::pair<double,double>& y) {
  return (x.first < y.first);
}


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


/* ******************************************************************
* Constructor.
****************************************************************** */
CycleDriver::CycleDriver(Teuchos::RCP<Teuchos::ParameterList> glist,
                         Teuchos::RCP<Amanzi::State>& S,
                         const Comm_ptr_type& comm,
                         Amanzi::ObservationData& observations_data) :
    S_(S),
    glist_(glist),
    comm_(comm),
    observations_data_(observations_data),
    restart_requested_(false) {

  // create and start the global timer
  CoordinatorInit_();

  vo_ = Teuchos::rcp(new VerboseObject("CycleDriver", glist_->sublist("cycle driver")));
}


/* ******************************************************************
* High-level initialization.
****************************************************************** */
void CycleDriver::CoordinatorInit_() {
  coordinator_list_ = Teuchos::sublist(glist_, "cycle driver");
  ReadParameterList_();

  // create the global solution vector
  soln_ = Teuchos::rcp(new TreeVector());
}


/* ******************************************************************
* Create the pk tree root node which then creates the rest of the tree.
****************************************************************** */
void CycleDriver::Init_PK(int time_pr_id) {
  PKFactory pk_factory;

  Teuchos::RCP<Teuchos::ParameterList> time_periods_list = Teuchos::sublist(coordinator_list_, "time periods", true);

  std::ostringstream ss;
  ss << time_pr_id;

  std::string tp_list_name = "TP " + ss.str();
  Teuchos::RCP<Teuchos::ParameterList> tp_list =Teuchos::sublist( time_periods_list, tp_list_name.data(), true);
  Teuchos::ParameterList pk_tree_list = tp_list->sublist("PK tree");

  if (pk_tree_list.numParams() != 1) {
    Errors::Message message("CycleDriver: PK tree list should contain exactly one root node list");
    Exceptions::amanzi_throw(message);
  }
  Teuchos::ParameterList::ConstIterator pk_item = pk_tree_list.begin();
  const std::string &pk_name = pk_tree_list.name(pk_item);

  if (!pk_tree_list.isSublist(pk_name)) {
    Errors::Message message("CycleDriver: PK tree list does not have node \"" + pk_name + "\".");
    Exceptions::amanzi_throw(message);
  }

  pk_ = pk_factory.CreatePK(pk_name, pk_tree_list, glist_, S_, soln_);
}


/* ******************************************************************
* Setup PK first follwed by State's setup.
****************************************************************** */
void CycleDriver::Setup() {
  // Set up the states, creating all data structures.

  // create the observations
  if (glist_->isSublist("observation data")) {
    Teuchos::RCP<Teuchos::ParameterList> obs_list = Teuchos::sublist(glist_, "observation data");
    Teuchos::RCP<Teuchos::ParameterList> units_list = Teuchos::sublist(glist_, "units");
    observations_ = Teuchos::rcp(new Amanzi::FlexibleObservations(
        coordinator_list_, obs_list, units_list, observations_data_, S_));

    if (coordinator_list_->isParameter("component names")) {
      std::vector<std::string> names = coordinator_list_->get<Teuchos::Array<std::string> >("component names").toVector();
      std::vector<double> mol_masses = coordinator_list_->get<Teuchos::Array<double> >("component molar masses").toVector();
      int num_liquid = coordinator_list_->get<int>("number of liquid components", names.size());

      observations_->RegisterComponentNames(names, mol_masses, num_liquid);
    }
  }

  // create the checkpointing
  if (glist_->isSublist("checkpoint data")) {
    Teuchos::ParameterList& chkp_plist = glist_->sublist("checkpoint data");
    checkpoint_ = Teuchos::rcp(new Amanzi::Checkpoint(chkp_plist, *S_));
  }
  else {
    checkpoint_ = Teuchos::rcp(new Amanzi::Checkpoint());
  }

  // create the walkabout
  if (glist_->isSublist("walkabout data")) {
    Teuchos::ParameterList& walk_plist = glist_->sublist("walkabout data");
    walkabout_ = Teuchos::rcp(new Amanzi::WalkaboutCheckpoint(walk_plist, *S_));
  }
  else {
    walkabout_ = Teuchos::rcp(new Amanzi::WalkaboutCheckpoint());
  }

  // vis successful steps
  bool surface_done = false;
  for (auto mesh = S_->mesh_begin(); mesh != S_->mesh_end(); ++mesh) {
    if (mesh->first == "surface_3d") {
      // pass
    } else if ((mesh->first == "surface") && surface_done) {
      // pass
    } else {
      // vis successful steps
      std::string plist_name = "visualization data " + mesh->first;
      // in the case of just a domain mesh, we want to allow no name.
      if ((mesh->first == "domain") && !glist_->isSublist(plist_name)) {
        plist_name = "visualization data";
      }

      if (glist_->isSublist(plist_name)) {
        Teuchos::ParameterList& vis_plist = glist_->sublist(plist_name);
        Teuchos::RCP<Visualization> vis = Teuchos::rcp(new Visualization(vis_plist));
        vis->set_mesh(mesh->second.first);
        vis->set_name(mesh->first);
        vis->CreateFiles();
        visualization_.push_back(vis);
      }

      // vis unsuccessful steps
      std::string fail_plist_name = "visualization data " + mesh->first + " Failed Steps";
      // in the case of just a domain mesh, we want to allow no name.
      if ((mesh->first == "domain") && !glist_->isSublist(fail_plist_name)) {
        fail_plist_name = "visualization data failed steps";
      }

      if (glist_->isSublist(fail_plist_name)) {
        Teuchos::ParameterList& fail_vis_plist = glist_->sublist(fail_plist_name);
        Teuchos::RCP<Visualization> fail_vis =
          Teuchos::rcp(new Visualization(fail_vis_plist));
        fail_vis->set_mesh(mesh->second.first);
        fail_vis->CreateFiles();
        failed_visualization_.push_back(fail_vis);
      }
    }
    std::string plist_name = "mesh info " + mesh->first;
    // in the case of just a domain mesh, we want to allow no name.
    if ((mesh->first == "domain") && !glist_->isSublist(plist_name)) {
      plist_name = "mesh info";
    }
    if (glist_->isSublist(plist_name)) {
      auto& mesh_info_list = glist_->sublist(plist_name);
      Teuchos::RCP<Amanzi::MeshInfo> mesh_info = Teuchos::rcp(new Amanzi::MeshInfo(mesh_info_list, *S_));
      mesh_info->WriteMeshCentroids(mesh->first, *(mesh->second.first));
    }
  }

  pk_->Setup();
  pk_->set_tags(Tags::CURRENT, Tags::NEXT);
  S_->Require<double>("dt", Tags::NEXT, "dt");
  S_->Setup();

  // create the time step manager
  tsm_ = Teuchos::rcp(new TimeStepManager(glist_->sublist("cycle driver")));

  // set up the TSM
  // -- register visualization times
  for (std::vector<Teuchos::RCP<Visualization> >::iterator vis=visualization_.begin();
       vis!=visualization_.end(); ++vis) {
    (*vis)->RegisterWithTimeStepManager(tsm_.ptr());
  }
  // -- register checkpoint times
  if (checkpoint_ != Teuchos::null) 
  checkpoint_->RegisterWithTimeStepManager(tsm_.ptr());
  // -- register observation times
  if (observations_ != Teuchos::null) 
    observations_->RegisterWithTimeStepManager(tsm_.ptr());
  // -- register the final time
  // register reset_times
  for(std::vector<std::pair<double,double> >::const_iterator it = reset_info_.begin();
      it != reset_info_.end(); ++it) tsm_->RegisterTimeEvent(it->first);


  for (int i=0; i<num_time_periods_; i++) {
    tsm_->RegisterTimeEvent(tp_end_[i]);
    tsm_->RegisterTimeEvent(tp_start_[i] + tp_dt_[i]);
  } 


  if (vo_->os_OK(Teuchos::VERB_MEDIUM)) {
    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "Setup is complete." << std::endl;
  }
}


/* ******************************************************************
* Initialize State followed by initialization of PK.
****************************************************************** */
void CycleDriver::Initialize() {
 
  S_->GetW<double>("dt", Tags::DEFAULT, "dt") = tp_dt_[0];
  S_->GetRecordW("dt", "dt").set_initialized();

  // Initialize the state (initializes all dependent variables).
  S_->InitializeFields();
  S_->InitializeEvaluators();

  // Initialize the process kernels and verify
  pk_->Initialize();
  S_->CheckAllFieldsInitialized();

  // S_->WriteDependencyGraph();
  S_->InitializeIOFlags(); 

  // commit the initial conditions.
  if (!restart_requested_) {
    pk_->CommitStep(S_->get_time(), S_->get_time(), Tags::DEFAULT);
  }
}


/* ******************************************************************
* Force checkpoint at the end of simulation.
* Only do if the checkpoint was not already written, or we would be writing
* the same file twice.
* This really should be removed, but for now is left to help stupid developers.
****************************************************************** */
void CycleDriver::Finalize() {
  if (!checkpoint_->DumpRequested(S_->get_cycle(), S_->get_time())) {
    pk_->CalculateDiagnostics(Tags::DEFAULT);
    checkpoint_->Write(*S_, 0.0, true, &observations_data_);
  }
}


/* ******************************************************************
* Report the memory high water mark (using ru_maxrss)
* this should be called at the very end of a simulation
****************************************************************** */
void CycleDriver::ReportMemory() {
  if (vo_->os_OK(Teuchos::VERB_MEDIUM)) {
    double global_ncells(0.0);
    double local_ncells(0.0);
    for (State::mesh_iterator mesh = S_->mesh_begin(); mesh != S_->mesh_end(); ++mesh) {
      Epetra_Map cell_map = (mesh->second.first)->cell_map(false);
      global_ncells += cell_map.NumGlobalElements();
      local_ncells += cell_map.NumMyElements();
    }    

    double mem = Amanzi::rss_usage();
    
    double percell(mem);
    if (local_ncells > 0) {
      percell = mem/local_ncells;
    }

    double max_percell(0.0);
    double min_percell(0.0);
    comm_->MinAll(&percell,&min_percell,1);
    comm_->MaxAll(&percell,&max_percell,1);

    double total_mem(0.0);
    double max_mem(0.0);
    double min_mem(0.0);
    comm_->SumAll(&mem,&total_mem,1);
    comm_->MinAll(&mem,&min_mem,1);
    comm_->MaxAll(&mem,&max_mem,1);

    Teuchos::OSTab tab = vo_->getOSTab();
    std::ios save(NULL);
    save.copyfmt(*vo_->os());
    *vo_->os() << "======================================================================" << std::endl;
    *vo_->os() << "Simulation made " << S_->get_cycle() << " cycles.\n";
    *vo_->os() << "All meshes combined have " << global_ncells << " cells.\n";
    *vo_->os() << "Memory usage (high water mark):\n";
    *vo_->os() << std::fixed << std::setprecision(1);
    *vo_->os() << "  Maximum per core:   " << std::setw(7) << max_mem 
               << " MBytes,  maximum per cell: " << std::setw(7) << max_percell*1024*1024 
               << " Bytes" << std::endl;
    *vo_->os() << "  Minumum per core:   " << std::setw(7) << min_mem 
               << " MBytes,  minimum per cell: " << std::setw(7) << min_percell*1024*1024 
               << " Bytes" << std::endl;
    *vo_->os() << "  Total:              " << std::setw(7) << total_mem 
               << " MBytes,  total per cell:   " << std::setw(7) << total_mem/global_ncells*1024*1024 
               << " Bytes" << std::endl;
    vo_->os()->copyfmt(save);
  }

  double doubles_count(0.0);
  for (auto it = S_->data_begin(); it != S_->data_end(); ++it) {
    if (S_->GetRecord(it->first).ValidType<CompositeVector>()) {
      doubles_count += static_cast<double>(S_->Get<CompositeVector>(it->first).GetLocalElementCount());
    }
  }
  double global_doubles_count(0.0);
  double min_doubles_count(0.0);
  double max_doubles_count(0.0);
  comm_->SumAll(&doubles_count,&global_doubles_count,1);
  comm_->MinAll(&doubles_count,&min_doubles_count,1);
  comm_->MaxAll(&doubles_count,&max_doubles_count,1);

  if (vo_->os_OK(Teuchos::VERB_LOW)) {
    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "Doubles allocated in state fields " << std::endl;
    *vo_->os() << "  Maximum per core:   " << std::setw(7)
               << max_doubles_count*8/1024/1024 << " MBytes" << std::endl;
    *vo_->os() << "  Minimum per core:   " << std::setw(7)
               << min_doubles_count*8/1024/1024 << " MBytes" << std::endl; 
    *vo_->os() << "  Total:              " << std::setw(7)
               << global_doubles_count*8/1024/1024 << " MBytes" << std::endl;
  }
}


/* ******************************************************************
* Read control variables for all time periods. 
****************************************************************** */
void CycleDriver::ReadParameterList_() {
  max_dt_ = coordinator_list_->get<double>("max time step size", 1.0e+99);
  min_dt_ = coordinator_list_->get<double>("min time step size", 1.0e-12);
  cycle0_ = coordinator_list_->get<int>("start cycle", 0);
  cycle1_ = coordinator_list_->get<int>("end cycle", -1);

  Teuchos::ParameterList time_periods_list = coordinator_list_->sublist("time periods");

  num_time_periods_ = time_periods_list.numParams();
  Teuchos::ParameterList::ConstIterator item;
  tp_start_.resize(num_time_periods_);
  tp_end_.resize(num_time_periods_);
  tp_dt_.resize(num_time_periods_);
  tp_max_cycle_.resize(num_time_periods_);
  tp_max_dt_.resize(num_time_periods_);

  int i = 0;
  for (item = time_periods_list.begin(); item !=time_periods_list.end(); ++item) {
    const std::string & tp_name = time_periods_list.name(item);
    tp_start_[i] = time_periods_list.sublist(tp_name).get<double>("start period time");
    tp_end_[i] = time_periods_list.sublist(tp_name).get<double>("end period time");
    tp_dt_[i] = time_periods_list.sublist(tp_name).get<double>("initial time step", 1.0);
    tp_max_dt_[i] = time_periods_list.sublist(tp_name).get<double>("maximum time step", 1.0e+99);
    tp_max_cycle_[i] = time_periods_list.sublist(tp_name).get<int>("maximum cycle number", -1);
    i++;
  }

  // restart control
  // are we restarting from a file?
  // first assume we're not
  restart_requested_ = false;

  if (coordinator_list_->isSublist("restart")) {
    restart_requested_ = true;

    Teuchos::ParameterList restart_list = coordinator_list_->sublist("restart");
    restart_filename_ = restart_list.get<std::string>("file name");

    // make sure that the restart file actually exists, if not throw an error
    boost::filesystem::path restart_from_filename_path(restart_filename_);
    if (!boost::filesystem::exists(restart_from_filename_path)) {
      Errors::Message msg;
      msg << "CycleDriver: restart file \"" << restart_filename_ 
          << "\" does not exist or is not a regular file.";
      Exceptions::amanzi_throw(msg);
    }
  }

  if (coordinator_list_->isSublist("time period control")) {
    Teuchos::ParameterList& tpc_list = coordinator_list_->sublist("time period control");
    Teuchos::Array<double> reset_times = tpc_list.get<Teuchos::Array<double> >("start times");
    Teuchos::Array<double> reset_times_dt = tpc_list.get<Teuchos::Array<double> >("initial time step");   
    AMANZI_ASSERT(reset_times.size() == reset_times_dt.size());

    {
      Teuchos::Array<double>::const_iterator it_tim;
      Teuchos::Array<double>::const_iterator it_dt;
      for (it_tim = reset_times.begin(), it_dt = reset_times_dt.begin();
           it_tim != reset_times.end();
         ++it_tim, ++it_dt) {
        reset_info_.push_back(std::make_pair(*it_tim, *it_dt));
      }
    }  

    if (tpc_list.isParameter("maximum time step")) {
      Teuchos::Array<double> reset_max_dt = tpc_list.get<Teuchos::Array<double> >("maximum time step");
      AMANZI_ASSERT(reset_times.size() == reset_max_dt.size());

      Teuchos::Array<double>::const_iterator it_tim;
      Teuchos::Array<double>::const_iterator it_max;
      for (it_tim = reset_times.begin(), it_max = reset_max_dt.begin();
           it_tim != reset_times.end();
           ++it_tim, ++it_max) {
        reset_max_.push_back(std::make_pair(*it_tim, *it_max));
      }  
    }

    // now we sort in ascending order by time
    std::sort(reset_info_.begin(), reset_info_.end(), reset_info_compfunc);
    std::sort(reset_max_.begin(),  reset_max_.end(),  reset_info_compfunc);
  }

  // verification (move this to state ?)
  S_->GetMeshPartition("materials");
}


/* ******************************************************************
* Acquire the chosen timestep size
*******************************************************************/
double CycleDriver::get_dt(bool after_failure) {
  // get the physical step size
  double dt;
 
  dt = pk_->get_dt();

  std::vector<std::pair<double,double> >::iterator it;
  std::vector<std::pair<double,double> >::iterator it_max;

  for (it = reset_info_.begin(), it_max = reset_max_.begin(); it != reset_info_.end(); ++it, ++it_max) {    
    if (S_->get_time() == it->first) {
      if (reset_max_.size() > 0) {
        max_dt_ = it_max->second;
      }

      if (dt < it->second) {
        pk_->set_dt(dt);
      } else {
        dt = it->second;
        pk_->set_dt(dt);
        after_failure = true;
      }
      
      it = reset_info_.erase(it);
      if (reset_max_.size() > 0) 
        it_max = reset_max_.erase(it_max);

      break;
    }
  }

  // check if the step size has gotten too small
  if (dt < min_dt_) {
    Errors::Message message("CycleDriver: error, timestep too small");
    Exceptions::amanzi_throw(message);
  }

  if (S_->get_time() > 0) {
    if (dt / S_->get_time() < 1e-14) {
      Errors::Message message("CycleDriver: error, timestep too small with respect to current time");
      Exceptions::amanzi_throw(message);
    }
  }


  // ask the step manager if this step is ok
  dt = tsm_->TimeStep(S_->get_time(), dt, after_failure);


  // cap the max step size
  if (dt > max_dt_) {
    dt = max_dt_;   
    Utils::Units units("molar");
    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "Resetting time step to the maximum allowed of " << units.OutputTime(dt) << "\n";
  }

  return dt;
}


/* ******************************************************************
* Time step management.
****************************************************************** */
void CycleDriver::set_dt(double dt) {
  double dt_;

  // check if the step size has gotten too small
  if (dt < min_dt_) {
    Errors::Message message("CycleDriver: error, timestep too small");
    Exceptions::amanzi_throw(message);
  }

  // cap the max step size
  if (dt > max_dt_) {
    dt_ = max_dt_;
  }

  // ask the step manager if this step is ok
  dt_ = tsm_->TimeStep(S_->get_time() + dt, dt);

  // set the physical step size
  pk_->set_dt(dt_);
}


/* ******************************************************************
* This is used by CLM
****************************************************************** */
double CycleDriver::Advance(double dt) {

  bool advance = true;
  bool fail = false;
  bool reinit = false;
  double dt_new;

  if (tp_end_[time_period_id_] == tp_start_[time_period_id_]) 
    advance = false;

  Teuchos::OSTab tab = vo_->getOSTab();
  
  if (advance) {
    std::vector<std::pair<double,double> >::const_iterator it;
    for (it = reset_info_.begin(); it != reset_info_.end(); ++it) {
      if (it->first == S_->get_time()) break;
    }

    if (it != reset_info_.end()) {
      if (vo_->os_OK(Teuchos::VERB_MEDIUM)) {
        *vo_->os() << vo_->color("blue") << " Reinitializing PKs due to BCs or sources/sinks" 
                   << vo_->reset() << std::endl;
      }
      reinit = true;
    }      

    fail = pk_->AdvanceStep(S_->get_time(), S_->get_time() + dt, reinit);
  }

  if (!fail) {
    pk_->CommitStep(S_->get_time(), S_->get_time() + dt, Tags::DEFAULT);
    // advance the iteration count and timestep size
    if (advance) {
      S_->advance_cycle();
      S_->advance_time(dt);
    }

    bool force_vis(false);
    bool force_check(false);
    bool force_obser(false);

    if (abs(S_->get_time() - tp_end_[time_period_id_]) < 1e-10) {
      force_vis = true;
      force_check = true;                       
      force_obser = true;
      S_->set_position(TIME_PERIOD_END);
    }

    dt_new = get_dt(fail);

    if (!reset_info_.empty())
        if (S_->get_time() == reset_info_.front().first)
            force_check = true;

    // make vis, observations, and checkpoints in this order
    // Amanzi::timer_manager.start("I/O");
    if (advance) {
      pk_->CalculateDiagnostics(Tags::DEFAULT);
      Visualize(force_vis);
      Observations(force_obser, true);
      WriteCheckpoint(dt_new, force_check);   // write Checkpoint with new dt
      WriteWalkabout(force_check);
    }
    // Amanzi::timer_manager.start("I/O");

    if (vo_->os_OK(Teuchos::VERB_MEDIUM)) {
      Utils::Units units("molar");
      *vo_->os() << "New time = " << units.OutputTime(S_->get_time()) << std::endl;
    }
  } else {
    dt_new = get_dt(fail);
    // Failed the timestep.  
    // Potentially write out failed timestep for debugging
    for (auto& vis : failed_visualization_) {
      WriteVis(*vis, *S_);
    }
    // The timestep sizes have been updated, so copy back old soln and try again.
    // NOT YET IMPLEMENTED, requires PKs to deal with failure.  Fortunately
    // transport and chemistry never fail, so we shouldn't break things.
    // Otherwise this would be very broken, as flow could succeed, but
    // transport fail, and we wouldn't have a way of backing up. --ETC
  }
  return dt_new;
}


/* ******************************************************************
* Make observations
****************************************************************** */
void CycleDriver::Observations(bool force, bool integrate) {
  if (observations_ != Teuchos::null) {
    // integrate continuous observations in time and save results in internal variables
    if (integrate) observations_->MakeContinuousObservations(*S_);

    if (observations_->DumpRequested(S_->get_cycle(), S_->get_time()) || force) {
      // continuous observations are not updated here
      int n = observations_->MakeObservations(*S_);
      Teuchos::OSTab tab = vo_->getOSTab();
      *vo_->os() << "Cycle " << S_->get_cycle() << ": writing observations... " << n << std::endl;
    }
  }
}


/* ******************************************************************
* Write visualization file with extension (cycle + tag) if requested.
****************************************************************** */
void CycleDriver::Visualize(bool force, const Tag& tag) {
  bool dump = force;
  if (!dump) {
    for (auto vis = visualization_.begin(); vis != visualization_.end(); ++vis) {
      if ((*vis)->DumpRequested(S_->get_cycle(), S_->get_time())) {
        dump = true;
      }
    }
  }

  // if (dump || force) pk_->CalculateDiagnostics();
  
  for (const auto& vis : visualization_) {
    if (force || vis->DumpRequested(S_->get_cycle(), S_->get_time())) {
      vis->set_tag(tag);
      WriteVis(*vis, *S_);
      Teuchos::OSTab tab = vo_->getOSTab();
      *vo_->os() << "writing visualization file: " << vis->get_name() << std::endl;
    }
  }
}


/* ******************************************************************
* Write a checkpoint file if requested.
****************************************************************** */
void CycleDriver::WriteCheckpoint(double dt, bool force) {
  if (force || checkpoint_->DumpRequested(S_->get_cycle(), S_->get_time())) {
    bool final = false;

    if (fabs(S_->get_time() - tp_end_[num_time_periods_-1]) < 1e-6) {
      final = true;
    }

    // checkpoint uses dt from State, but we do not want to modify state 
    // in this place and reset dt temporarily
    double dt_save = S_->Get<double>("dt", Tags::DEFAULT);
    S_->GetW<double>("dt", Tags::DEFAULT, "dt") = dt;

    checkpoint_->Write(*S_, dt, final, &observations_data_);

    S_->GetW<double>("dt", Tags::DEFAULT, "dt") = dt_save;
    
    if (vo_->os_OK(Teuchos::VERB_LOW)) {
      Teuchos::OSTab tab = vo_->getOSTab();
      *vo_->os() << "writing checkpoint file" << std::endl;
    }
  }
}


void CycleDriver::WriteWalkabout(bool force) {
  if (walkabout_ != Teuchos::null) {
    if (walkabout_->DumpRequested(S_->get_cycle(), S_->get_time()) || force) {
      if (!walkabout_->is_disabled())
         *vo_->os() << "Cycle " << S_->get_cycle() << ": writing walkabout file" << std::endl;
      walkabout_->WriteDataFile(S_, pk_);
    }
  }
}


/* ******************************************************************
* timestep loop.
****************************************************************** */
Teuchos::RCP<State> CycleDriver::Go() {

  time_period_id_ = 0;
  int position = 0;
  double restart_time = 0.;

  double dt;
  double restart_dT(1.0e99);

  if (!restart_requested_) {  // No restart

    Init_PK(time_period_id_);

    Setup();

    // start at time t = t0 and initialize data.
    S_->set_time(tp_start_[time_period_id_]);
    S_->set_cycle(cycle0_);
    S_->set_position(TIME_PERIOD_START);
    
    Initialize();

    dt = tp_dt_[time_period_id_];
    max_dt_ = tp_max_dt_[time_period_id_];
    dt = tsm_->TimeStep(S_->get_time(), dt);
    pk_->set_dt(dt);

  } else {

    // Read restart file
    restart_time = ReadCheckpointInitialTime(comm_, restart_filename_);

    position = ReadCheckpointPosition(comm_, restart_filename_);
    ReadCheckpointObservations(comm_, restart_filename_, observations_data_);
    for (int i = 0; i < num_time_periods_; i++) {
      if (restart_time - tp_end_[i] > -1e-10) 
        time_period_id_++;
    }  

    if (position == TIME_PERIOD_END) 
      if (time_period_id_>0) 
        time_period_id_--;   


    Init_PK(time_period_id_); 
    Setup();
    // Only field which are in State are initialize from the input file
    // to initialize field which are not in the restart file
    S_->InitializeFields();
    S_->InitializeEvaluators();
   
    S_->GetMeshPartition("materials");
    
    // re-initialize the state object
    ReadCheckpoint(comm_, *S_, restart_filename_);
    restart_dT = S_->Get<double>("dt");

    cycle0_ = S_->get_cycle();
    for (std::vector<std::pair<double,double> >::iterator it = reset_info_.begin();
          it != reset_info_.end(); ++it) {
      if (it->first <= S_->get_time()) it = reset_info_.erase(it);
      if (it == reset_info_.end() ) break;
    }

    if (vo_->os_OK(Teuchos::VERB_LOW)) {
      Teuchos::OSTab tab = vo_->getOSTab();
      *vo_->os() << "Restarted from checkpoint file: " << restart_filename_ << std::endl;
    }

    if (position == TIME_PERIOD_END && time_period_id_ < num_time_periods_ - 1) {
      time_period_id_++;
      ResetDriver(time_period_id_); 
      restart_dT = tp_dt_[time_period_id_];
    }
    else {
      pk_->Initialize();
    }     
    max_dt_ = tp_max_dt_[time_period_id_];

    S_->set_initial_time(S_->get_time());
    dt = tsm_->TimeStep(S_->get_time(), restart_dT);
    pk_->set_dt(dt);
  }

  // enfoce consistent physics after initialization
  // this is optional but helps with statistics
  S_->GetW<double>("dt", Tags::NEXT, "dt") = dt;
  S_->GetRecordW("dt", Tags::NEXT, "dt").set_initialized();

  S_->InitializeEvaluators();
  S_->CheckAllFieldsInitialized();

  // visualization at IC
  // Amanzi::timer_manager.start("I/O");
  // after initialization of State and PK we know all fields and clean of
  S_->InitializeIOFlags(); 

  pk_->CalculateDiagnostics(Tags::DEFAULT);
  Visualize();
  Observations();
  WriteCheckpoint(dt);
  WriteStateStatistics(*S_, *vo_);

  // Amanzi::timer_manager.stop("I/O");
  if (vo_->os_OK(Teuchos::VERB_MEDIUM)) {
    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "\nSimulation end time: " << tp_end_[time_period_id_] << " sec." << std::endl;
    *vo_->os() << "CPU time stamp: " << vo_->clock() << std::endl;
  }

  // iterate process kernels
  {
#if !DEBUG_MODE
  try {
#endif

    while (time_period_id_ < num_time_periods_) {
      int start_cycle_num = S_->get_cycle();
      //  do 
      while ((S_->get_time() < tp_end_[time_period_id_]) && 
             ((tp_max_cycle_[time_period_id_] == -1) || 
              (S_->get_cycle() - start_cycle_num < tp_max_cycle_[time_period_id_])))
      {
        if (vo_->os_OK(Teuchos::VERB_MEDIUM)) {
          if (S_->get_cycle() % 100 == 0 && S_->get_cycle() > 0) {
            WriteStateStatistics(*S_, *vo_);
            Teuchos::OSTab tab = vo_->getOSTab();
            *vo_->os() << "\nSimulation end time: " << tp_end_[time_period_id_] << " sec." << std::endl;
            *vo_->os() << "CPU time stamp: " << vo_->clock() << std::endl;
          }
          Utils::Units units("molar");
          Teuchos::OSTab tab = vo_->getOSTab();
          *vo_->os() << "\nCycle " << S_->get_cycle()
                     << ": time = " << units.OutputTime(S_->get_time())
                     << ", dt = " << units.OutputTime(dt) << "\n";
        }
        S_->GetW<double>("dt", Tags::DEFAULT, "dt") = dt;
        S_->set_initial_time(S_->get_time());
        S_->set_final_time(S_->get_time() + dt);
        S_->set_intermediate_time(S_->get_time());
        S_->set_position(TIME_PERIOD_INSIDE);

        dt = Advance(dt);
        // dt = get_dt(fail);
      }  // while not finished


      time_period_id_++;
      if (time_period_id_ < num_time_periods_) {
        WriteStateStatistics(*S_, *vo_);
        ResetDriver(time_period_id_); 
        dt = get_dt(false);
      }      
    }
  }
#if !DEBUG_MODE
  catch (Exceptions::Amanzi_exception &e) {
    // write one more vis for help debugging
    S_->advance_cycle();
    Visualize(true);  // force vis

    // flush observations to make sure they are saved
    observations_->Flush();

    // catch errors to dump two checkpoints -- one as a "last good" checkpoint
    // and one as a "debugging data" checkpoint.
    checkpoint_->set_filebasename("error_checkpoint");
    checkpoint_->Write(*S_, dt);
    throw e;
  }
#endif
  
  // finalizing simulation
  WriteStateStatistics(*S_, *vo_);
  ReportMemory();
  // Finalize();
 
  return S_;
} 


/* ******************************************************************
* TBW.
****************************************************************** */
void CycleDriver::ResetDriver(int time_pr_id) {

  if (vo_->os_OK(Teuchos::VERB_LOW)) {
    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "Reseting CD: TP " << time_pr_id - 1 << " -> TP " << time_pr_id << "." << std::endl;
  }

  Teuchos::RCP<AmanziMesh::Mesh> mesh =
    Teuchos::rcp_const_cast<AmanziMesh::Mesh>(S_->GetMesh("domain"));
 
  S_old_ = S_;

  Teuchos::ParameterList state_plist = glist_->sublist("state");
  S_ = Teuchos::rcp(new Amanzi::State(state_plist));

  for (auto it = S_old_->mesh_begin(); it != S_old_->mesh_end(); ++it) {
    S_->RegisterMesh(it->first, it->second.first);
  }    
  
  // delete the old global solution vector
  pk_ = Teuchos::null;
  soln_ = Teuchos::null;

  // create the global solution vector
  soln_ = Teuchos::rcp(new TreeVector());
  
  // create new pk
  Init_PK(time_pr_id);

  // register observation times with the time step manager
  //if (observations_ != Teuchos::null) observations_->RegisterWithTimeStepManager(tsm_);

  // Setup
  pk_->Setup();
  pk_->set_tags(Tags::CURRENT, Tags::NEXT);

  S_->Require<double>("dt", Tags::NEXT, "dt");
  S_->Setup();

  S_->set_cycle(S_old_->get_cycle());
  S_->set_time(tp_start_[time_pr_id]); 
  S_->set_position(TIME_PERIOD_START);

  S_->GetW<double>("dt", Tags::NEXT, "dt") = tp_dt_[time_pr_id];
  S_->GetRecordW("dt", Tags::NEXT, "dt").set_initialized();

  // Initialize
  S_->InitializeFields();
  S_->InitializeEvaluators();

  // Initialize the state from the old state.
  S_->Initialize(*S_old_);

  // Initialize the process kernels and verify
  pk_->Initialize();
  S_->CheckAllFieldsInitialized();

  S_->GetMeshPartition("materials");

  pk_->CalculateDiagnostics(Tags::DEFAULT);
  Observations();

  pk_->set_dt(tp_dt_[time_pr_id]);
  max_dt_ = tp_max_dt_[time_pr_id];
  
  // conditional i/o after initialization is performed only when 
  // new fields are added/removed to/from the state
  auto fields_old = StateVisFields(*S_old_);
  auto fields = StateVisFields(*S_);
  if (fields_old != fields) Visualize(true, make_tag("ic"));

  S_old_ = Teuchos::null;
}

}  // namespace Amanzi

