/*
  This is the MPC component of the Amanzi code. 

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

#include <iostream>
#include <unistd.h>
#include <sys/resource.h>
#include "errors.hh"

#include "Teuchos_VerboseObjectParameterListHelpers.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_TimeMonitor.hpp"


#define BOOST_FILESYSTEM_NO_DEPRECATED
#include "boost/filesystem/operations.hpp"
#include "boost/filesystem/path.hpp"

#include "TimeStepManager.hh"
#include "visualization.hh"
#include "checkpoint.hh"
#include "ObservationData.hh"
#include "Unstructured_observations.hh"
#include "State.hh"
#include "PK.hh"
#include "TreeVector.hh"
#include "PK_Factory.hh"
#include "TimerManager.hh"
#include "CycleDriver.hh"

#define DEBUG_MODE 1

namespace Amanzi {

bool reset_info_compfunc(std::pair<double,double> x, std::pair<double,double> y) {
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


CycleDriver::CycleDriver(Teuchos::RCP<Teuchos::ParameterList> glist_,
                         Teuchos::RCP<State>& S,
                         Epetra_MpiComm* comm,
                         Amanzi::ObservationData& output_observations) :
    parameter_list_(glist_),
    S_(S),
    comm_(comm),
    output_observations_(output_observations),
    restart_requested_(false) {

  // create and start the global timer
  coordinator_init_();

  vo_ = Teuchos::rcp(new VerboseObject("CycleDriver", parameter_list_->sublist("Cycle Driver")));
};


void CycleDriver::coordinator_init_() {
  coordinator_list_ = Teuchos::sublist(parameter_list_, "Cycle Driver");
  read_parameter_list_();

  // create the global solution vector
  soln_ = Teuchos::rcp(new TreeVector());
}


void CycleDriver::init_pk(int time_pr_id){
  // create the pk tree root node (which then creates the rest of the tree)
  PKFactory pk_factory;

  Teuchos::RCP<Teuchos::ParameterList> time_periods_list = Teuchos::sublist(coordinator_list_, "time periods", true);

  std::ostringstream ss; ss << time_pr_id;
  std::string tp_list_name = "TP "+ ss.str();
  Teuchos::RCP<Teuchos::ParameterList> tp_list =Teuchos::sublist( time_periods_list, tp_list_name.data(), true);
  Teuchos::ParameterList pk_tree_list = tp_list->sublist("PK Tree");
  if (pk_tree_list.numParams() != 1) {
    Errors::Message message("CycleDriver: PK Tree list should contain exactly one root node list");
    Exceptions::amanzi_throw(message);
  }
  Teuchos::ParameterList::ConstIterator pk_item = pk_tree_list.begin();
  const std::string &pk_name = pk_tree_list.name(pk_item);

  if (!pk_tree_list.isSublist(pk_name)) {
    Errors::Message message("CycleDriver: PK Tree list does not have node \"" + pk_name + "\".");
    Exceptions::amanzi_throw(message);
  }

  pk_ = pk_factory.CreatePK(pk_tree_list.sublist(pk_name), parameter_list_, S_, soln_);
}


void CycleDriver::setup() {
  // Set up the states, creating all data structures.

  // create the observations
  if (parameter_list_->isSublist("Observation Data")) {
    Teuchos::ParameterList observation_plist = parameter_list_->sublist("Observation Data");
    observations_ = Teuchos::rcp(new Amanzi::Unstructured_observations(observation_plist, output_observations_, comm_));
    //std::cout<<*coordinator_list_<<"\n";
    if (coordinator_list_->isParameter("component names")) {
      Teuchos::Array<std::string> comp_names =
          coordinator_list_->get<Teuchos::Array<std::string> >("component names");
      observations_->RegisterComponentNames(comp_names.toVector());
    }
  }

  // create the checkpointing
  if (parameter_list_->isSublist("Checkpoint Data")) {
    Teuchos::ParameterList& chkp_plist = parameter_list_->sublist("Checkpoint Data");
    checkpoint_ = Teuchos::rcp(new Amanzi::Checkpoint(chkp_plist, comm_));
  }
  else{
    checkpoint_ = Teuchos::rcp(new Amanzi::Checkpoint());
  }

  // create the time step manager
  tsm_ = Teuchos::ptr(new TimeStepManager());

  pk_->Setup();

  S_->RequireScalar("dt", "coordinator");

  S_->Setup();

  // S_->InitializeFields();
  // S_->InitializeEvaluators();

  if (vo_->os_OK(Teuchos::VERB_MEDIUM)) {
    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "Setup is complete." << std::endl;
  }
}


void CycleDriver::initialize() {
  // register observation times with the time step manager
 
  *S_->GetScalarData("dt", "coordinator") = tp_dt_[0];
  S_->GetField("dt", "coordinator")->set_initialized();

  // Initialize the state (initializes all dependent variables).
  S_->InitializeFields();
  S_->InitializeEvaluators();

  // Initialize the process kernels (initializes all independent variables)
  pk_->Initialize();

  // Final checks.
  S_->CheckNotEvaluatedFieldsInitialized();
  S_->CheckAllFieldsInitialized();


  // S_->WriteDependencyGraph();

  S_->GetMeshPartition("materials");

  // commit the initial conditions.
  // pk_->CommitStep(t0_-get_dt(), get_dt());
  if (!restart_requested_) {
    pk_->CommitStep(S_->time(), S_->time());
    // visualize();
    // checkpoint(*S_->GetScalarData("dt", "coordinator"));
  }

  // // vis for the state
  // // HACK to vis with a surrogate surface mesh.  This needs serious re-design. --etc
  bool surface_done = false;
  // if (S_->HasMesh("surface") && S_->HasMesh("surface_3d")) {
  //   Teuchos::RCP<const AmanziMesh::Mesh> surface_3d = S_->GetMesh("surface_3d");
  //   Teuchos::RCP<const AmanziMesh::Mesh> surface = S_->GetMesh("surface");

  //   // vis successful timesteps
  //   std::string plist_name = "Visualization Data surface";
  //   Teuchos::ParameterList& vis_plist = parameter_list_->sublist(plist_name);
  //   Teuchos::RCP<Visualization> vis = Teuchos::rcp(new Visualization(vis_plist, comm_));
  //   vis->set_mesh(surface_3d);
  //   vis->CreateFiles();
  //   vis->set_mesh(surface);
  //   visualization_.push_back(vis);
  //   surface_done = true;

  //   // vis unsuccesful timesteps
  //   std::string fail_plist_name = "Visualization Data surface Failed Steps";
  //   if (parameter_list_->isSublist(fail_plist_name)) {
  //     Teuchos::ParameterList& fail_vis_plist = parameter_list_->sublist(fail_plist_name);
  //     Teuchos::RCP<Visualization> fail_vis = Teuchos::rcp(new Visualization(fail_vis_plist, comm_));
  //     fail_vis->set_mesh(surface_3d);
  //     fail_vis->CreateFiles();
  //     fail_vis->set_mesh(surface);
  //     failed_visualization_.push_back(fail_vis);
  //   }
  // }

  for (State::mesh_iterator mesh=S_->mesh_begin();
       mesh!=S_->mesh_end(); ++mesh) {
    if (mesh->first == "surface_3d") {
      // pass
    } else if ((mesh->first == "surface") && surface_done) {
      // pass
    } else {
      // vis successful steps
      std::string plist_name = "Visualization Data "+mesh->first;
      // in the case of just a domain mesh, we want to allow no name.
      if ((mesh->first == "domain") && !parameter_list_->isSublist(plist_name)) {
        plist_name = "Visualization Data";
      }

      if (parameter_list_->isSublist(plist_name)) {
        Teuchos::ParameterList& vis_plist = parameter_list_->sublist(plist_name);
        Teuchos::RCP<Visualization> vis =
          Teuchos::rcp(new Visualization(vis_plist, comm_));
        vis->set_mesh(mesh->second.first);
        vis->CreateFiles();
        visualization_.push_back(vis);
      }

      // vis unsuccessful steps
      std::string fail_plist_name = "Visualization Data "+mesh->first+" Failed Steps";
      // in the case of just a domain mesh, we want to allow no name.
      if ((mesh->first == "domain") && !parameter_list_->isSublist(fail_plist_name)) {
        fail_plist_name = "Visualization Data Failed Steps";
      }

      if (parameter_list_->isSublist(fail_plist_name)) {
        Teuchos::ParameterList& fail_vis_plist = parameter_list_->sublist(fail_plist_name);
        Teuchos::RCP<Visualization> fail_vis =
          Teuchos::rcp(new Visualization(fail_vis_plist, comm_));
        fail_vis->set_mesh(mesh->second.first);
        fail_vis->CreateFiles();
        failed_visualization_.push_back(fail_vis);
      }
    }
  }

  // set up the TSM
  // -- register visualization times
  for (std::vector<Teuchos::RCP<Visualization> >::iterator vis=visualization_.begin();
       vis!=visualization_.end(); ++vis) {
    (*vis)->RegisterWithTimeStepManager(tsm_.ptr());
  }
  // -- register checkpoint times
  checkpoint_->RegisterWithTimeStepManager(tsm_.ptr());
  // -- register observation times
  if (observations_ != Teuchos::null) 
    observations_->RegisterWithTimeStepManager(tsm_.ptr());
  // -- register the final time
  // register reset_times
  for(std::vector<std::pair<double,double> >::const_iterator it = reset_info_.begin();
      it != reset_info_.end(); ++it) tsm_->RegisterTimeEvent(it->first);

  for (int i=0;i<num_time_periods_; i++) tsm_->RegisterTimeEvent(tp_end_[i]);
  
  //tsm_->RegisterTimeEvent(t1_);

}


void CycleDriver::finalize() {
  // Force checkpoint at the end of simulation.
  // Only do if the checkpoint was not already written, or we would be writing
  // the same file twice.
  // This really should be removed, but for now is left to help stupid developers.
  if (!checkpoint_->DumpRequested(S_->cycle(), S_->time())) {
    pk_->CalculateDiagnostics();
    WriteCheckpoint(checkpoint_.ptr(), S_.ptr(), 0.0);
  }
}


// double rss_usage() { // return ru_maxrss in MBytes
// #if (defined(__unix__) || defined(__unix) || defined(unix) || defined(__APPLE__) || defined(__MACH__))
//   struct rusage usage;
//   getrusage(RUSAGE_SELF, &usage);
// #if (defined(__APPLE__) || defined(__MACH__))
//   return static_cast<double>(usage.ru_maxrss)/1024.0/1024.0;
// #else
//   return static_cast<double>(usage.ru_maxrss)/1024.0;
// #endif
// #else
//   return 0.0;
// #endif
// }


void CycleDriver::report_memory() {
  // report the memory high water mark (using ru_maxrss)
  // this should be called at the very end of a simulation
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
    *vo_->os() << "======================================================================" << std::endl;
    *vo_->os() << "All meshes combined have " << global_ncells << " cells." << std::endl;
    *vo_->os() << "Memory usage (high water mark):" << std::endl;
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
  }

  
  double doubles_count(0.0);
  for (State::field_iterator field=S_->field_begin(); field!=S_->field_end(); ++field) {
    doubles_count += static_cast<double>(field->second->GetLocalElementCount());
  }
  double global_doubles_count(0.0);
  double min_doubles_count(0.0);
  double max_doubles_count(0.0);
  comm_->SumAll(&doubles_count,&global_doubles_count,1);
  comm_->MinAll(&doubles_count,&min_doubles_count,1);
  comm_->MaxAll(&doubles_count,&max_doubles_count,1);

  Teuchos::OSTab tab = vo_->getOSTab();
  *vo_->os() << "Doubles allocated in state fields " << std::endl;
  *vo_->os() << "  Maximum per core:   " << std::setw(7)
             << max_doubles_count*8/1024/1024 << " MBytes" << std::endl;
  *vo_->os() << "  Minimum per core:   " << std::setw(7)
             << min_doubles_count*8/1024/1024 << " MBytes" << std::endl; 
  *vo_->os() << "  Total:              " << std::setw(7)
             << global_doubles_count*8/1024/1024 << " MBytes" << std::endl;
}


void CycleDriver::read_parameter_list_() {
  //  std::cout<<*coordinator_list_<<"\n";
  // t0_ = coordinator_list_->get<double>("start time");
  // t1_ = coordinator_list_->get<double>("end time");
  // std::string t0_units = coordinator_list_->get<std::string>("start time units", "s");
  // std::string t1_units = coordinator_list_->get<std::string>("end time units", "s");

  // if (t0_units == "s") {
  //   // internal units in s
  // } else if (t0_units == "d") { // days
  //   t0_ = t0_ * 24.0*3600.0;
  // } else if (t0_units == "yr") { // years
  //   t0_ = t0_ * 365.25*24.0*3600.0;
  // } else {
  //   Errors::Message message("CycleDriver: error, invalid start time units");
  //   Exceptions::amanzi_throw(message);
  // }

  // if (t1_units == "s") {
  //   // internal units in s
  // } else if (t1_units == "d") { // days
  //   t1_ = t1_ * 24.0*3600.0;
  // } else if (t1_units == "yr") { // years
  //   t1_ = t1_ * 365.25*24.0*3600.0;
  // } else {
  //   Errors::Message message("CycleDriver: error, invalid end time units");
  //   Exceptions::amanzi_throw(message);
  // }

  max_dt_ = coordinator_list_->get<double>("max time step size", 1.0e99);
  min_dt_ = coordinator_list_->get<double>("min time step size", 1.0e-12);
  cycle0_ = coordinator_list_->get<int>("start cycle",0);
  cycle1_ = coordinator_list_->get<int>("end cycle",-1);

  Teuchos::ParameterList time_periods_list = coordinator_list_->sublist("time periods");

  num_time_periods_ = time_periods_list.numParams();
  Teuchos::ParameterList::ConstIterator item;
  tp_start_.resize(num_time_periods_);
  tp_end_.resize(num_time_periods_);
  tp_dt_.resize(num_time_periods_);
  tp_max_cycle_.resize(num_time_periods_);

  int i = 0;
  for (item = time_periods_list.begin(); item !=time_periods_list.end(); ++item){
    const std::string & tp_name = time_periods_list.name(item);
    tp_start_[i] = time_periods_list.sublist(tp_name).get<double>("start period time");
    tp_end_[i] = time_periods_list.sublist(tp_name).get<double>("end period time");
    tp_dt_[i] = time_periods_list.sublist(tp_name).get<double>("initial time step", 1);
    tp_max_cycle_[i] = time_periods_list.sublist(tp_name).get<int>("maximum cycle number", -1);
   
    std::string t_units = time_periods_list.sublist(tp_name).get<std::string>("start time units", "s");
    if (t_units == "s") {
      // internal units in s
    } else if (t_units == "d") {  // days
      tp_start_[i] = tp_start_[i] * 24.0*3600.0;
    } else if (t_units == "yr") {  // years
      tp_start_[i] = tp_start_[i] * 365.25*24.0*3600.0;
    } else {
      Errors::Message message("CycleDriver: error, invalid start time units");
      Exceptions::amanzi_throw(message);
    }
    t_units = time_periods_list.sublist(tp_name).get<std::string>("end time units", "s");
    if (t_units == "s") {
      // internal units in s
    } else if (t_units == "d") {  // days
      tp_end_[i] = tp_end_[i] * 24.0*3600.0;
    } else if (t_units == "yr") {  // years
      tp_end_[i] = tp_end_[i] * 365.25*24.0*3600.0;
    } else {
      Errors::Message message("CycleDriver: error, invalid end time units");
      Exceptions::amanzi_throw(message);
    }
    t_units = time_periods_list.sublist(tp_name).get<std::string>("initial time step units", "s");
    if (t_units == "s") {
      // internal units in s
    } else if (t_units == "d") {  // days
      tp_dt_[i] = tp_dt_[i] * 24.0*3600.0;
    } else if (t_units == "yr") {  // years
      tp_dt_[i] = tp_dt_[i] * 365.25*24.0*3600.0;
    } else {
      Errors::Message message("CycleDriver: error, invalid initial time step time units");
      Exceptions::amanzi_throw(message);
    }
    i++;
  }

  // restart control
  // are we restarting from a file?
  // first assume we're not
  restart_requested_ = false;

  if (coordinator_list_->isSublist("Restart from Checkpoint Data File")) {
    restart_requested_ = ! coordinator_list_->sublist("Restart from Checkpoint Data File")
        .get<bool>("initialize from checkpoint data file and do not restart",false);

    if (restart_requested_) {
      Teuchos::ParameterList restart_list = coordinator_list_->sublist("Restart from Checkpoint Data File");
      restart_filename_ = restart_list.get<std::string>("Checkpoint Data File Name");

      // make sure that the restart file actually exists, if not throw an error
      boost::filesystem::path restart_from_filename_path(restart_filename_);
      if (!boost::filesystem::exists(restart_from_filename_path)) {
        Errors::Message message("CycleDriver::the specified restart file does not exist or is not a regular file.");
        Exceptions::amanzi_throw(message);
      }
    }

    // if (restart_requested_) {
    //   if (vo_->os_OK(Teuchos::VERB_LOW)) {
    //     *vo_->os() << "Restarting from checkpoint file: " << restart_filename_ << std::endl;
    //   }
    // } else {
    //   if (vo_->os_OK(Teuchos::VERB_LOW)) {
    //     *vo_->os() << "Initializing data from checkpoint file: " << restart_filename_ << std::endl
    //                << "    (Ignoring all initial conditions.)" << std::endl;
    //   }
    // }
  }
}


// -----------------------------------------------------------------------------
// acquire the chosen timestep size
// -----------------------------------------------------------------------------
double CycleDriver::get_dt() {
  // get the physical step size
  double dt = pk_->get_dt();

  // check if the step size has gotten too small
  if (dt < min_dt_) {
    Errors::Message message("CycleDriver: error, timestep too small");
    Exceptions::amanzi_throw(message);
  }

  // cap the max step size
  if (dt > max_dt_) {
    dt = max_dt_;
  }

  // ask the step manager if this step is ok
  dt = tsm_->TimeStep(S_->time(), dt);
  return dt;
}


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
  dt_ = tsm_->TimeStep(S_->time() + dt, dt);

  // set the physical step size
  pk_->set_dt(dt_);
}


// This is used by CLM
bool CycleDriver::advance(double dt) {
  Teuchos::OSTab tab = vo_->getOSTab();
  //bool fail = pk_->AdvanceStep(S_->last_time(), S_->time());
  bool fail = pk_->AdvanceStep(S_->time(), S_->time()+dt);

  if (!fail) {
    pk_->CommitStep(S_->last_time(), S_->time());
    // advance the iteration count and timestep size
    S_->advance_cycle();
    S_->advance_time(dt);

    if (vo_->os_OK(Teuchos::VERB_MEDIUM)) {
      *vo_->os() << "New time(y) = "<< S_->time() / (60*60*24*365.25);
      *vo_->os() << std::endl;
    }

    bool force_vis(false);
    bool force_check(false);
    bool force_obser(false);

    if (abs(S_->time() - tp_end_[time_period_id_]) < 1e-7){
      force_vis = true;
      force_check = true;                       
      force_obser = true;
    }

    if (!reset_info_.empty())
        if (S_->time() == reset_info_.front().first)
            force_check = true;

    // make observations, vis, and checkpoints
    //Amanzi::timer_manager.start("I/O");
    observations(force_obser);
    visualize(force_vis);
    checkpoint(dt, force_check);
    //Amanzi::timer_manager.start("I/O");

    // we're done with this time step, advance the state
    // NOT YET IMPLEMENTED, requires PKs to deal with this in CommitStep().
    // This breaks coupling in general, but fixing it requires changes to both
    // the PKs and the State. --ETC

  } else {
    // Failed the timestep.  
    // Potentially write out failed timestep for debugging
    for (std::vector<Teuchos::RCP<Visualization> >::iterator vis=failed_visualization_.begin();
         vis!=failed_visualization_.end(); ++vis) {
      WriteVis((*vis).ptr(), S_.ptr());
    }
    // The timestep sizes have been updated, so copy back old soln and try again.
    // NOT YET IMPLEMENTED, requires PKs to deal with failure.  Fortunately
    // transport and chemistry never fail, so we shouldn't break things.
    // Otherwise this would be very broken, as flow could succeed, but
    // transport fail, and we wouldn't have a way of backing up. --ETC
  }
  return fail;
}


void CycleDriver::observations(bool force) {
  if (observations_ != Teuchos::null) {
    if (observations_->DumpRequested(S_->cycle(), S_->time()) || force) {
      pk_->CalculateDiagnostics();
      observations_->MakeObservations(*S_);
    }
  }
}


void CycleDriver::visualize(bool force) {
  // write visualization if requested
  bool dump = force;
  if (!dump) {
    for (std::vector<Teuchos::RCP<Visualization> >::iterator vis=visualization_.begin();
         vis!=visualization_.end(); ++vis) {
      if ((*vis)->DumpRequested(S_->cycle(), S_->time())) {
        dump = true;
      }
    }
  }

  if (dump || force) {
    pk_->CalculateDiagnostics();
  }

  for (std::vector<Teuchos::RCP<Visualization> >::iterator vis=visualization_.begin();
       vis!=visualization_.end(); ++vis) {
    if (force || (*vis)->DumpRequested(S_->cycle(), S_->time())) {
      WriteVis((*vis).ptr(), S_.ptr());

      Teuchos::OSTab tab = vo_->getOSTab();
      *vo_->os() << "Cycle " << S_->cycle() << ": writing visualization file" << std::endl;
    }
  }
}


void CycleDriver::checkpoint(double dt, bool force) {
  if (force || checkpoint_->DumpRequested(S_->cycle(), S_->time())) {
    WriteCheckpoint(checkpoint_.ptr(), S_.ptr(), dt);
    pk_->CalculateDiagnostics();

    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "Cycle " << S_->cycle() << ": writing checkpoint file" << std::endl;
  }
}


// -----------------------------------------------------------------------------
// timestep loop
// -----------------------------------------------------------------------------
void CycleDriver::go() {
  time_period_id_ = 0;
  init_pk(time_period_id_);

  // start at time t = t0 and initialize the state.
  setup();
  initialize();

  S_->set_time(tp_start_[time_period_id_]);
  S_->set_cycle(cycle0_);

  double dt;
  double restart_dT(1.0e99);
  if (restart_requested_) {
    // re-initialize the state object
    restart_dT = ReadCheckpoint(comm_, Teuchos::ptr(&*S_), restart_filename_);
    cycle0_ = S_->cycle();
    for (std::vector<std::pair<double,double> >::iterator it = reset_info_.begin();
          it != reset_info_.end(); ++it) {
      if (it->first < S_->time()) it = reset_info_.erase(it);
      if (it == reset_info_.end() ) break;
    }
    S_->set_initial_time(S_->time());

    // Restart time after the first time period create new PK
    for (int i=0;i<num_time_periods_;i++){
      if (tp_end_[i] <=S_->time()) time_period_id_++;
    }
    if (time_period_id_ > 0){
        reset_driver(time_period_id_); 
    }
    
    dt = restart_dT;
    pk_->set_dt(dt);
  }
  else {    
    dt = tp_dt_[time_period_id_];
    pk_->set_dt(dt);
  }

  
  *S_->GetScalarData("dt", "coordinator") = dt;
  S_->GetField("dt","coordinator")->set_initialized();

  // visualization at IC
  //Amanzi::timer_manager.start("I/O");
  pk_->CalculateDiagnostics();
  visualize();
  checkpoint(dt);
  //Amanzi::timer_manager.stop("I/O");
 
  // iterate process kernels
  {
#if !DEBUG_MODE
  try {
#endif
    bool fail = false;
    while (time_period_id_ < num_time_periods_){
      int start_cycle_num = S_->cycle();
      do {
        if (vo_->os_OK(Teuchos::VERB_MEDIUM)) {
          Teuchos::OSTab tab = vo_->getOSTab();
          *vo_->os() << "\nCycle " << S_->cycle()
                     << ": time(y) = " << S_->time() / (60*60*24*365.25)
                     << ", dt(y) = " << dt / (60*60*24*365.25) << std::endl;
        }
        *S_->GetScalarData("dt", "coordinator") = dt;
	S_->set_initial_time(S_->time());
	S_->set_final_time(S_->time() + dt);
        fail = advance(dt);
        dt = get_dt();
      }  // while not finished
      while ((S_->time() < tp_end_[time_period_id_]) &&  ((tp_max_cycle_[time_period_id_] == -1) 
                                     || (S_->cycle() - start_cycle_num <= tp_max_cycle_[time_period_id_])));
      time_period_id_++;
      if (time_period_id_ < num_time_periods_){
        reset_driver(time_period_id_); 
        dt = get_dt();
      }      
    }

#if !DEBUG_MODE
  }

  catch (Exceptions::Amanzi_exception &e) {
    // write one more vis for help debugging
    S_->advance_cycle();
    visualize(true); // force vis

    // flush observations to make sure they are saved
    observations_->Flush();

    // catch errors to dump two checkpoints -- one as a "last good" checkpoint
    // and one as a "debugging data" checkpoint.
    checkpoint_->set_filebasename("error_checkpoint");
    WriteCheckpoint(checkpoint_.ptr(), S_.ptr(), dt);
    throw e;
  }
#endif
  }
  
  report_memory();
  //finalize();
} 


void CycleDriver::reset_driver(int time_pr_id) {
  Teuchos::RCP<AmanziMesh::Mesh> mesh = Teuchos::rcp_const_cast<AmanziMesh::Mesh>(S_->GetMesh("domain"));
  S_old_ = S_;

  Teuchos::ParameterList state_plist = parameter_list_->sublist("State");
  S_ = Teuchos::rcp(new Amanzi::State(state_plist));
  S_->RegisterMesh("domain", mesh);

  //delete the old global solution vector
  // soln_ = Teuchos::null;
  // pk_ = Teuchos::null;
  
  //if (pk_.get()) delete pk_.get(); 
  pk_ = Teuchos::null;

  //if (soln_.get()) delete soln_.get(); 
  soln_ = Teuchos::null;

  // create the global solution vector
  soln_ = Teuchos::rcp(new TreeVector());
  
  // create new pk;
  init_pk(time_pr_id);

  // register observation times with the time step manager
  //if (observations_ != Teuchos::null) observations_->RegisterWithTimeStepManager(tsm_);

  // Setup
  pk_->Setup();

  S_->RequireScalar("dt", "coordinator");
  S_->Setup();

  // Initialize
  S_->InitializeFields();
  S_->InitializeEvaluators();

  // Initialize the state from the old state.
  S_->Initialize(S_old_);

  // Initialize the process kernels variables 
  pk_->Initialize();

  // Final checks
  S_->CheckNotEvaluatedFieldsInitialized();
  S_->CheckAllFieldsInitialized();

  S_->GetMeshPartition("materials");

  S_->set_cycle(S_old_->cycle());
  S_->set_time(tp_start_[time_pr_id]); 
  pk_->set_dt(tp_dt_[time_pr_id]);

  S_old_ = Teuchos::null;

  // WriteCheckpoint(checkpoint_.ptr(), S_.ptr(), tp_dt_[time_pr_id]);
  // exit(0);
}

}  // namespace Amanzi

