/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Implementation for the CycleDriver.  CycleDriver is basically just a class to hold
the cycle driver, which runs the overall, top level timestep loop.  It
instantiates states, ensures they are initialized, and runs the timestep loop
including Vis and restart/checkpoint dumps.  It contains one and only one PK
-- most likely this PK is an MPC of some type -- to do the actual work.
------------------------------------------------------------------------- */

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

#include "CycleDriver.hh"

#define DEBUG_MODE 1

namespace Amanzi {

CycleDriver::CycleDriver(Teuchos::ParameterList& parameter_list,
                         Teuchos::RCP<State>& S,
                         Epetra_MpiComm* comm,
                         Amanzi::ObservationData& output_observations) :
    parameter_list_(Teuchos::rcp(new Teuchos::ParameterList(parameter_list))),
    S_(S),
    comm_(comm),
    output_observations_(output_observations),
    restart_(false) {

  // create and start the global timer
  coordinator_init_();

  vo_ = Teuchos::rcp(new VerboseObject("CycleDriver", *parameter_list_));
};

void CycleDriver::coordinator_init_() {
  coordinator_list_ = Teuchos::sublist(parameter_list_, "Cycle Driver");
  read_parameter_list_();

  // create the global solution vector
  soln_ = Teuchos::rcp(new TreeVector());

  // create the pk tree root node (which then creates the rest of the tree)
  PKFactory pk_factory;
  Teuchos::RCP<Teuchos::ParameterList> pk_tree_list = Teuchos::sublist(parameter_list_, "PK Tree");
  if (pk_tree_list->numParams() == 0 || pk_tree_list->numParams() > 1) {
    Errors::Message message("CycleDriver: PK Tree list should contain exactly one root node list");
    Exceptions::amanzi_throw(message);
  }
  Teuchos::ParameterList::ConstIterator pk_item = pk_tree_list->begin();
  const std::string &pk_name = pk_tree_list->name(pk_item);
  if (!pk_tree_list->isSublist(pk_name)) {
    Errors::Message message("CycleDriver: PK Tree list should contain exactly one root node list");
    Exceptions::amanzi_throw(message);
  }

  pk_ = pk_factory.CreatePK(pk_tree_list->sublist(pk_name), parameter_list_, S_, soln_);
  pk_->Setup();

  // create the observations
  if (parameter_list_->isSublist("Observation Data")) {
    Teuchos::ParameterList observation_plist = parameter_list_->sublist("Observation Data");
    observations_ = Teuchos::rcp(new Amanzi::Unstructured_observations(observation_plist, output_observations_, comm_));

    if (coordinator_list_->isSublist("component names")) {
      Teuchos::Array<std::string> comp_names =
          coordinator_list_->get<Teuchos::Array<std::string> >("component names");
      observations_->RegisterComponentNames(comp_names.toVector());
    }
  }

  // create the checkpointing
  Teuchos::ParameterList& chkp_plist = parameter_list_->sublist("Checkpoint Data");
  checkpoint_ = Teuchos::rcp(new Checkpoint(chkp_plist, comm_));

  // create the time step manager
  tsm_ = Teuchos::rcp(new TimeStepManager());
}

void CycleDriver::setup() {
  // Set up the states, creating all data structures.
  S_->set_time(t0_);
  S_->set_cycle(cycle0_);
  S_->RequireScalar("dt", "coordinator");
  S_->Setup();
}

void CycleDriver::initialize() {
  // Restart from checkpoint, part 1.

  // This is crufty -- blame the BDF1 time integrator, whose solution history
  // needs to be updated to accept the new time as its initial time.

  // Note that if this is so, we can probably ignore some of the above
  // initialize() calls and the commit_state() call, but I'm afraid to try
  // that and break all the PKs.
  // Currently not a true restart -- for a true restart this should also get:
  // -- timestep size dt
  // -- BDF history to allow projection to continue correctly.
  if (restart_) {
    t0_ = ReadCheckpointInitialTime(comm_, restart_filename_);
    S_->set_time(t0_);
  }

  // Restart from checkpoint, part 2.
  if (restart_) {
    ReadCheckpoint(comm_, S_.ptr(), restart_filename_);
    t0_ = S_->time();
    cycle0_ = S_->cycle();

    DeformCheckpointMesh(S_.ptr());
  }

  // Initialize the process kernels (initializes all independent variables)
  pk_->Initialize();
  *S_->GetScalarData("dt", "coordinator") = 0.;
  S_->GetField("dt","coordinator")->set_initialized();

  // Initialize the state (initializes all dependent variables).
  S_->Initialize();

  // commit the initial conditions.
  pk_->CommitStep(0., 0.);

  // vis for the state
  // HACK to vis with a surrogate surface mesh.  This needs serious re-design. --etc
  bool surface_done = false;
  if (S_->HasMesh("surface") && S_->HasMesh("surface_3d")) {
    Teuchos::RCP<const AmanziMesh::Mesh> surface_3d = S_->GetMesh("surface_3d");
    Teuchos::RCP<const AmanziMesh::Mesh> surface = S_->GetMesh("surface");

    // vis successful timesteps
    std::string plist_name = "Visualization Data surface";
    Teuchos::ParameterList& vis_plist = parameter_list_->sublist(plist_name);
    Teuchos::RCP<Visualization> vis = Teuchos::rcp(new Visualization(vis_plist, comm_));
    vis->set_mesh(surface_3d);
    vis->CreateFiles();
    vis->set_mesh(surface);
    visualization_.push_back(vis);
    surface_done = true;

    // vis unsuccesful timesteps
    std::string fail_plist_name = "Visualization Data surface Failed Steps";
    if (parameter_list_->isSublist(fail_plist_name)) {
      Teuchos::ParameterList& fail_vis_plist = parameter_list_->sublist(fail_plist_name);
      Teuchos::RCP<Visualization> fail_vis = Teuchos::rcp(new Visualization(fail_vis_plist, comm_));
      fail_vis->set_mesh(surface_3d);
      fail_vis->CreateFiles();
      fail_vis->set_mesh(surface);
      failed_visualization_.push_back(fail_vis);
    }
  }

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

  // make observations
  observations_->MakeObservations(*S_);

  S_->set_time(t0_); // in case steady state solve changed this
  S_->set_cycle(cycle0_);

  // set up the TSM
  // -- register visualization times
  for (std::vector<Teuchos::RCP<Visualization> >::iterator vis=visualization_.begin();
       vis!=visualization_.end(); ++vis) {
    (*vis)->RegisterWithTimeStepManager(tsm_.ptr());
  }

  // -- register checkpoint times
  checkpoint_->RegisterWithTimeStepManager(tsm_.ptr());

  // -- register observation times
  observations_->RegisterWithTimeStepManager(tsm_.ptr());

  // -- register the final time
  tsm_->RegisterTimeEvent(t1_);

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

    double mem = rss_usage();
    
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
  t0_ = coordinator_list_->get<double>("start time");
  t1_ = coordinator_list_->get<double>("end time");
  std::string t0_units = coordinator_list_->get<std::string>("start time units", "s");
  std::string t1_units = coordinator_list_->get<std::string>("end time units", "s");

  if (t0_units == "s") {
    // internal units in s
  } else if (t0_units == "d") { // days
    t0_ = t0_ * 24.0*3600.0;
  } else if (t0_units == "yr") { // years
    t0_ = t0_ * 365.25*24.0*3600.0;
  } else {
    Errors::Message message("CycleDriver: error, invalid start time units");
    Exceptions::amanzi_throw(message);
  }

  if (t1_units == "s") {
    // internal units in s
  } else if (t1_units == "d") { // days
    t1_ = t1_ * 24.0*3600.0;
  } else if (t1_units == "yr") { // years
    t1_ = t1_ * 365.25*24.0*3600.0;
  } else {
    Errors::Message message("CycleDriver: error, invalid end time units");
    Exceptions::amanzi_throw(message);
  }

  max_dt_ = coordinator_list_->get<double>("max time step size", 1.0e99);
  min_dt_ = coordinator_list_->get<double>("min time step size", 1.0e-12);
  cycle0_ = coordinator_list_->get<int>("start cycle",0);
  cycle1_ = coordinator_list_->get<int>("end cycle",-1);

  // restart control
  restart_ = coordinator_list_->isSublist("Restart from Checkpoint Data File");
  if (restart_) {
    Teuchos::ParameterList restart_list = coordinator_list_->sublist("Restart from Checkpoint Data File");
    restart_filename_ = coordinator_list_->get<std::string>("Checkpoint Data File Name");

    // make sure that the restart file actually exists, if not throw an error
    boost::filesystem::path restart_from_filename_path(restart_filename_);
    if (!boost::filesystem::exists(restart_from_filename_path)) {
      Errors::Message message("MPC: the specified restart file does not exist or is not a regular file.");
      Exceptions::amanzi_throw(message);
    }
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
  dt = tsm_->TimeStep(S_->time() + dt, dt);
  return dt;
}



// This is used by CLM
bool CycleDriver::advance(double dt) {
  S_->advance_time(dt);
  bool fail = pk_->AdvanceStep(S_->last_time(), S_->time());

  if (!fail) {
    // advance the iteration count and timestep size
    S_->advance_cycle();

    // make observations, vis, and checkpoints
    observations_->MakeObservations(*S_);
    visualize();
    checkpoint(dt);

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

  if (dump) {
    pk_->CalculateDiagnostics();
  }

  for (std::vector<Teuchos::RCP<Visualization> >::iterator vis=visualization_.begin();
       vis!=visualization_.end(); ++vis) {
    if (force || (*vis)->DumpRequested(S_->cycle(), S_->time())) {
      WriteVis((*vis).ptr(), S_.ptr());
    }
  }
}

void CycleDriver::checkpoint(double dt, bool force) {
  if (force || checkpoint_->DumpRequested(S_->cycle(), S_->time())) {
    WriteCheckpoint(checkpoint_.ptr(), S_.ptr(), dt);
  }
}


// -----------------------------------------------------------------------------
// timestep loop
// -----------------------------------------------------------------------------
void CycleDriver::go() {
  // start at time t = t0 and initialize the state.
  setup();
  initialize();

  // get the intial timestep -- note, this would have to be fixed for a true restart
  double dt = get_dt();

  // visualization at IC
  visualize();
  checkpoint(dt);

  // iterate process kernels
  {
#if !DEBUG_MODE
  try {
#endif
    bool fail = false;
    while ((S_->time() < t1_) &&
           ((cycle1_ == -1) || (S_->cycle() <= cycle1_))) {
      if (vo_->os_OK(Teuchos::VERB_MEDIUM)) {
        Teuchos::OSTab tab = vo_->getOSTab();
        *vo_->os() << "======================================================================"
                  << std::endl << std::endl;
        *vo_->os() << "Cycle = " << S_->cycle();
        *vo_->os() << ",  Time [days] = "<< S_->time() / (60*60*24);
        *vo_->os() << ",  dt [days] = " << dt / (60*60*24)  << std::endl;
        *vo_->os() << "----------------------------------------------------------------------"
                  << std::endl;
      }

      *S_->GetScalarData("dt", "coordinator") = dt;
      fail = advance(dt);
      dt = get_dt();

    } // while not finished


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
  finalize();

} // cycle driver

} // close namespace Amanzi
