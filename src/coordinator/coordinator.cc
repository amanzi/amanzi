/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Implementation for the Coordinator.  Coordinator is basically just a class to hold
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

#include "TimeStepManager.hh"
#include "Visualization.hh"
#include "checkpoint.hh"
#include "UnstructuredObservations.hh"
#include "State.hh"
#include "PK.hh"
#include "TreeVector.hh"
#include "PK_Factory.hh"
//#include "pk_factory_ats.hh"

#include "coordinator.hh"

#define DEBUG_MODE 1

namespace ATS {

Coordinator::Coordinator(Teuchos::ParameterList& parameter_list,
                         Teuchos::RCP<Amanzi::State>& S,
                         Epetra_MpiComm* comm ) :
    parameter_list_(Teuchos::rcp(new Teuchos::ParameterList(parameter_list))),
    S_(S),
    comm_(comm),
    restart_(false) {

  // create and start the global timer
  timer_ = Teuchos::rcp(new Teuchos::Time("wallclock_monitor",true));
  setup_timer_ = Teuchos::TimeMonitor::getNewCounter("setup");
  cycle_timer_ = Teuchos::TimeMonitor::getNewCounter("cycle");
  coordinator_init();

  vo_ = Teuchos::rcp(new Amanzi::VerboseObject("Coordinator", *parameter_list_));
};

void Coordinator::coordinator_init() {
  coordinator_list_ = Teuchos::sublist(parameter_list_, "cycle driver");
  read_parameter_list();

  // create the top level PK
  Teuchos::RCP<Teuchos::ParameterList> pks_list = Teuchos::sublist(parameter_list_, "PKs");
  Teuchos::ParameterList pk_tree_list = coordinator_list_->sublist("PK tree");
  if (pk_tree_list.numParams() != 1) {
    Errors::Message message("CycleDriver: PK tree list should contain exactly one root node list");
    Exceptions::amanzi_throw(message);
  }
  Teuchos::ParameterList::ConstIterator pk_item = pk_tree_list.begin();
  const std::string &pk_name = pk_tree_list.name(pk_item);
  
  // create the solution
  soln_ = Teuchos::rcp(new Amanzi::TreeVector());

  // create the pk
  Amanzi::PKFactory pk_factory;
  Teuchos::RCP<Teuchos::ParameterList> pk_list = Teuchos::sublist(pks_list, pk_name);
  pk_list->set("PK name", pk_name);
  const std::string &pk_origin = pk_list -> get<std::string>("PK origin", "ATS");


  pk_ = pk_factory.CreatePK(pk_tree_list.sublist(pk_name), parameter_list_, S_, soln_);

  // create the checkpointing
  Teuchos::ParameterList& chkp_plist = parameter_list_->sublist("checkpoint");
  checkpoint_ = Teuchos::rcp(new Amanzi::Checkpoint(chkp_plist, comm_));

  // create the observations
  Teuchos::ParameterList& observation_plist = parameter_list_->sublist("observations");
  observations_ = Teuchos::rcp(new Amanzi::UnstructuredObservations(observation_plist,
          Teuchos::null, comm_));

  // check whether meshes are deformable, and if so require a nodal position
  for (Amanzi::State::mesh_iterator mesh=S_->mesh_begin();
       mesh!=S_->mesh_end(); ++mesh) {
    if (S_->IsDeformableMesh(mesh->first)) {
      std::string node_key = std::string("vertex_coordinate_")+mesh->first;
      S_->RequireField(node_key)->SetMesh(mesh->second.first)->SetGhosted()
          ->AddComponent("node", Amanzi::AmanziMesh::NODE, mesh->second.first->space_dimension());
    }
  }

  // create the time step manager
  tsm_ = Teuchos::rcp(new Amanzi::TimeStepManager());
  
}

void Coordinator::setup() {
  // Set up the states, creating all data structures.
  S_->set_time(t0_);
  S_->set_cycle(cycle0_);
  S_->RequireScalar("dt", "coordinator");

  pk_->Setup(S_.ptr());  
  S_->Setup();
}

void Coordinator::initialize() {
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
    t0_ = Amanzi::ReadCheckpointInitialTime(comm_, restart_filename_);
    S_->set_time(t0_);
  }

  // Restart from checkpoint, part 2.
  if (restart_) {
    ReadCheckpoint(comm_, S_.ptr(), restart_filename_);
    t0_ = S_->time();
    cycle0_ = S_->cycle();

    DeformCheckpointMesh(S_.ptr());
  }

  // Initialize the state (initializes all dependent variables).
  //S_->Initialize();
  *S_->GetScalarData("dt", "coordinator") = 0.;
  S_->GetField("dt","coordinator")->set_initialized();

  S_->InitializeFields();

  // Initialize the process kernels (initializes all independent variables)
  pk_->Initialize(S_.ptr());

 // Final checks.
  S_->CheckNotEvaluatedFieldsInitialized();

  S_->InitializeEvaluators();


  S_->CheckAllFieldsInitialized();

  S_->WriteStatistics(vo_);


  // commit the initial conditions.
  pk_->CommitStep(0., 0., S_);

  // vis for the state
  // HACK to vis with a surrogate surface mesh.  This needs serious re-design. --etc
  bool surface_done = false;
  if (S_->HasMesh("surface") && S_->HasMesh("surface_3d")) {
    Teuchos::RCP<const Amanzi::AmanziMesh::Mesh> surface_3d = S_->GetMesh("surface_3d");
    Teuchos::RCP<const Amanzi::AmanziMesh::Mesh> surface = S_->GetMesh("surface");

    // vis successful timesteps
    std::string plist_name = "visualization surface";
    Teuchos::ParameterList& vis_plist = parameter_list_->sublist(plist_name);
/*
    Teuchos::RCP<Amanzi::Visualization> vis_2d = Teuchos::rcp(new Visualization(vis_plist, comm_));
    vis_2d->set_mesh(surface);
    vis_2d->CreateFiles();
    vis_2d->set_mesh(surface);
    visualization_.push_back(vis_2d);
*/    
    Teuchos::RCP<Amanzi::Visualization> vis = Teuchos::rcp(new Amanzi::Visualization(vis_plist, comm_));
    vis->set_mesh(surface_3d);
    vis->CreateFiles();
    vis->set_mesh(surface);
    visualization_.push_back(vis);
    surface_done = true;

    // vis unsuccesful timesteps
    std::string fail_plist_name = "visualization surface failed steps";
    if (parameter_list_->isSublist(fail_plist_name)) {
      Teuchos::ParameterList& fail_vis_plist = parameter_list_->sublist(fail_plist_name);
      Teuchos::RCP<Amanzi::Visualization> fail_vis = Teuchos::rcp(new Amanzi::Visualization(fail_vis_plist, comm_));
      fail_vis->set_mesh(surface_3d);
      fail_vis->CreateFiles();
      fail_vis->set_mesh(surface);
      failed_visualization_.push_back(fail_vis);
    }
  }

  for (Amanzi::State::mesh_iterator mesh=S_->mesh_begin();
       mesh!=S_->mesh_end(); ++mesh) {
    if (mesh->first == "surface_3d") {
      // pass
    } else if ((mesh->first == "surface") && surface_done) {
      // pass
    } else {
      // vis successful steps
      std::string plist_name = "visualization "+mesh->first;
      // in the case of just a domain mesh, we want to allow no name.
      if ((mesh->first == "domain") && !parameter_list_->isSublist(plist_name)) {
        plist_name = "visualization";
      }

      if (parameter_list_->isSublist(plist_name)) {
        Teuchos::ParameterList& vis_plist = parameter_list_->sublist(plist_name);
        Teuchos::RCP<Amanzi::Visualization> vis =
          Teuchos::rcp(new Amanzi::Visualization(vis_plist, comm_));
        vis->set_mesh(mesh->second.first);
        vis->CreateFiles();
        visualization_.push_back(vis);
      }

      // vis unsuccessful steps
      std::string fail_plist_name = "visualization "+mesh->first+" failed steps";
      // in the case of just a domain mesh, we want to allow no name.
      if ((mesh->first == "domain") && !parameter_list_->isSublist(fail_plist_name)) {
        fail_plist_name = "visualization failed steps";
      }

      if (parameter_list_->isSublist(fail_plist_name)) {
        Teuchos::ParameterList& fail_vis_plist = parameter_list_->sublist(fail_plist_name);
        Teuchos::RCP<Amanzi::Visualization> fail_vis =
          Teuchos::rcp(new Amanzi::Visualization(fail_vis_plist, comm_));
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
  for (std::vector<Teuchos::RCP<Amanzi::Visualization> >::iterator vis=visualization_.begin();
       vis!=visualization_.end(); ++vis) {
    (*vis)->RegisterWithTimeStepManager(tsm_.ptr());
  }

  // -- register checkpoint times
  checkpoint_->RegisterWithTimeStepManager(tsm_.ptr());

  // -- register observation times
  observations_->RegisterWithTimeStepManager(tsm_.ptr());

  // -- register the final time
  tsm_->RegisterTimeEvent(t1_);

  // -- register any intermediate requested times
  if (coordinator_list_->isSublist("required times")) {
    Teuchos::ParameterList& sublist = coordinator_list_->sublist("required times");
    Amanzi::IOEvent pause_times(sublist);
    pause_times.RegisterWithTimeStepManager(tsm_.ptr());
  }

  // Create an intermediate state that will store the updated solution until
  // we know it has succeeded.
  S_next_ = Teuchos::rcp(new Amanzi::State(*S_));
  *S_next_ = *S_;
  S_inter_ = Teuchos::rcp(new Amanzi::State(*S_));
  *S_inter_ = *S_;

  // set the states in the PKs
  //Teuchos::RCP<const State> cS = S_; // ensure PKs get const reference state
  pk_->set_states(S_, S_inter_, S_next_); // note this does not allow subcycling
}

void Coordinator::finalize() {
  // Force checkpoint at the end of simulation.
  // Only do if the checkpoint was not already written, or we would be writing
  // the same file twice.
  // This really should be removed, but for now is left to help stupid developers.
  if (!checkpoint_->DumpRequested(S_next_->cycle(), S_next_->time())) {
    pk_->CalculateDiagnostics(S_next_);
    WriteCheckpoint(checkpoint_.ptr(), S_next_.ptr(), 0.0);
  }

  // flush observations to make sure they are saved
  observations_->Flush();
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


void Coordinator::report_memory() {
  // report the memory high water mark (using ru_maxrss)
  // this should be called at the very end of a simulation
  if (vo_->os_OK(Teuchos::VERB_MEDIUM)) {
    double global_ncells(0.0);
    double local_ncells(0.0);
    for (Amanzi::State::mesh_iterator mesh = S_->mesh_begin(); mesh != S_->mesh_end(); ++mesh) {
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
  for (Amanzi::State::field_iterator field=S_->field_begin(); field!=S_->field_end(); ++field) {
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



void Coordinator::read_parameter_list() {
  t0_ = coordinator_list_->get<double>("start time");
  t1_ = coordinator_list_->get<double>("end time");
  std::string t0_units = coordinator_list_->get<std::string>("start time units", "s");
  std::string t1_units = coordinator_list_->get<std::string>("end time units", "s");

  if (t0_units == "s") {
    // internal units in s
  } else if (t0_units == "d") { // days
    t0_ = t0_ * 24.0*3600.0;
  } else if (t0_units == "yr") { // years
    t0_ = t0_ * 365.*24.0*3600.0;
  } else {
    Errors::Message message("Coordinator: error, invalid start time units");
    Exceptions::amanzi_throw(message);
  }

  if (t1_units == "s") {
    // internal units in s
  } else if (t1_units == "d") { // days
    t1_ = t1_ * 24.0*3600.0;
  } else if (t1_units == "yr") { // years
    t1_ = t1_ * 365.25*24.0*3600.0;
  } else {
    Errors::Message message("Coordinator: error, invalid end time units");
    Exceptions::amanzi_throw(message);
  }

  max_dt_ = coordinator_list_->get<double>("max time step size", 1.0e99);
  min_dt_ = coordinator_list_->get<double>("min time step size", 1.0e-12);
  cycle0_ = coordinator_list_->get<int>("start cycle",0);
  cycle1_ = coordinator_list_->get<int>("end cycle",-1);
  duration_ = coordinator_list_->get<double>("wallclock duration [hrs]", -1.0);

  // restart control
  restart_ = coordinator_list_->isParameter("restart from checkpoint file");
  if (restart_) {
    restart_filename_ = coordinator_list_->get<std::string>("restart from checkpoint file");
    // likely should ensure the file exists here? --etc
  }
}



// -----------------------------------------------------------------------------
// acquire the chosen timestep size
// -----------------------------------------------------------------------------
double Coordinator::get_dt(bool after_fail) {
  // get the physical step size
  double dt = pk_->get_dt();

  // check if the step size has gotten too small
  if (dt < min_dt_) {
    Errors::Message message("Coordinator: error, timestep too small");
    Exceptions::amanzi_throw(message);
  }

  // cap the max step size
  if (dt > max_dt_) {
    dt = max_dt_;
  }

  // ask the step manager if this step is ok
  dt = tsm_->TimeStep(S_next_->time(), dt, after_fail);
  return dt;
}


bool Coordinator::advance(double t_old, double t_new) {
  double dt = t_new - t_old;

  S_next_->advance_time(dt);
  bool fail = pk_->AdvanceStep(t_old, t_new, false);
  fail |= !pk_->ValidStep();

  // advance the iteration count and timestep size
  S_next_->advance_cycle();

  if (!fail) {
    // commit the state
    pk_->CommitStep(t_old, t_new, S_next_);
    
    // make observations, vis, and checkpoints
    observations_->MakeObservations(*S_next_);
    visualize();
    checkpoint(dt);

    // we're done with this time step, copy the state
    *S_ = *S_next_;
    *S_inter_ = *S_next_;

  } else {
    // Failed the timestep.  
    // Potentially write out failed timestep for debugging
    for (std::vector<Teuchos::RCP<Amanzi::Visualization> >::iterator vis=failed_visualization_.begin();
         vis!=failed_visualization_.end(); ++vis) {
      WriteVis((*vis).ptr(), S_next_.ptr());
    }

    // The timestep sizes have been updated, so copy back old soln and try again.
    *S_next_ = *S_;

    // check whether meshes are deformable, and if so, recover the old coordinates
    for (Amanzi::State::mesh_iterator mesh=S_->mesh_begin();
         mesh!=S_->mesh_end(); ++mesh) {
      if (S_->IsDeformableMesh(mesh->first)) {
        // collect the old coordinates
        std::string node_key = std::string("vertex_coordinate_")+mesh->first;
        Teuchos::RCP<const Amanzi::CompositeVector> vc_vec = S_->GetFieldData(node_key);
        vc_vec->ScatterMasterToGhosted();
        const Epetra_MultiVector& vc = *vc_vec->ViewComponent("node", true);
        std::vector<int> node_ids(vc.MyLength());
        Amanzi::AmanziGeometry::Point_List old_positions(vc.MyLength());
        for (int n=0;n!=vc.MyLength();++n) {
          node_ids[n] = n;
          if (mesh->second.first->space_dimension() == 2) {
            old_positions[n] = Amanzi::AmanziGeometry::Point(vc[0][n], vc[1][n]);
          } else {
            old_positions[n] = Amanzi::AmanziGeometry::Point(vc[0][n], vc[1][n], vc[2][n]);
          }
        }

        // undeform the mesh
        Amanzi::AmanziGeometry::Point_List final_positions;
        mesh->second.first->deform(node_ids, old_positions, false, &final_positions);
      }
    }
  }
  return fail;
}

void Coordinator::visualize(bool force) {
  // write visualization if requested
  bool dump = force;
  if (!dump) {
    for (std::vector<Teuchos::RCP<Amanzi::Visualization> >::iterator vis=visualization_.begin();
         vis!=visualization_.end(); ++vis) {
      if ((*vis)->DumpRequested(S_next_->cycle(), S_next_->time())) {
        dump = true;
      }
    }
  }

  if (dump) {
    pk_->CalculateDiagnostics(S_next_);
  }

  for (std::vector<Teuchos::RCP<Amanzi::Visualization> >::iterator vis=visualization_.begin();
       vis!=visualization_.end(); ++vis) {
    if (force || (*vis)->DumpRequested(S_next_->cycle(), S_next_->time())) {
      WriteVis((*vis).ptr(), S_next_.ptr());
    }
  }
}

void Coordinator::checkpoint(double dt, bool force) {
  if (force || checkpoint_->DumpRequested(S_next_->cycle(), S_next_->time())) {
    WriteCheckpoint(checkpoint_.ptr(), S_next_.ptr(), dt);
  }
}


// -----------------------------------------------------------------------------
// timestep loop
// -----------------------------------------------------------------------------
void Coordinator::cycle_driver() {
  // wallclock duration -- in seconds
  const double duration(duration_ * 3600);

  // start at time t = t0 and initialize the state.
  {
    Teuchos::TimeMonitor monitor(*setup_timer_);
    setup();   
    initialize();
   
  }

  //  exit(0);

  // get the intial timestep -- note, this would have to be fixed for a true restart
  double dt = get_dt(false);

  // visualization at IC
  visualize();
  checkpoint(dt);



  // iterate process kernels
  {
    Teuchos::TimeMonitor cycle_monitor(*cycle_timer_);
#if !DEBUG_MODE
  try {
#endif
    bool fail = false;
    while ((S_->time() < t1_) &&
           ((cycle1_ == -1) || (S_->cycle() <= cycle1_)) &&
           (duration_ < 0 || timer_->totalElapsedTime(true) < duration)) {
      if (vo_->os_OK(Teuchos::VERB_MEDIUM)) {
        Teuchos::OSTab tab = vo_->getOSTab();
        *vo_->os() << "======================================================================"
                  << std::endl << std::endl;
        *vo_->os() << "Cycle = " << S_->cycle();
        *vo_->os() << std::setprecision(15) << ",  Time [days] = "<< S_->time() / (60*60*24);
        *vo_->os() << ",  dt [days] = " << dt / (60*60*24)  << std::endl;
        *vo_->os() << "----------------------------------------------------------------------"
                  << std::endl;
      }

      *S_->GetScalarData("dt", "coordinator") = dt;
      *S_inter_->GetScalarData("dt", "coordinator") = dt;
      *S_next_->GetScalarData("dt", "coordinator") = dt;

      S_->set_initial_time(S_->time());
      S_->set_final_time(S_->time() + dt);
      S_->set_intermediate_time(S_->time());

      fail = advance(S_->time(), S_->time() + dt);
      //S_->WriteStatistics(vo_);  
      dt = get_dt(fail);

    } // while not finished


#if !DEBUG_MODE
  }

  catch (Amanzi::Exceptions::Amanzi_exception &e) {
    // write one more vis for help debugging
    S_next_->advance_cycle();
    visualize(true); // force vis

    // flush observations to make sure they are saved
    observations_->Flush();

    // catch errors to dump two checkpoints -- one as a "last good" checkpoint
    // and one as a "debugging data" checkpoint.
    checkpoint_->set_filebasename("last_good_checkpoint");
    WriteCheckpoint(checkpoint_.ptr(), S_.ptr(), dt);
    checkpoint_->set_filebasename("error_checkpoint");
    WriteCheckpoint(checkpoint_.ptr(), S_next_.ptr(), dt);
    throw e;
  }
#endif
  }


  // finalizing simulation                                                                                                                                                                                                               
  S_->WriteStatistics(vo_);  
  report_memory();
  Teuchos::TimeMonitor::summarize(*vo_->os());

  finalize();

} // cycle driver

} // close namespace Amanzi
