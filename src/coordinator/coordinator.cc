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
#include "errors.hh"

#include "global_verbosity.hh"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"

#include "time_step_manager.hh"
#include "visualization.hh"
#include "checkpoint.hh"
#include "state.hh"
#include "PK.hh"
#include "tree_vector.hh"
#include "pk_factory.hh"

#include "coordinator.hh"

#define DEBUG_MODE 1

namespace Amanzi {

Coordinator::Coordinator(Teuchos::ParameterList parameter_list,
                         Teuchos::RCP<State>& S,
                         Epetra_MpiComm* comm ) :
  parameter_list_(parameter_list),
  S_(S),
  comm_(comm),
  restart_(false) {
  coordinator_init();

  setLinePrefix("Coordinator");
  setDefaultVerbLevel(ATS::VerbosityLevel::level_);
  Teuchos::readVerboseObjectSublist(&parameter_list_,this);
  // get the fancy output ??
  verbosity_ = getVerbLevel();
  out_ = getOStream();

  // the time step manager coordinates all non-physical timesteps
  tsm_ = Teuchos::rcp(new TimeStepManager());
};

void Coordinator::coordinator_init() {
  coordinator_plist_ = parameter_list_.sublist("coordinator");
  read_parameter_list();

  // checkpointing for the state
  Teuchos::ParameterList chkp_plist = parameter_list_.sublist("checkpoint");
  checkpoint_ = Teuchos::rcp(new Checkpoint(chkp_plist, comm_));

  // create the top level PK
  Teuchos::ParameterList pks_list = parameter_list_.sublist("PKs");
  Teuchos::ParameterList::ConstIterator pk_item = pks_list.begin();
  const std::string &pk_name = pks_list.name(pk_item);

  // -- create the solution
  soln_ = Teuchos::rcp(new TreeVector(pk_name));

  // -- create the pk
  PKFactory pk_factory;
  Teuchos::ParameterList pk_list = pks_list.sublist(pk_name);
  pk_list.set("PK name", pk_name);
  pk_ = pk_factory.CreatePK(pk_list, soln_);
  pk_->setup(S_.ptr());

  // // create the observations
  // Teuchos::ParameterList observation_plist = parameter_list_.sublist("Observation");
  // observations_ = Teuchos::rcp(new UnstructuredObservations(observation_plist,
  //         output_observations_));
}

void Coordinator::initialize() {
  // Set up the state, creating all data structures.
  S_->Setup();

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

  // Initialize the process kernels (initializes all independent variables)
  pk_->initialize(S_.ptr());

  // Initialize the state (initializes all dependent variables).
  S_->Initialize();

  // commit the initial conditions.
  pk_->commit_state(0., S_);

  // Restart from checkpoint, part 2.
  if (restart_) {
    ReadCheckpoint(comm_, S_.ptr(), restart_filename_);
    t0_ = S_->time();
    cycle0_ = S_->cycle();
    
    DeformCheckpointMesh(S_.ptr());
  }

  // vis for the state
  // HACK to vis with a surrogate surface mesh.  This needs serious re-design. --etc
  bool surface_done = false;
  if (S_->HasMesh("surface") && S_->HasMesh("surface_3d")) {
    Teuchos::RCP<const AmanziMesh::Mesh> surface_3d = S_->GetMesh("surface_3d");
    Teuchos::RCP<const AmanziMesh::Mesh> surface = S_->GetMesh("surface");

    std::string plist_name = "visualization surface";
    Teuchos::ParameterList& vis_plist = parameter_list_.sublist(plist_name);
    Teuchos::RCP<Visualization> vis = Teuchos::rcp(new Visualization(vis_plist, comm_));
    vis->set_mesh(surface_3d);
    vis->CreateFiles();
    vis->set_mesh(surface);
    visualization_.push_back(vis);
    S_->RemoveMesh("surface_3d");
    surface_done = true;
  }

  for (State::mesh_iterator mesh=S_->mesh_begin();
       mesh!=S_->mesh_end(); ++mesh) {
    if (mesh->first == "surface_3d") {
      // pass
    } else if ((mesh->first == "surface") && surface_done) {
      // pass
    } else {
      std::string plist_name = "visualization "+mesh->first;
      // in the case of just a domain mesh, we want to allow no name.
      if ((mesh->first == "domain") && !parameter_list_.isSublist(plist_name)) {
        plist_name = "visualization";
      }

      Teuchos::ParameterList& vis_plist = parameter_list_.sublist(plist_name);
      Teuchos::RCP<Visualization> vis =
        Teuchos::rcp(new Visualization(vis_plist, comm_));
      vis->set_mesh(mesh->second);
      vis->CreateFiles();
      visualization_.push_back(vis);
    }
  }
}

void Coordinator::read_parameter_list() {
  t0_ = coordinator_plist_.get<double>("start time");
  t1_ = coordinator_plist_.get<double>("end time");
  string t0_units = coordinator_plist_.get<string>("start time units", "s");
  string t1_units = coordinator_plist_.get<string>("end time units", "s");

  if (t0_units == "s") {
    // internal units in s
  } else if (t0_units == "d") { // days
    t0_ = t0_ * 24.0*3600.0;
  } else if (t0_units == "yr") { // years
    t0_ = t0_ * 365.25*24.0*3600.0;
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

  max_dt_ = coordinator_plist_.get<double>("max time step size", 1.0e99);
  min_dt_ = coordinator_plist_.get<double>("min time step size", 1.0e-12);
  cycle0_ = coordinator_plist_.get<int>("start cycle",0);
  cycle1_ = coordinator_plist_.get<int>("end cycle",-1);

  // restart control
  restart_ = coordinator_plist_.isParameter("restart from checkpoint file");
  if (restart_) {
    restart_filename_ = coordinator_plist_.get<std::string>("restart from checkpoint file");
    // likely should ensure the file exists here? --etc
  }
}



// -----------------------------------------------------------------------------
// acquire the chosen timestep size
// -----------------------------------------------------------------------------
double Coordinator::get_dt() {
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
  dt = tsm_->TimeStep(S_next_->time(), dt);
  return dt;
}


// -----------------------------------------------------------------------------
// timestep loop
// -----------------------------------------------------------------------------
void Coordinator::cycle_driver() {
  // start at time t = t0 and initialize the state.  In a flow steady-state
  // problem, this should include advancing flow to steady state (which should
  // be done by flow_pk->initialize_state(S)
  S_->set_time(t0_);
  S_->set_cycle(cycle0_);
  initialize();
  S_->set_time(t0_); // in case steady state solve changed this
  S_->set_cycle(cycle0_);

  // register times with the tsm_
  // -- register visualization times
  for (std::vector<Teuchos::RCP<Visualization> >::iterator vis=visualization_.begin();
       vis!=visualization_.end(); ++vis) {
    (*vis)->RegisterWithTimeStepManager(tsm_.ptr());
  }

  // -- register observation times
  //if (observations_) observations_->register_with_time_step_manager(TSM);
  // -- register the final time
  tsm_->RegisterTimeEvent(t1_);

  // we need to create an intermediate state that will store the updated
  // solution until we know it has succeeded
  S_next_ = Teuchos::rcp(new State(*S_));
  *S_next_ = *S_;
  S_inter_ = Teuchos::rcp(new State(*S_));
  *S_inter_ = *S_;

  // set the states in the PKs
  //Teuchos::RCP<const State> cS = S_; // ensure PKs get const reference state
  pk_->set_states(S_, S_inter_, S_next_); // note this does not allow subcycling

  // get the intial timestep -- note, this would have to be fixed for a true restart
  double dt = get_dt();

  // make observations
  //  observations_->MakeObservations(*S_);

  // write visualization if requested at IC
  pk_->calculate_diagnostics(S_);
  for (std::vector<Teuchos::RCP<Visualization> >::iterator vis=visualization_.begin();
       vis!=visualization_.end(); ++vis) {
    WriteVis((*vis).ptr(), S_.ptr());
  }
  WriteCheckpoint(checkpoint_.ptr(), S_.ptr(), dt);

  // iterate process kernels
#if !DEBUG_MODE
  try {
#endif
    bool fail = false;
    while ((S_->time() < t1_) && ((cycle1_ == -1) || (S_->cycle() <= cycle1_))) {
      if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_MEDIUM, true)) {
        Teuchos::OSTab tab = getOSTab();
        *out_ << "======================================================================"
                  << std::endl << std::endl;
        *out_ << "Cycle = " << S_->cycle();
        *out_ << ",  Time [days] = "<< S_->time() / (60*60*24);
        *out_ << ",  dt [days] = " << dt / (60*60*24)  << std::endl;
        *out_ << "----------------------------------------------------------------------"
                  << std::endl;
      }

      S_next_->advance_time(dt);
      fail = pk_->advance(dt);

      // advance the iteration count and timestep size
      S_next_->advance_cycle();
      dt = get_dt();

      if (!fail) {
        // make observations
        //      observations_->MakeObservations(*S_next_);

        // write visualization if requested
        // this needs to be fixed...
        bool dump = false;
        for (std::vector<Teuchos::RCP<Visualization> >::iterator vis=visualization_.begin();
             vis!=visualization_.end(); ++vis) {
          if ((*vis)->DumpRequested(S_next_->cycle(), S_next_->time())) {
            dump = true;
          }
        }
        if (dump) {
          pk_->calculate_diagnostics(S_next_);
        }

        for (std::vector<Teuchos::RCP<Visualization> >::iterator vis=visualization_.begin();
             vis!=visualization_.end(); ++vis) {
          if ((*vis)->DumpRequested(S_next_->cycle(), S_next_->time())) {
            WriteVis((*vis).ptr(), S_next_.ptr());
          }
        }

        if (checkpoint_->DumpRequested(S_next_->cycle(), S_next_->time())) {
          WriteCheckpoint(checkpoint_.ptr(), S_next_.ptr(), dt);
        }

        // write restart dump if requested
        // restart->dump_state(*S_next_);

        // we're done with this time step, copy the state
        *S_ = *S_next_;
        *S_inter_ = *S_next_;

      } else {
        // Failed the timestep.  The timestep sizes have been updated, so we can
        // try again.
        *S_next_ = *S_;
      }
    } // while not finished


#if !DEBUG_MODE
  }

  catch (Exceptions::Amanzi_exception &e) {
    // catch errors to dump two checkpoints -- one as a "last good" checkpoint
    // and one as a "debugging data" checkpoint.
    checkpoint_->set_filebasename("last_good_checkpoint");
    WriteCheckpoint(checkpoint_.ptr(), S_.ptr(), dt);
    checkpoint_->set_filebasename("error_checkpoint");
    WriteCheckpoint(checkpoint_.ptr(), S_next_.ptr(), dt);
    throw e;
  }
#endif

  // Force checkpoint at the end of simulation.
  // Only do if the checkpoint was not already written, or we would be writing
  // the same file twice.
  // This really should be removed, but for now is left to help stupid developers.
  if (!checkpoint_->DumpRequested(S_next_->cycle(), S_next_->time())) {
    pk_->calculate_diagnostics(S_next_);
    WriteCheckpoint(checkpoint_.ptr(), S_next_.ptr(), dt);
  }

  // dump observations
  //  output_observations_.print(std::cout);
} // cycle driver

} // close namespace Amanzi
