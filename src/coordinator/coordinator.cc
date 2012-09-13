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

#include "coordinator.hh"

namespace Amanzi {

Coordinator::Coordinator(Teuchos::ParameterList parameter_list,
                         Teuchos::RCP<Amanzi::AmanziMesh::Mesh>& mesh,
                         Epetra_MpiComm* comm ) :
  parameter_list_(parameter_list),
  mesh_(mesh),
  comm_(comm) {
  coordinator_init();
};

void Coordinator::coordinator_init() {
  coordinator_plist_ = parameter_list_.sublist("coordinator");
  read_parameter_list();

  // create the state object
  Teuchos::ParameterList state_plist = parameter_list_.sublist("state");
  S_ = Teuchos::rcp(new State(state_plist));
  S_->RegisterDomainMesh(mesh_);

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
  pk_->setup(S_);

  // // create the observations
  // Teuchos::ParameterList observation_plist = parameter_list_.sublist("Observation");
  // observations_ = Teuchos::rcp(new UnstructuredObservations(observation_plist,
  //         output_observations_));
}

void Coordinator::initialize() {
  // Set up the state, creating all data structures.
  S_->Setup();

  // Initialize the process kernels (initializes all independent variables)
  pk_->initialize(S_);

  // Initialize the state (initializes all dependent variables).
  S_->Initialize();

  // commit the initial conditions.
  pk_->commit_state(0., S_);

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
      vis_plist.set<std::string>("file name base",mesh->first);
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
  end_cycle_ = coordinator_plist_.get<int>("end cycle",-1);
}


// -----------------------------------------------------------------------------
// timestep loop
// -----------------------------------------------------------------------------
void Coordinator::cycle_driver() {
  // start at time t = t0 and initialize the state.  In a flow steady-state
  // problem, this should include advancing flow to steady state (which should
  // be done by flow_pk->initialize_state(S)
  S_->set_time(t0_);
  S_->set_cycle(0);
  initialize();
  S_->set_time(t0_); // in case steady state solve changed this
  S_->set_cycle(0);

  // the time step manager coordinates all non-physical timesteps
  Teuchos::Ptr<TimeStepManager> tsm = Teuchos::ptr(new TimeStepManager());

  // register times with the tsm
  // -- register visualization times
  for (std::vector<Teuchos::RCP<Visualization> >::iterator vis=visualization_.begin();
       vis!=visualization_.end(); ++vis) {
    (*vis)->RegisterWithTimeStepManager(tsm);
  }

  // -- register observation times
  //if (observations_) observations_->register_with_time_step_manager(TSM);
  // -- register the final time
  tsm->RegisterTimeEvent(t1_);

  // make observations
  //  observations_->MakeObservations(*S_);

  // write visualization if requested at IC
  pk_->calculate_diagnostics(S_);
  for (std::vector<Teuchos::RCP<Visualization> >::iterator vis=visualization_.begin();
       vis!=visualization_.end(); ++vis) {
    WriteVis((*vis).ptr(), S_.ptr());
  }

  // we need to create an intermediate state that will store the updated
  // solution until we know it has succeeded
  S_next_ = Teuchos::rcp(new State(*S_));
  *S_next_ = *S_;
  S_inter_ = Teuchos::rcp(new State(*S_));
  *S_inter_ = *S_;

  // set the states in the PKs
  //Teuchos::RCP<const State> cS = S_; // ensure PKs get const reference state
  pk_->set_states(S_, S_inter_, S_next_); // note this does not allow subcycling

  // iterate process kernels
  double dt;
  bool fail = false;
  while ((S_->time() < t1_) && ((end_cycle_ == -1) || (S_->cycle() <= end_cycle_))) {
    // get the physical step size
    dt = pk_->get_dt();

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
    dt = tsm->TimeStep(S_->time(), dt);

    if (comm_->MyPID() == 0) {
      std::cout << "======================================================================"
                << std::endl << std::endl;
      std::cout << "Cycle = " << S_->cycle();
      std::cout << ",  Time [days] = "<< S_->time() / (60*60*24);
      std::cout << ",  dt [days] = " << dt / (60*60*24)  << std::endl;
      std::cout << "----------------------------------------------------------------------"
                << std::endl;
    }

    S_next_->advance_time(dt);
    fail = pk_->advance(dt);

    // advance the iteration count
    S_next_->advance_cycle();

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
        WriteCheckpoint(checkpoint_.ptr(), S_next_.ptr());
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

  // force visualization and checkpoint at the end of simulation
  // this needs to be fixed -- should not force, but ask if we want to checkpoint/vis at end
  S_next_->advance_cycle(); // hackery to make the vis stop whining
  pk_->calculate_diagnostics(S_next_);

  for (std::vector<Teuchos::RCP<Visualization> >::iterator vis=visualization_.begin();
       vis!=visualization_.end(); ++vis) {
    WriteVis((*vis).ptr(), S_next_.ptr());
  }

  WriteCheckpoint(checkpoint_.ptr(), S_next_.ptr());

  // dump observations
  //  output_observations_.print(std::cout);
} // cycle driver

} // close namespace Amanzi
