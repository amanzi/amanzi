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
  coordinator_plist_ = parameter_list_.sublist("Coordinator");
  read_parameter_list();

  // create the state object
  Teuchos::ParameterList state_plist = parameter_list_.sublist("State");
  S_ = Teuchos::rcp(new State(state_plist));
  S_->RegisterDomainMesh(mesh_);

  // checkpointing for the state
  Teuchos::ParameterList chkp_plist = parameter_list_.sublist("Checkpoint");
  checkpoint_ = Teuchos::rcp(new Checkpoint(chkp_plist, comm_));

  // create the top level PK
  Teuchos::ParameterList pks_list = parameter_list_.sublist("PKs");
  Teuchos::ParameterList::ConstIterator pk_item = pks_list.begin();

  const std::string &pk_name = pks_list.name(pk_item);
  const Teuchos::ParameterEntry &pk_value = pks_list.entry(pk_item);

  // -- create the solution
  std::cout << "Coordinator creating TreeVec with name: " << pk_name << std::endl;
  soln_ = Teuchos::rcp(new TreeVector(pk_name));

  // -- create the pk
  PKFactory pk_factory;
  std::cout << "Coordinator creating PK with name: " << pk_name << std::endl;
  pk_ = pk_factory.CreatePK(pks_list.sublist(pk_name), S_, soln_);
  pk_->set_name(pk_name);

  // // create the observations
  // Teuchos::ParameterList observation_plist = parameter_list_.sublist("Observation");
  // observations_ = Teuchos::rcp(new UnstructuredObservations(observation_plist,
  //         output_observations_));
}

void Coordinator::initialize() {
  // initialize the state (which should initialize all independent variables)
  S_->Initialize();

  // initialize the process kernels (which should initialize all dependent variables)
  pk_->initialize(S_);

  // commit state to get secondary variables intiailized
  pk_->commit_state(0.0, S_);
  pk_->calculate_diagnostics(S_);

  // Check that all fields have now been initialized or die.
  S_->CheckAllInitialized();

  // vis for the state
  // HACK to vis with a surrogate surface mesh.  This needs serious re-design. --etc
  bool surface_done = false;
  if (S_->has_mesh("surface") && S_->has_mesh("surface_3d")) {
    Teuchos::RCP<const AmanziMesh::Mesh> surface_3d = S_->Mesh("surface_3d");
    Teuchos::RCP<const AmanziMesh::Mesh> surface = S_->Mesh("surface");

    std::string plist_name = "Visualization surface";
    Teuchos::ParameterList& vis_plist = parameter_list_.sublist(plist_name);
      vis_plist.set<std::string>("File Name Base","surface");
    Teuchos::RCP<Visualization> vis =
      Teuchos::rcp(new Visualization(vis_plist, comm_));
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
      std::string plist_name = "Visualization "+mesh->first;
      Teuchos::ParameterList& vis_plist = parameter_list_.sublist(plist_name);
      vis_plist.set<std::string>("File Name Base",mesh->first);
      Teuchos::RCP<Visualization> vis =
        Teuchos::rcp(new Visualization(vis_plist, comm_));
      vis->set_mesh(mesh->second);
      vis->CreateFiles();
      visualization_.push_back(vis);
    }
  }
  parameter_list_.print(std::cout);
}

void Coordinator::read_parameter_list() {
  t0_ = coordinator_plist_.get<double>("Start Time");
  t1_ = coordinator_plist_.get<double>("End Time");
  string t0_units = coordinator_plist_.get<string>("Start Time Units", "s");
  string t1_units = coordinator_plist_.get<string>("End Time Units", "s");

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

  max_dt_ = coordinator_plist_.get<double>("Max Time Step Size", 1.0e99);
  min_dt_ = coordinator_plist_.get<double>("Min Time Step Size", 1.0e-12);
  end_cycle_ = coordinator_plist_.get<int>("End Cycle",-1);
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
  TimeStepManager tsm;
  // register times with the tsm
  // -- register visualization times
  for (std::vector<Teuchos::RCP<Visualization> >::iterator vis=visualization_.begin();
       vis!=visualization_.end(); ++vis) {
    (*vis)->RegisterWithTimeStepManager(tsm);
  }

  // -- register observation times
  //if (observations_) observations_->register_with_time_step_manager(TSM);
  // -- register the final time
  tsm.RegisterTimeEvent(t1_);

  // make observations
  //  observations_->MakeObservations(*S_);

  // write visualization if requested at IC
  for (std::vector<Teuchos::RCP<Visualization> >::iterator vis=visualization_.begin();
       vis!=visualization_.end(); ++vis) {
    S_->WriteVis((*vis).ptr());
  }

  // we need to create an intermediate state that will store the updated
  // solution until we know it has succeeded
  S_next_ = Teuchos::rcp(new State(*S_));
  *S_next_ = *S_;

  // set the states in the PKs
  Teuchos::RCP<const State> cS = S_; // ensure PKs get const reference state
  pk_->set_states(cS, S_, S_next_); // note this does not allow subcycling

  // iterate process kernels
  double dt;
  bool fail = false;
  while ((S_->time() < t1_) && ((end_cycle_ == -1) || (S_->cycle() <= end_cycle_))) {
    // get the physical step size
    dt = pk_->get_dt();
    std::cout << "TIME STPPR: Got PK recommendation of " << dt << std::endl;

    // check if the step size has gotten too small
    if (dt < min_dt_) {
      Errors::Message message("Coordinator: error, timestep too small");
      Exceptions::amanzi_throw(message);
    }

    // cap the max step size
    if (dt > max_dt_) {
      std::cout << "TIME STPPR: Capping step size at max of " << max_dt_ << std::endl;
      dt = max_dt_;
    }

    // ask the step manager if this step is ok
    dt = tsm.TimeStep(S_->time(), dt);
    std::cout << "TIME STPPR Manager: dt is " << dt << std::endl;

    std::cout << "Cycle = " << S_->cycle();
    std::cout << ",  Time [days] = "<< S_->time() / (60*60*24);
    std::cout << ",  dt [days] = " << dt / (60*60*24)  << std::endl;

    // advance
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
          S_next_->WriteVis((*vis).ptr());
        }
      }

      S_next_->WriteCheckpoint(checkpoint_);

      // write restart dump if requested
      // restart->dump_state(*S_next_);

      // we're done with this time step, copy the state
      *S_ = *S_next_;
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
    S_next_->WriteVis((*vis).ptr());
  }

  S_next_->WriteCheckpoint(checkpoint_, true);

  // dump observations
  //  output_observations_.print(std::cout);
} // cycle driver

} // close namespace Amanzi
