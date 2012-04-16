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
  Teuchos::ParameterList state_parameter_list = parameter_list_.sublist("State");
  S_ = Teuchos::rcp(new State( state_parameter_list, mesh_));

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

//   // create the visualization object
//   if (parameter_list_.isSublist("Visualization Data")) {
//       Teuchos::ParameterList vis_parameter_list =
//         parameter_list_.sublist("Visualization Data");
//       visualization_ = new Amanzi::Vis(vis_parameter_list, comm_);
//       visualization_->create_files(*mesh_);
//   } else {
//     visualization_ = new Amanzi::Vis();
//   }
}

void Coordinator::initialize() {
  // initialize the state (which should initialize all independent variables)
  S_->Initialize();

  // initialize the process kernels (which should initialize all dependent variables)
  pk_->initialize(S_);

  // commit state to get secondary variables intiailized
  pk_->commit_state(0.0, S_);

  // Check that all fields have now been initialized or die.
  S_->CheckAllInitialized();
}

void Coordinator::read_parameter_list() {
  t0_ = coordinator_plist_.get<double>("Start Time");
  t1_ = coordinator_plist_.get<double>("End Time");
  max_dt_ = coordinator_plist_.get<double>("Max Time Step Size", 1.0e99);
  min_dt_ = coordinator_plist_.get<double>("Min Time Step Size", 1.0e-17);
  end_cycle_ = coordinator_plist_.get<int>("End Cycle",-1);
}

void Coordinator::cycle_driver () {

  // start at time t = t0 and initialize the state.  In a flow steady-state
  // problem, this should include advancing flow to steady state (which should
  // be done by flow_pk->initialize_state(S)
  S_->set_time(t0_);
  initialize();
  S_->set_time(t0_); // in case steady state solve changed this

  // make observations
  //  observations_->MakeObservations(*S_);

  // write visualization if requested at IC
  //  S_->WriteVis(*visualization_);

  // we need to create an intermediate state that will store the updated
  // solution until we know it has succeeded
  S_next_ = Teuchos::rcp(new State(*S_));

  // set the states in the PKs
  Teuchos::RCP<const State> cS = S_; // ensure PKs get const reference state
  pk_->set_states(cS, S_, S_next_);

  // iterate process kernels
  double dt;
  bool fail = false;
  int iter = 0;
  while ((S_->time() < t1_) && ((end_cycle_ == -1) || (iter <= end_cycle_))) {
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

    // check if we are within reach of a vis step, restart step, or end of
    // simulation, and adjust timestep appropriately
    if (S_->time() + dt > t1_) {
      dt = t1_ - S_->time();
    }

    std::cout << "Cycle = " << iter;
    std::cout << ",  Time = "<< S_->time() / (60*60*24);
    std::cout << ",  dt = " << dt / (60*60*24)  << std::endl;

    // advance
    fail = pk_->advance(dt);
    if (!fail) {
      // update the new state with the new solution
      pk_->solution_to_state(soln_, S_next_);
      pk_->commit_state(dt, S_next_);

      // advance the iteration count
      ++iter;

      // make observations
      //      observations_->MakeObservations(*S_next_);

      // write visualization if requested
      // S_
      // S_next_->WriteVis(*visualization_);

      // write restart dump if requested
      // restart->dump_state(*S_next_);

      // we're done with this time step, copy the state
      *S_ = *S_next_;
    } else {
      // Failed the timestep.  Do nothing, the timestep sizes have been
      // updated, so we can try again.
    }
  } // while not finished

  // dump observations
  //  output_observations_.print(std::cout);
} // cycle driver

} // close namespace Amanzi
