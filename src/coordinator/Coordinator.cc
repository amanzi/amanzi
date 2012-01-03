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
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"
#include "Epetra_Comm.h"
#include "Epetra_MpiComm.h"

// TODO: We are using depreciated parts of boost::filesystem
#define BOOST_FILESYSTEM_VERSION 2
#include "boost/filesystem/operations.hpp"
#include "boost/filesystem/path.hpp"

#include "Coordinator.hh"

namespace Amanzi {

Coordinator::Coordinator(Teuchos::ParameterList &parameter_list,
         Teuchos::RCP<Amanzi::AmanziMesh::Mesh> &mesh_maps,
         Epetra_MpiComm* comm,
         Amanzi::ObservationData& output_observations):
  parameter_list_(parameter_list),
  mesh_maps_(mesh_maps),
  comm_(comm),
  output_observations_(output_observations) {
  coordinator_init();
};

void Coordinator::coordinator_init() {
  // set the line prefix for output
  this->setLinePrefix("Amanzi::Coordinator  ");
  // make sure that the line prefix is printed
  this->getOStream()->setShowLinePrefix(true);

  // Read the sublist for verbosity settings.
  Teuchos::readVerboseObjectSublist(&parameter_list_,this);

  using Teuchos::OSTab;
  Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
  Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
  OSTab tab = this->getOSTab(); // This sets the line prefix and adds one tab

  coordinator_plist_ = parameter_list_.sublist("Coordinator");
  read_parameter_list();

  // create the state object
  Teuchos::ParameterList state_parameter_list = parameter_list_.sublist("State");
  S_ = Teuchos::rcp(new State( state_parameter_list, mesh_maps_));

  // create the top level PK
  Teuchos::ParameterList pks_list = parameter_list_.sublist("PKs");
  Teuchos::ParameterList::ConstIterator pk_item = pks_list.begin();

  const std::string &pk_name = pks_list.name(pk_item);
  const Teuchos::ParameterEntry &pk_value = pks_list.entry(pk_item);

  // -- create the solution
  if(out.get() && includesVerbLevel(verbLevel,Teuchos::VERB_MEDIUM,true)) {
    *out << "Coordinator creating TreeVec with name: " << pk_name << std::endl;
  }
  soln_ = Teuchos::rcp(new TreeVector(pk_name));

  // -- create the pk
  if(out.get() && includesVerbLevel(verbLevel,Teuchos::VERB_MEDIUM,true)) {
    *out << "Coordinator creating PK with name: " << pk_name << std::endl;
  }
  pk_ = pk_factory_.create_pk(pks_list.sublist(pk_name), S_, soln_);
  pk_->set_name(pk_name);

  // create the observations
  Teuchos::ParameterList observation_plist = parameter_list_.sublist("Observation");
  observations_ = new Amanzi::Unstructured_observations(observation_plist,
          output_observations_);

  // create the visualization object
  if (parameter_list_.isSublist("Visualization Data")) {
      Teuchos::ParameterList vis_parameter_list =
        parameter_list_.sublist("Visualization Data");
      visualization_ = new Amanzi::Vis(vis_parameter_list, comm_);
      visualization_->create_files(*mesh_maps_);
  } else {
    visualization_ = new Amanzi::Vis();
  }
}

void Coordinator::initialize() {
  // initialize the state (which should initialize all independent variables)
  S_->initialize();

  // initialize the process kernels (which should initialize all dependent variables)
  pk_->initialize(S_);

  // Check that all fields have now been initialized or die.
  S_->check_all_initialized();
}

void Coordinator::read_parameter_list() {
  T0_ = coordinator_plist_.get<double>("Start Time");
  T1_ = coordinator_plist_.get<double>("End Time");
  end_cycle_ = coordinator_plist_.get<int>("End Cycle",-1);
}

void Coordinator::cycle_driver () {
  using Teuchos::OSTab;
  Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
  Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
  OSTab tab = this->getOSTab(); // This sets the line prefix and adds one tab

  // start at time T=T0, iteration 0;
  S_->set_time(T0_);
  int iter = 0;
  S_->set_cycle(iter);

  // initialize the state.  In a flow steady-state problem,
  // this should include advancing flow to steady state
  // (which should be done by flow_pk->initialize_state(S)
  initialize();
  S_->set_time(T0_); // in case steady state solve changed this

  // make observations
  observations_->make_observations(*S_);

  // write visualization if requested at IC
  S_->write_vis(*visualization_);

  // we need to create an intermediate state that will store the updated
  // solution until we know it has succeeded
  S_next_ = Teuchos::rcp(new State(*S_));

  // set the states in the PKs
  Teuchos::RCP<const State> cS = S_; // cast as const as passing non-const to const
                                     // not possible under the RCP
  pk_->set_states(cS, S_next_);

  // iterate process kernels
  double mpc_dT;
  bool fail = 0;
  while (  (S_->get_time() <= T1_)  &&   ((end_cycle_ == -1) || (iter <= end_cycle_)) ) {
    mpc_dT = pk_->get_dt();

    if(out.get() && includesVerbLevel(verbLevel,Teuchos::VERB_LOW,true)) {
      *out << "Cycle = " << iter;
      *out << ",  Time = "<< S_->get_time() / (60*60*24);
      *out << ",  dT = " << mpc_dT / (60*60*24)  << std::endl;
    }

    // advance
    fail = pk_->advance(mpc_dT);
    if (fail) {
      Errors::Message message("Coordinator: error advancing time");
      Exceptions::amanzi_throw(message);
    } else {
      // update the time in the state object
      S_->advance_time(mpc_dT);

      // we're done with this time step, copy the state
      *S_ = *S_next_;

      // advance the iteration count
      ++iter;
      S_->set_cycle(iter);

      // make observations
      observations_->make_observations(*S_);

      // write visualization if requested
      S_->write_vis(*visualization_);

      // write restart dump if requested
      // restart->dump_state(*S);
    } // if fail
  } // while not finished

  // dump observations
  output_observations_.print(std::cout);
} // cycle driver

} // close namespace Amanzi
