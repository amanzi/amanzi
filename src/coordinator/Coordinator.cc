/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

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
  this->setLinePrefix("Amanzi::Coordinator         ");
  // make sure that the line prefix is printed
  this->getOStream()->setShowLinePrefix(true);

  // Read the sublist for verbosity settings.
  Teuchos::readVerboseObjectSublist(&parameter_list,this);

  using Teuchos::OSTab;
  Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
  Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
  OSTab tab = this->getOSTab(); // This sets the line prefix and adds one tab

  coordinator_plist_ = parameter_list.sublist("Coordinator");
  read_parameter_list();

  // create the state object
  Teuchos::ParameterList state_parameter_list = parameter_list.sublist("State");
  S = Teuchos::rcp(new State( state_parameter_list, mesh_maps));

  // create the PKs
  Teuchos::ParameterList pks_list = parameter_list.sublist("PKs");
  Teuchos::ParameterList::ConstIterator pk_item = pks_list.begin();

  const std::string &pk_name = pks_list.name(pk_item);
  const Teuchos::ParameterEntry &pk_value = pks_list.entry(pk_item);
  pk_ = pk_factory_.create_pk(pks_list.sublist(pk_name), S_, solution_);

  // create the observations
  Teuchos::ParameterList observation_plist = parameter_list.sublist("Observation");
  observations = new Amanzi::Unstructured_observations(observation_plist,
                                                       output_observations);

  // create the visualization object
  if (parameter_list.isSublist("Visualization Data")) {
      Teuchos::ParameterList vis_parameter_list =
        parameter_list.sublist("Visualization Data");
      visualization = new Amanzi::Vis(vis_parameter_list, comm);
      visualization->create_files(*mesh_maps);
  } else {
    visualization = new Amanzi::Vis();
  }
}

void Coordinator::initialize() {
  // initialize the state (which should initialize all independent variables)
  S_->initialize();

  // initialize the process kernels (which should initialize all dependent variables)
  pk_->initialize(S_, solution_);
}

void Coordinator::read_parameter_list() {
  T0 = coordinator_plist_.get<double>("Start Time");
  T1 = coordinator_plist_.get<double>("End Time");
  end_cycle = coordinator_plist_.get<int>("End Cycle",-1);
}

void Coordinator::cycle_driver () {
  using Teuchos::OSTab;
  Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
  Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
  OSTab tab = this->getOSTab(); // This sets the line prefix and adds one tab

  // start at time T=T0, iteration 0;
  S_->set_time(T0);
  int iter = 0;
  S_->set_cycle(iter);

  // initialize the state.  In a flow steady-state problem,
  // this should include advancing flow to steady state
  // (which should be done by flow_pk->initialize_state(S)
  initialize();
  S->set_time(T0); // in case steady state solve changed this

  // make observations
  observations->make_observations(*S_);

  // write visualization if requested at IC
  S_->write_vis(*visualization_);

  // we need to create an intermediate state that will store the 
  // updated solution until we know it has succeeded
  Teuchos::RCP<State> S_new = Teuchos::rcp(new State(*S_));

  // iterate process kernels
  double mpc_dT;
  bool fail = 0;
  while (  (S->get_time() <= T1)  &&   ((end_cycle == -1) || (iter <= end_cycle)) ) {
    mpc_dT = get_dT();

    if(out.get() && includesVerbLevel(verbLevel,Teuchos::VERB_LOW,true)) {
      *out << "Cycle = " << iter;
      *out << ",  Time = "<< S->get_time() / (60*60*24);
      *out << ",  dT = " << mpc_dT / (60*60*24)  << std::endl;
    }

    // advance
    fail = advance_transient(mpc_dT, S, S_new);
    if (fail) {
      Errors::Message message("Coordinator: error advancing time");
      Exceptions::amanzi_throw(message);
    } else {
      // update the time in the state object
      S->advance_time(mpc_dT);

      // we're done with this time step
      S = S_new;

      // advance the iteration count
      ++iter;
      S->set_cycle(iter);

      // make observations
      observations->make_observations(*S);

      // write visualization if requested
      S->write_vis(*visualization);

      // write restart dump if requested
      // restart->dump_state(*S);
    } // if fail
  } // while not finished

  // dump observations
  output_observations.print(std::cout);
} // cycle driver

} // close namespace Amanzi
