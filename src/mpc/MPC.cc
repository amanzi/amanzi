/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

#include <string>

#include "errors.hh"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"
#include "Epetra_Comm.h"
#include "Epetra_MpiComm.h"

#include "MPC.hh"
#include "State.hh"
// TODO: We are using depreciated parts of boost::filesystem
#define BOOST_FILESYSTEM_VERSION 2
#include "boost/filesystem/operations.hpp"
#include "boost/filesystem/path.hpp"


namespace Amanzi
{

MPC::MPC(Teuchos::ParameterList &parameter_list_,
         Teuchos::RCP<Amanzi::AmanziMesh::Mesh> &mesh_maps_,
         Epetra_MpiComm* comm_,
         Amanzi::ObservationData& output_observations_):
  parameter_list(parameter_list_),
  mesh_maps(mesh_maps_),
  comm(comm_),
  output_observations(output_observations_) {
  mpc_init();
  }


void MPC::mpc_init() {
  // set the line prefix for output
  this->setLinePrefix("Amanzi::MPC         ");
  // make sure that the line prefix is printed
  this->getOStream()->setShowLinePrefix(true);

  // Read the sublist for verbosity settings.
  Teuchos::readVerboseObjectSublist(&parameter_list,this);

  using Teuchos::OSTab;
  Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
  Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
  OSTab tab = this->getOSTab(); // This sets the line prefix and adds one tab

  mpc_parameter_list =  parameter_list.sublist("MPC");
  read_parameter_list();

  // create the state object
  Teuchos::ParameterList state_parameter_list = parameter_list.sublist("State");
  S = Teuchos::rcp(new State( state_parameter_list, mesh_maps));

  // create the PKs
  Teuchos::ParameterList pks_list = parameter_list.sublist("PKs");
  for (Teuchos::ParameterList::ConstIterator i = pks_list.begin();
       i != pks_list.end(); ++i) {

    const std::string &name_i  = pks_list.name(i);
    const Teuchos::ParameterEntry  &entry_i = pks_list.entry(i);

    if (entry_i.isList()) {
      pks.push_back(pk_factory.create_pk(pks_list.sublist(name_i), S));
    }
  }

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

double MPC::get_dT() {
  double dt = 1.e99;
  for (std::vector< Teuchos::RCP<PK> >::iterator pk = pks.begin(); pk != pks.end(); ++pk) {
    dt = std::min(dt, (*pk)->get_dT());
  }
  return dt;
}

bool MPC::advance_transient(double dt, Teuchos::RCP<State> &S0, Teuchos::RCP<State> &S1) {
  bool fail = 0;
  for (std::vector< Teuchos::RCP<PK> >::iterator pk = pks.begin(); pk != pks.end(); ++pk) {
    fail = (*pk)->advance_transient(dt, S0, S1);
    if (fail) {
      return fail;
    } else {
      (*pk)->commit_state(dt, S1);
    }
  }
  return fail;
}

void MPC::initialize(Teuchos::ParameterList& plist, Teuchos::RCP<State> &S) {
  // initialize the state (which should initialize all independent variables)
  S->initialize();

  // initialize the process kernels (which should initialize all dependent variables)
  //Teuchos::ParameterList pks_list = parameter_list.sublist("PKs");
  for (std::vector< Teuchos::RCP<PK> >::iterator pk = pks.begin(); pk != pks.end(); ++pk) {
    (*pk)->initialize(S);
  }
}

void MPC::commit_state(double dt, Teuchos::RCP<State> &S) {
  for (std::vector< Teuchos::RCP<PK> >::iterator pk = pks.begin(); pk != pks.end(); ++pk) {
    (*pk)->commit_state(dt, S);
  }
}

void MPC::read_parameter_list() {
  T0 = mpc_parameter_list.get<double>("Start Time");
  T1 = mpc_parameter_list.get<double>("End Time");
  end_cycle = mpc_parameter_list.get<int>("End Cycle",-1);
}


void MPC::cycle_driver () {
  using Teuchos::OSTab;
  Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
  Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
  OSTab tab = this->getOSTab(); // This sets the line prefix and adds one tab

  // start at time T=T0, iteration 0;
  S->set_time(T0);
  int iter = 0;
  S->set_cycle(iter);

  // initialize the state.  In a flow steady-state problem,
  // this should include advancing flow to steady state
  // (which should be done by flow_pk->initialize_state(S)
  initialize(parameter_list, S);
  S->set_time(T0); // in case steady state solve changed this

  // make observations
  observations->make_observations(*S);

  // write visualization if requested at IC
  S->write_vis(*visualization);

  // we need to create an intermediate state that will store the 
  // updated solution until we know it has succeeded
  Teuchos::RCP<State> S_new = Teuchos::rcp(new State(*S));

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
      Errors::Message message("MPC: error advancing time");
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
