#include <utility>

#include "errors.hh"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"
#include "Epetra_Comm.h"
#include "Epetra_MpiComm.h"
#include "MPC.hpp"
#include "State.hpp"
#include "chemistry_state.hh"
#include "chemistry_pk.hh"
#include "Flow_State.hpp"
#include "Darcy_PK.hpp"
#include "Richards_PK.hpp"
#include "Transport_State.hpp"
#include "Transport_PK.hpp"
// TODO: We are using depreciated parts of boost::filesystem
#define BOOST_FILESYSTEM_VERSION 2
#include "boost/filesystem/operations.hpp"
#include "boost/filesystem/path.hpp"
#include "time_step_manager.hh"

namespace Amanzi {

using amanzi::chemistry::Chemistry_State;
using amanzi::chemistry::Chemistry_PK;
using amanzi::chemistry::ChemistryException;


/* *******************************************************************/
MPC::MPC(Teuchos::ParameterList parameter_list_,
         Teuchos::RCP<Amanzi::AmanziMesh::Mesh> mesh_maps_,
         Epetra_MpiComm* comm_,
         Amanzi::ObservationData& output_observations_):
    parameter_list(parameter_list_),
    mesh_maps(mesh_maps_),
    chemistry_enabled(false),
    comm(comm_),
    output_observations(output_observations_),
    transport_subcycling(0)
{
  mpc_init();
}


/* *******************************************************************/
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

  // let users selectively disable individual process kernels
  // to allow for testing of the process kernels separately
  transport_enabled =
      (mpc_parameter_list.get<string>("disable Transport_PK","no") == "no");

  if (mpc_parameter_list.get<string>("Chemistry Model","Off") != "Off") {
    chemistry_enabled = true;
  }

  flow_enabled =
      (mpc_parameter_list.get<string>("disable Flow_PK","no") == "no");

  if(out.get() && includesVerbLevel(verbLevel,Teuchos::VERB_LOW,true)) {
    *out << "The following process kernels are enabled: ";
    if (flow_enabled) *out << "Flow ";
    if (transport_enabled) *out << "Transport ";
    if (chemistry_enabled) *out << "Chemistry ";
    *out << std::endl;
  }

  if (transport_enabled || flow_enabled || chemistry_enabled) {
    Teuchos::ParameterList state_parameter_list = parameter_list.sublist("State");
    S = Teuchos::rcp(new State(state_parameter_list, mesh_maps));
  }

  //  
  // create auxilary state objects for the process models
  //

  // chemistry...
  if (chemistry_enabled) {
    if (parameter_list.isSublist("Chemistry")) {
      try {
        CS = Teuchos::rcp( new Chemistry_State( S ) );

        Teuchos::ParameterList chemistry_parameter_list =
            parameter_list.sublist("Chemistry");
 
        CPK = Teuchos::rcp( new Chemistry_PK(chemistry_parameter_list, CS) );
      } catch (const ChemistryException& chem_error) {
        *out << "MPC::mpc_init() : The chemistry process kernel constructor returned an error: "
                  << std::endl << chem_error.what() << std::endl;
        amanzi_throw(chem_error);
      }
    } else {
      Errors::Message message("MPC::mpc_init() : XML input file must contain a \'Chemistry\' section if the chemistry process kernel is enabled.");
      Exceptions::amanzi_throw(message);
    }
  }

  // transport...
  if (transport_enabled) {
    TS = Teuchos::rcp(new AmanziTransport::Transport_State(*S));
    Teuchos::ParameterList transport_parameter_list = parameter_list.sublist("Transport");

    bool subcycling = parameter_list.sublist("MPC").get<bool>("transport subcycling", false);
    transport_subcycling = (subcycling) ? 1 : 0;

    TPK = Teuchos::rcp(new AmanziTransport::Transport_PK(transport_parameter_list, TS));
    TPK->InitPK();
  }

  // flow...
  if (flow_enabled) {
    FS = Teuchos::rcp(new AmanziFlow::Flow_State(S));

    flow_model = mpc_parameter_list.get<string>("Flow model", "Darcy");
    if (flow_model == "Darcy") {
      FPK = Teuchos::rcp(new AmanziFlow::Darcy_PK(parameter_list, FS));
    } else if (flow_model == "Richards") {
      FPK = Teuchos::rcp(new AmanziFlow::Richards_PK(parameter_list, FS));
    } else if (flow_model == "Steady State Richards") {
      FPK = Teuchos::rcp(new AmanziFlow::Richards_PK(parameter_list, FS));
    } else {
      cout << "MPC: unknown flow model: " << flow_model << endl;
      throw std::exception();
    }
    FPK->InitPK();
  }
  if (flow_model == "Steady State Richards") {
    *out << "Flow will be off during the transient phase" << std::endl;
  }
  // done creating auxilary state objects and  process models

  // create the observations
  if (parameter_list.isSublist("Observation Data")) {
    Teuchos::ParameterList observation_plist = parameter_list.sublist("Observation Data"); 
    observations = new Amanzi::Unstructured_observations(observation_plist, output_observations);
  } else {
    observations = NULL;
  }

  // create the visualization object
  if (parameter_list.isSublist("Visualization Data"))  {
    Teuchos::ParameterList vis_parameter_list = parameter_list.sublist("Visualization Data");
    visualization = new Amanzi::Vis(vis_parameter_list, comm);
    visualization->create_files(*mesh_maps);
  } else {  // create a dummy vis object
    visualization = new Amanzi::Vis();
  }
  

  // create the restart object
  if (parameter_list.isSublist("Checkpoint Data")) {
    Teuchos::ParameterList checkpoint_parameter_list = parameter_list.sublist("Checkpoint Data");
    restart = new Amanzi::Restart(checkpoint_parameter_list, comm);
  } else {
    restart = new Amanzi::Restart();
  }

  // are we restarting from a file?
  // first assume we're not
  restart_requested = false;

  // then check if indeed we are
  if (mpc_parameter_list.isSublist("Restart from Checkpoint Data File")) {
    restart_requested = true;

    Teuchos::ParameterList& restart_parameter_list =
        mpc_parameter_list.sublist("Restart from Checkpoint Data File");

    restart_from_filename = restart_parameter_list.get<string>("Checkpoint Data File Name");
  }
}


/* *******************************************************************/
void MPC::read_parameter_list()  {
  end_cycle = mpc_parameter_list.get<int>("End Cycle",-1);
  
  Teuchos::ParameterList& ti_list =  mpc_parameter_list.sublist("Time Integration Mode");
  if (ti_list.isSublist("Initialize To Steady")) {
    ti_mode = INIT_TO_STEADY;
    
    Teuchos::ParameterList& init_to_steady_list = ti_list.sublist("Initialize To Steady");

    T0 = init_to_steady_list.get<double>("Start");
    Tswitch = init_to_steady_list.get<double>("Switch");
    T1 = init_to_steady_list.get<double>("End");

    dTsteady = init_to_steady_list.get<double>("Steady Initial Time Step");
    dTtransient = init_to_steady_list.get<double>("Transient Initial Time Step");
  } else if ( ti_list.isSublist("Steady")) {
    ti_mode = STEADY;
    
    Teuchos::ParameterList& steady_list = ti_list.sublist("Steady");

    T0 = steady_list.get<double>("Start");
    T1 = steady_list.get<double>("End");
    dTsteady = steady_list.get<double>("Initial Time Step");
  } else if ( ti_list.isSublist("Transient") ) {
    ti_mode = TRANSIENT;

    Teuchos::ParameterList& transient_list = ti_list.sublist("Transient");

    T0 = transient_list.get<double>("Start");
    T1 = transient_list.get<double>("End");
    dTtransient =  transient_list.get<double>("Initial Time Step");

  } else {
    Errors::Message message("MPC: no valid Time Integration Mode was specified, you must specify exactly one of Initialize To Steady, Steady, or Transient.");
    Exceptions::amanzi_throw(message);    
  }

  if (mpc_parameter_list.isSublist("Time Period Control")) {
    Teuchos::ParameterList& tpc_list =  mpc_parameter_list.sublist("Time Period Control"); 

    reset_times_    = tpc_list.get<Teuchos::Array<double> >("Start Times");
    reset_times_dt_ = tpc_list.get<Teuchos::Array<double> >("Initial Time Step");
    
    if (reset_times_.size() != reset_times_dt_.size()) {
      Errors::Message message("You must specify the same number of Reset Times and Initial Time Steps under Time Period Control");
      Exceptions::amanzi_throw(message);
    }    
  }
}


/* *******************************************************************/
void MPC::cycle_driver() {

  // create the time step manager
  Amanzi::TimeStepManager TSM;
  // register visualization times with the time step manager
  visualization->register_with_time_step_manager(TSM);
  // register observation times with the time step manager
  if (observations) observations->register_with_time_step_manager(TSM);
  // register reset_times
  TSM.RegisterTimeEvent(reset_times_.toVector());
  // if this is an init to steady run, register the switchover time
  if (ti_mode == INIT_TO_STEADY) TSM.RegisterTimeEvent(Tswitch);
  // register the final time
  TSM.RegisterTimeEvent(T1);



  enum time_step_limiter_type {FLOW_LIMITS, TRANSPORT_LIMITS, CHEMISTRY_LIMITS, MPC_LIMITS};
  time_step_limiter_type tslimiter;

  using Teuchos::OSTab;
  Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
  Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
  OSTab tab = this->getOSTab(); // This sets the line prefix and adds one tab

  //TSM.print(*out,0.0, 20.0); *out << std::endl;

  if (transport_enabled || flow_enabled || chemistry_enabled) {
    S->set_time(T0);  // start at time T=T0;
  }

  if (chemistry_enabled) {
    try {
      // these are the vectors that chemistry will populate with
      // the names for the auxillary output vectors and the
      // names of components
      std::vector<string> compnames;

      // total view needs this to be outside the constructor
      CPK->InitializeChemistry();
      CPK->set_chemistry_output_names(&auxnames);
      CPK->set_component_names(&compnames);

      // set the names in the visualization object
      S->set_compnames(compnames);

    } catch (const ChemistryException& chem_error) {
      std::cout << "MPC: Chemistry_PK.InitializeChemistry returned an error "
                << std::endl << chem_error.what() << std::endl;
      Exceptions::amanzi_throw(chem_error);
    }
  }

  int iter = 0;  // set the iteration counter to zero
  S->set_cycle(iter);

  // read the checkpoint file as requested
  if (restart_requested == true) {
    // re-initialize the state object
    restart->read_state(*S, restart_from_filename);
    iter = S->get_cycle();
    
    if (!reset_times_.empty()) {
      while (reset_times_.front()<S->get_time()) {
	reset_times_.erase(reset_times_.begin());
	reset_times_dt_.erase(reset_times_dt_.begin());
      }
    }
  } else { // no restart, we will call the PKs to allow them to init their auxilary data
    FPK->InitializeAuxiliaryData();
  }



  // write visualization output as requested
  if (chemistry_enabled) {
    // get the auxillary data from chemistry
    Teuchos::RCP<Epetra_MultiVector> aux = CPK->get_extra_chemistry_output_data();
    // write visualization data for timestep
    S->write_vis(*visualization, aux, auxnames, true, true);
  } else {
    // always write the initial visualization dump
    S->write_vis(*visualization, false, true);
  }

  // write a restart dump if requested (determined in dump_state)
  restart->dump_state(*S);

  
  if (flow_enabled) {
    if (ti_mode == STEADY) {
      FPK->InitSteadyState(S->get_time(), dTsteady);
    } else if ( ti_mode == TRANSIENT && flow_model !=std::string("Steady State Richards")) {
      FPK->InitTransient(S->get_time(), dTtransient);
    } else if ( ti_mode == INIT_TO_STEADY && flow_model == std::string("Steady State Saturated")) {
      if (!restart_requested) {
	FPK->InitializeSteadySaturated();
	FPK->CommitStateForTransport(FS);
	FPK->CommitState(FS);
	S->advance_time(Tswitch-T0);
      }
    } else if ( ti_mode == INIT_TO_STEADY ) {
      if (S->get_time() < Tswitch) {
	FPK->InitSteadyState(S->get_time(), dTsteady);
      } else {
	if (flow_model !=std::string("Steady State Richards")) {
	  FPK->InitTransient(S->get_time(), dTtransient);
	}
      }
    }
  }

  if (flow_enabled || transport_enabled || chemistry_enabled) {
    if (observations) observations->make_observations(*S);

    // we need to create an EpetraMulitVector that will store the
    // intermediate value for the total component concentration
    total_component_concentration_star =
        Teuchos::rcp(new Epetra_MultiVector(*S->get_total_component_concentration()));

    // then start time stepping
    while ((S->get_time() < T1) && ((end_cycle == -1) || (iter <= end_cycle))) {
      // determine the time step we are now going to take
      double chemistry_dT = 1e+99, transport_dT = 1e+99, flow_dT = 1e+99;
      double mpc_dT = 1e+99, limiter_dT = 1e+99, observation_dT = 1e+99;
      
      // Update our reset times (delete the next one if we just did it)
      if (!reset_times_.empty()) {
	if (S->get_last_time()>=reset_times_.front()) {
          reset_times_.erase(reset_times_.begin());
	  reset_times_dt_.erase(reset_times_dt_.begin());
	}
      }

      // catch the switchover time to transient
      if (flow_enabled) {
	if (ti_mode == INIT_TO_STEADY && S->get_last_time() < Tswitch && S->get_time() >= Tswitch) {
	  if(out.get() && includesVerbLevel(verbLevel,Teuchos::VERB_LOW,true)) {
	    *out << "Steady state computation complete... now running in transient mode." << std::endl;
	  }
	  // only init the transient problem if we need to 
	  if (flow_model != std::string("Steady State Richards")) {
	    FPK->InitTransient(S->get_time(), dTtransient);
	  }
	}
      }
      
      // find the flow time step
      if (flow_enabled) {
	// only if we are actually running with flow

	if ((ti_mode == STEADY) || 
	    (ti_mode == TRANSIENT && flow_model != std::string("Steady State Richards")) ||
	    (ti_mode == INIT_TO_STEADY && (flow_model != std::string("Steady State Richards") || S->get_time() < Tswitch)) ||
	    (ti_mode == INIT_TO_STEADY && (flow_model != std::string("Steady State Saturated")))) {
	  flow_dT = FPK->CalculateFlowDt();	  
	}
      }

      if (ti_mode == TRANSIENT || (ti_mode == INIT_TO_STEADY && S->get_time() >= Tswitch)) {
        if (transport_enabled) {
          double transport_dT_tmp = TPK->EstimateTransportDt();
          if (transport_subcycling == 0) transport_dT = transport_dT_tmp;
        }
        if (chemistry_enabled) {
          chemistry_dT = CPK->max_time_step();
        }
      }

      // take the mpc time step as the min of all suggested time steps 
      mpc_dT = std::min(std::min(flow_dT, transport_dT), chemistry_dT);

      // take the mpc time step as the min of the last limiter and itself 
      mpc_dT = TSM.TimeStep(S->get_time(), mpc_dT);

      // figure out who limits the time step
      if (mpc_dT == flow_dT) {
	tslimiter = FLOW_LIMITS;
      } else if (mpc_dT == transport_dT) {
	tslimiter = TRANSPORT_LIMITS;
      } else if (mpc_dT == chemistry_dT) {
	tslimiter = CHEMISTRY_LIMITS;
      } else {
	tslimiter = MPC_LIMITS;
      }
      
      // make sure we reset the timestep at switchover time
      if (ti_mode == INIT_TO_STEADY && S->get_time() >= Tswitch && S->get_last_time() < Tswitch) {
	// mpc_dT = std::min( mpc_dT, dTtransient );
	tslimiter = MPC_LIMITS; 
      }

      // make sure that if we are currently on a reset time, to reset the time step
      if (! ti_mode == STEADY) {
        if (!reset_times_.empty()) {
          // this is probably iffy...
          if (S->get_time() == reset_times_.front()) {
            *out << "Resetting the time integrator at time = " << S->get_time() << std::endl;
            mpc_dT = reset_times_dt_.front();
	    mpc_dT = TSM.TimeStep(S->get_time(), mpc_dT);
            tslimiter = MPC_LIMITS;
	    // now reset the flow time integrator..
	    FPK->InitTransient(S->get_time(), mpc_dT);
	  }
        }
      }

      // steady flow is special, it might redo a time step, so we print
      // time step info after we've advanced steady flow
      // first advance flow
      if (flow_enabled) {
	if ((ti_mode == STEADY) ||
	    (ti_mode == TRANSIENT && flow_model != std::string("Steady State Richards")) ||
	    (ti_mode == INIT_TO_STEADY && (flow_model != std::string("Steady State Richards") || S->get_time() < Tswitch)) ||
	    (ti_mode == INIT_TO_STEADY && (flow_model != std::string("Steady State Saturated")))) {
	  bool redo(false);
	  do {
	    redo = false;
	    try {
	      FPK->Advance(mpc_dT);
	    } 
	    catch (int itr) {
	      mpc_dT = 0.5*mpc_dT;
	      redo = true;
	      tslimiter = FLOW_LIMITS;
	      *out << "will repeat time step with smaller dT = " << mpc_dT << std::endl;
	    }
	  } while (redo);
	  FPK->CommitStateForTransport(FS);
	}
      }

      // write some info about the time step we are about to take
      // first determine what we will write about the time step limiter
      std::string limitstring("");
      switch (tslimiter) {
      case(MPC_LIMITS): 
	limitstring = std::string("(mpc limits timestep)");
	break;
      case (TRANSPORT_LIMITS): 
	limitstring = std::string("(transport limits timestep)");
	break;
      case (CHEMISTRY_LIMITS): 
	limitstring = std::string("(chemistry limits timestep)");
	break;
      case (FLOW_LIMITS): 
	break;
      }

      if(out.get() && includesVerbLevel(verbLevel,Teuchos::VERB_LOW,true)) {
	*out << std::setprecision(6);
        *out << "Cycle = " << iter;
        *out << ",  Time(years) = "<< S->get_time() / (365.25*60*60*24);
        *out << ",  dT(years) = " << mpc_dT / (365.25*60*60*24);
	*out << " " << limitstring;
	//*out << " " << S->get_time() << " " << mpc_dT;
        *out << std::endl;
      }
      // ==============================================================

      // then advance transport and chemistry
      if (ti_mode == TRANSIENT || (ti_mode == INIT_TO_STEADY && S->get_time() >= Tswitch) ) {
        if (transport_enabled) {
          TPK->Advance(mpc_dT);
          if (TPK->get_transport_status() == AmanziTransport::TRANSPORT_STATE_COMPLETE) {
            // get the transport state and commit it to the state
            Teuchos::RCP<AmanziTransport::Transport_State> TS_next = TPK->transport_state_next();
            *total_component_concentration_star = *TS_next->total_component_concentration();
          } else {
            Errors::Message message("MPC: error... Transport_PK.advance returned an error status");
            Exceptions::amanzi_throw(message);
          }
        } else { // if we're not advancing transport we still need to prepare for chemistry
          *total_component_concentration_star = *S->get_total_component_concentration();
        }
      
        if (chemistry_enabled) {
          try {
            CPK->advance(mpc_dT, total_component_concentration_star);
            S->update_total_component_concentration(CPK->get_total_component_concentration());
          } catch (const ChemistryException& chem_error) {
            std::ostringstream error_message;
            error_message << "MPC: error... Chemistry_PK.advance returned an error status";
            error_message << chem_error.what();
            Errors::Message message(error_message.str());
            Exceptions::amanzi_throw(message);
          }
        } else {
          S->update_total_component_concentration(*total_component_concentration_star);
        }
      }
      
      // update the time in the state object
      S->advance_time(mpc_dT);
      if (FPK->flow_status() == AmanziFlow::FLOW_STATUS_STEADY_STATE_COMPLETE) S->set_time(Tswitch);

      // ===========================================================
      // we're done with this time step, commit the state
      // in the process kernels

      FPK->CommitState(FS);
      if (ti_mode == TRANSIENT || (ti_mode == INIT_TO_STEADY && S->get_time() >= Tswitch) ) {
        if (transport_enabled) TPK->CommitState(TS);
        if (chemistry_enabled) CPK->commit_state(CS, mpc_dT);
      }

      // advance the iteration count
      iter++;
      S->set_cycle(iter);

      // make observations
      if (observations) observations->make_observations(*S);

      // write visualization if requested
      bool force(false);
      if (abs(S->get_time() - T1) < 1e-7) { 
	force = true;
      }

      if (ti_mode == INIT_TO_STEADY) 
        if (abs(S->get_time() - Tswitch) < 1e-7) {
          force = true;
        }

      if (chemistry_enabled) {
        // get the auxillary data
        Teuchos::RCP<Epetra_MultiVector> aux = CPK->get_extra_chemistry_output_data();
        
        // write visualization data for timestep if requested
        S->write_vis(*visualization, aux, auxnames, true, force);
      } else {
        S->write_vis(*visualization, false, force);
      }

      // write restart dump if requested
      restart->dump_state(*S, force);
    }
  }

  // some final output
  if (out.get() && includesVerbLevel(verbLevel,Teuchos::VERB_LOW,true))
  {
    *out << "Cycle = " << iter;
    *out << ",  Time(years) = "<< S->get_time()/ (365.25*60*60*24);
    *out << std::endl;
  }
}


/* *******************************************************************/
double MPC::time_step_limiter (double T, double dT, double T_end) {
  double time_remaining = T_end - T;

  if (time_remaining < 0.0) {
    Errors::Message message("MPC: time step limiter logic error, T_end must be greater than T.");
    Exceptions::amanzi_throw(message);
  }

  if (dT >= time_remaining) {
    return time_remaining;
  } else if ( dT > 0.75*time_remaining ) {
    return 0.5*time_remaining;
  } else {
    return dT;
  }
}

}  // namespace Amanzi

