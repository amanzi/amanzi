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

#define BOOST_FILESYSTEM_NO_DEPRECATED
#include "boost/filesystem/operations.hpp"
#include "boost/filesystem/path.hpp"
#include "time_step_manager.hh"

// make sure that we use default parameters
// that are consistent with the input translator
#include "InputParserIS-defaults.hh"

#include "TimerManager.hh"

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
        std::ostringstream error_message;
        error_message << "MPC:mpc_init(): error... Chemistry_PK() returned an error status: ";
        error_message << chem_error.what();
        Errors::Message message(error_message.str());
        Exceptions::amanzi_throw(message);
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

  // transport and chemistry...
  chem_trans_dt_ratio = CHEM_TRANS_DT_RATIO;
  if (transport_enabled && chemistry_enabled) {
    chem_trans_dt_ratio = parameter_list.sublist("MPC").get<double>("max chemistry to transport timestep ratio",CHEM_TRANS_DT_RATIO);
  }

  // flow...
  if (flow_enabled) {
    FS = Teuchos::rcp(new AmanziFlow::Flow_State(S));

    flow_model = mpc_parameter_list.get<string>("Flow model", "Darcy");
    if (flow_model == "Darcy") {
      FPK = Teuchos::rcp(new AmanziFlow::Darcy_PK(parameter_list, FS));
    } else if (flow_model == "Steady State Saturated") {
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

    // make sure that the restart file actually exists, if not throw an error
    boost::filesystem::path restart_from_filename_path(restart_from_filename);
    if (!boost::filesystem::exists(restart_from_filename_path)) {
      Errors::Message message("MPC: the specified restart file does not exist or is not a regular file.");
      Exceptions::amanzi_throw(message);
    }
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

    do_picard_ = init_to_steady_list.get<bool>("Use Picard",false);
  } else if ( ti_list.isSublist("Steady")) {
    ti_mode = STEADY;

    Teuchos::ParameterList& steady_list = ti_list.sublist("Steady");

    T0 = steady_list.get<double>("Start");
    T1 = steady_list.get<double>("End");
    dTsteady = steady_list.get<double>("Initial Time Step");

    do_picard_ = steady_list.get<bool>("Use Picard",false);
  } else if ( ti_list.isSublist("Transient") ) {
    ti_mode = TRANSIENT;

    Teuchos::ParameterList& transient_list = ti_list.sublist("Transient");

    T0 = transient_list.get<double>("Start");
    T1 = transient_list.get<double>("End");
    dTtransient =  transient_list.get<double>("Initial Time Step");

    do_picard_ = false;
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

  // start timers
  Amanzi::timer_manager.add("Chemistry PK", Amanzi::Timer::ACCUMULATE);
  Amanzi::timer_manager.add("Flow PK", Amanzi::Timer::ACCUMULATE);
  Amanzi::timer_manager.add("Transport PK", Amanzi::Timer::ACCUMULATE);
  Amanzi::timer_manager.add("I/O", Amanzi::Timer::ACCUMULATE);

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
    S->set_intermediate_time(Tswitch);
  }

  if (chemistry_enabled) {
    Amanzi::timer_manager.start("Chemistry PK");
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
      std::ostringstream error_message;
      error_message << "MPC:mpc_init(): error... Chemistry_PK.InitializeChemistry returned an error status: ";
      error_message << chem_error.what();
      Errors::Message message(error_message.str());
      Exceptions::amanzi_throw(message);
    }
    Amanzi::timer_manager.stop("Chemistry PK");
  }


  if (chemistry_enabled) {
    // create stor for chemistry data
    int number_of_secondaries(0);
    if (S->secondary_activity_coeff() != Teuchos::null) {
      number_of_secondaries = S->secondary_activity_coeff()->NumVectors();
    }
    chem_data_ = Teuchos::rcp( new chemistry_data (mesh_maps->cell_map(false),
                                                   S->get_total_component_concentration()->NumVectors(),
                                                   S->number_of_minerals(),
                                                   number_of_secondaries,
                                                   S->number_of_ion_exchange_sites(),
                                                   S->number_of_sorption_sites(),
                                                   S->using_sorption(),
                                                   S->use_sorption_isotherms()) );
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
  } else { // no restart, we will call the PKs to allow them to init their auxilary data and massage initial conditions
    Amanzi::timer_manager.start("Flow PK");
    if (flow_enabled) FPK->InitializeAuxiliaryData();
    if (do_picard_) {
      FPK->InitPicard(S->get_time());
      FPK->CommitState(FS);
    }
    Amanzi::timer_manager.stop("Flow PK");
  }



  // write visualization output as requested
  Amanzi::timer_manager.start("I/O");
  if (chemistry_enabled) {

    // get the auxillary data from chemistry
    Teuchos::RCP<Epetra_MultiVector> aux = CPK->get_extra_chemistry_output_data();
    // write visualization data for timestep
    S->write_vis(*visualization, aux, auxnames, chemistry_enabled, true);
  } else {
    // always write the initial visualization dump
    S->write_vis(*visualization, chemistry_enabled, true);
  }

  // write a restart dump if requested (determined in dump_state)
  restart->dump_state(*S);
  Amanzi::timer_manager.stop("I/O");

  Amanzi::timer_manager.start("Flow PK");
  if (flow_enabled) {
    if (ti_mode == STEADY  && flow_model != std::string("Steady State Saturated")) {
      FPK->InitSteadyState(S->get_time(), dTsteady);
    } else if ( ti_mode == TRANSIENT && flow_model !=std::string("Steady State Richards")) {
      FPK->InitTransient(S->get_time(), dTtransient);
    } else if ( (ti_mode == INIT_TO_STEADY || ti_mode == STEADY) && flow_model == std::string("Steady State Saturated")) {
      if (!restart_requested) {
        FPK->InitSteadyState(S->get_time(), dTsteady);
        FPK->InitializeSteadySaturated();
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
  Amanzi::timer_manager.stop("Flow PK");

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
      Amanzi::timer_manager.start("Flow PK");
      if (flow_enabled) {
        if (ti_mode == INIT_TO_STEADY && S->get_last_time() < Tswitch && S->get_time() >= Tswitch) {
          if(out.get() && includesVerbLevel(verbLevel,Teuchos::VERB_LOW,true)) {
            *out << "Steady state computation complete... now running in transient mode." << std::endl;
          }
          // only init the transient problem if we need to
          if (flow_model != "Steady State Richards"  && flow_model != "Steady State Saturated" )  {
            FPK->InitTransient(S->get_time(), dTtransient);
          }
        }
      }

      // find the flow time step
      if (flow_enabled) {
        // only if we are actually running with flow

        if ((ti_mode == STEADY) ||
            (ti_mode == TRANSIENT && flow_model != std::string("Steady State Richards")) ||
            (ti_mode == INIT_TO_STEADY &&
             ( (flow_model == std::string("Steady State Richards") && S->get_time() >= Tswitch) ||
               (flow_model != std::string("Steady State Saturated"))))) {
          flow_dT = FPK->CalculateFlowDt();
        }
      }
      Amanzi::timer_manager.stop("Flow PK");

      if (ti_mode == TRANSIENT || (ti_mode == INIT_TO_STEADY && S->get_time() >= Tswitch)) {
        if (transport_enabled) {
          Amanzi::timer_manager.start("Transport PK");
          double transport_dT_tmp = TPK->EstimateTransportDt();
          if (transport_subcycling == 0) transport_dT = transport_dT_tmp;
          Amanzi::timer_manager.stop("Transport PK");
        }
        if (chemistry_enabled) {
          Amanzi::timer_manager.start("Chemistry PK");
          chemistry_dT = CPK->max_time_step();
          Amanzi::timer_manager.stop("Chemistry PK");
        }
      }

      // take the mpc time step as the min of all suggested time steps
      mpc_dT = std::min(flow_dT, transport_dT);

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
            *out << setprecision(5) << "Resetting the time integrator at time(y) = "
                 << std::fixed << S->get_time()/(365.25*24*60*60) << std::endl;
            mpc_dT = reset_times_dt_.front();
            mpc_dT = TSM.TimeStep(S->get_time(), mpc_dT);
            tslimiter = MPC_LIMITS;
            // now reset the flow time integrator..
            Amanzi::timer_manager.start("Flow PK");
            if (flow_enabled) FPK->InitTransient(S->get_time(), mpc_dT);
            Amanzi::timer_manager.stop("Flow PK");
          }
        }
      }
      // steady flow is special, it might redo a time step, so we print
      // time step info after we've advanced steady flow
      // first advance flow
      Amanzi::timer_manager.start("Flow PK");
      if (flow_enabled) {
        if ((ti_mode == STEADY) ||
            (ti_mode == TRANSIENT && flow_model != std::string("Steady State Richards")) ||
            (ti_mode == INIT_TO_STEADY &&
             ( (flow_model == std::string("Steady State Richards") && S->get_time() >= Tswitch) ||
               (flow_model != std::string("Steady State Saturated"))))) {
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
          FPK->CommitState(FS);
        }
        S->set_final_time(S->initial_time() + mpc_dT);
      }
      Amanzi::timer_manager.stop("Flow PK");
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
        *out << setprecision(5);
        *out << "Cycle = " << iter;
        *out << ",  Time(y) = "<< fixed << S->get_time() / (365.25*60*60*24);
        *out << ",  dT(y) = " << scientific << mpc_dT / (365.25*60*60*24);
        *out << " " << limitstring;
        *out << std::endl;
      }
      // ==============================================================

      // then advance transport and chemistry
      if (ti_mode == TRANSIENT || (ti_mode == INIT_TO_STEADY && S->get_time() >= Tswitch) ) {
        double tc_dT(mpc_dT);
        double c_dT(chemistry_dT);
        int ntc(1);


        if (chemistry_enabled) {
          Amanzi::timer_manager.start("Chemistry PK");
          // reduce chemistry time step according to the
          // ratio with transport time step that is specified
          // in the input file
          double t_dT(transport_dT);
          if (transport_enabled) {
            t_dT = TPK->EstimateTransportDt();
            double ratio(c_dT/t_dT);
            if (ratio > chem_trans_dt_ratio) {
              c_dT = chem_trans_dt_ratio * t_dT;
            }
          }
          if (mpc_dT > c_dT) {
            ntc = floor(mpc_dT/c_dT)+1;
            tc_dT = mpc_dT/static_cast<double>(ntc);
          }

          if(out.get() && includesVerbLevel(verbLevel,Teuchos::VERB_LOW,true)) {
            *out << "Subcycling info: MPC is taking " << ntc << " chemistry subcycling timesteps" << std::endl;
            if (transport_enabled) {
              *out << "  (chemistry sub cycling time step) / (transport time step) = " << tc_dT/t_dT << std::endl;
            }
            *out << "  chemistry subcycling timestep = " << tc_dT << std::endl;
          }
          Amanzi::timer_manager.stop("Chemistry PK");
        }

        // at this time we know the time step that we are going to use during subcycling: tc_dT

        if (chemistry_enabled) {
          // first store the chemistry state
          chem_data_->store(S->free_ion_concentrations(),
                            S->primary_activity_coeff(),
                            S->secondary_activity_coeff(),
                            S->mineral_volume_fractions(),
                            S->mineral_specific_surface_area(),
                            S->total_sorbed(),
                            S->sorption_sites(),
                            S->surface_complex_free_site_conc(),
                            S->ion_exchange_sites(),
                            S->ion_exchange_ref_cation_conc(),
                            S->isotherm_kd(),
                            S->isotherm_freundlich_n(),
                            S->isotherm_langmuir_b());
        }
        // store the total component concentration, so that we
        // can restore it in the case of a chemistry failure
        Epetra_MultiVector tcc_stor(*total_component_concentration_star);

        bool success(true);
        int tries(0);

        do {
          // try to subcycle with tc_dT, if that fails, we will cut that time step and try again
          try {

            // subcycling loop
            for (int iss = 0; iss<ntc; ++iss) {

              // first we do a transport step, or if transport is off, we simply prepare
              // total_component_concentration_star for the chemistry step
              if (transport_enabled) {
                Amanzi::timer_manager.start("Transport PK");
                TPK->Advance(tc_dT);
                if (TPK->get_transport_status() == AmanziTransport::TRANSPORT_STATE_COMPLETE) {
                  // get the transport state and commit it to the state
                  Teuchos::RCP<AmanziTransport::Transport_State> TS_next = TPK->transport_state_next();
                  *total_component_concentration_star = *TS_next->total_component_concentration();
                } else {
                  Errors::Message message("MPC: error... Transport_PK.advance returned an error status");
                  Exceptions::amanzi_throw(message);
                }
                Amanzi::timer_manager.stop("Transport PK");
              } else { // if we're not advancing transport we still need to prepare for chemistry
                *total_component_concentration_star = *S->get_total_component_concentration();
              }

              // second we do a chemistry step, or if chemistry is off, we simply update 
              // total_component_concentration in state
              if (chemistry_enabled) {

                if(out.get() && includesVerbLevel(verbLevel,Teuchos::VERB_LOW,true)) {
                  *out << "Chemistry PK: advancing, current subcycling time step = " << tc_dT << std::endl;
                }
                Amanzi::timer_manager.start("Chemistry PK");
                CPK->advance(tc_dT, total_component_concentration_star);

                Amanzi::timer_manager.stop("Chemistry PK");
                S->update_total_component_concentration(CPK->get_total_component_concentration());
              } else {
                S->update_total_component_concentration(*total_component_concentration_star);
              }

              // all went well, so we can advance intermediate time, and call commit state
              // for each pk
              S->set_intermediate_time(S->intermediate_time() + tc_dT);
              Amanzi::timer_manager.start("Transport PK");
              if (transport_enabled) TPK->CommitState(TS);
              Amanzi::timer_manager.stop("Transport PK");
              Amanzi::timer_manager.start("Chemistry PK");
              if (chemistry_enabled) CPK->commit_state(CS, tc_dT);
              Amanzi::timer_manager.stop("Chemistry PK");
            }
            success = true;
          } catch (const ChemistryException& chem_error) {
            
            // if the chemistry step failed, we back up to the beginning of
            // the chemistry subcycling loop, but to do that we must restore
            // a few things, such as the chemistry state, total component
            // concentration, and back up the intermediate time

            // decrease the chemistry subcycling timestep and adjust the
            // number of subcycles we need to take accordingly
            ntc = 2*ntc;
            tc_dT = 0.5 * tc_dT;

            // increase the retry count
            ++tries;

            // bail if we've cut the subcycling timestep too many times
            if (tries>=3) {
              Errors::Message message("MPC: cut chemistry subcycling time step too many times, bailing...");
              Exceptions::amanzi_throw(message);
            }

            // the the user know that we're backing up due to a chemistry failure
            if(out.get() && includesVerbLevel(verbLevel,Teuchos::VERB_LOW,true)) {
              *out << "Chemistry step failed, reducing chemistry subcycling time step." << std::endl;
              *out << "  new chemistry subcycling time step = " << tc_dT << std::endl;
            }

            // restore chemistry data to the beginning of the subcycling
            chem_data_->retrieve(S->free_ion_concentrations(),
                                 S->primary_activity_coeff(),
                                 S->secondary_activity_coeff(),
                                 S->mineral_volume_fractions(),
                                 S->mineral_specific_surface_area(),
                                 S->total_sorbed(),
                                 S->sorption_sites(),
                                 S->surface_complex_free_site_conc(),
                                 S->ion_exchange_sites(),
                                 S->ion_exchange_ref_cation_conc(),
                                 S->isotherm_kd(),
                                 S->isotherm_freundlich_n(),
                                 S->isotherm_langmuir_b());

            // restore the total component concentration to the beginning of chemistry subcycling
            S->update_total_component_concentration(tcc_stor);
            
            // reset the intermediate time to the beginning
            S->set_intermediate_time(S->initial_time());

            success = false;
          }

        } while (!success);

      }

      // update the times in the state object
      S->advance_time(mpc_dT);

      // ===========================================================
      // we're done with this time step, commit the state
      // in the process kernels

      // if (flow_enabled) FPK->CommitState(FS);
      if (ti_mode == TRANSIENT || (ti_mode == INIT_TO_STEADY && S->get_time() >= Tswitch) ) {
        Amanzi::timer_manager.start("Transport PK");
        if (transport_enabled) TPK->CommitState(TS);
        Amanzi::timer_manager.stop("Transport PK");

        Amanzi::timer_manager.start("Chemistry PK");
        if (chemistry_enabled) CPK->commit_state(CS, mpc_dT);
        Amanzi::timer_manager.stop("Chemistry PK");
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

      Amanzi::timer_manager.start("I/O");
      if (chemistry_enabled) {
        // get the auxillary data
        Teuchos::RCP<Epetra_MultiVector> aux = CPK->get_extra_chemistry_output_data();

        // write visualization data for timestep if requested
        S->write_vis(*visualization, aux, auxnames, chemistry_enabled, force);
      } else {
        S->write_vis(*visualization, chemistry_enabled, force);
      }

      // write restart dump if requested
      restart->dump_state(*S, force);
      Amanzi::timer_manager.stop("I/O");

    }
  }

  // write final visualization dump if no time stepping was done
  if (iter == 0) {
    ++iter;
    S->set_cycle(iter);
    Amanzi::timer_manager.start("I/O");
    S->write_vis(*visualization, false, true);
    Amanzi::timer_manager.stop("I/O");
  }


  // some final output
  if (out.get() && includesVerbLevel(verbLevel,Teuchos::VERB_LOW,true))
  {
    *out << "Cycle = " << iter;
    *out << ",  Time(years) = "<< S->get_time()/ (365.25*60*60*24);
    *out << std::endl;
  }


  // clean up
  delete visualization;
  delete restart;
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
