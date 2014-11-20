#include <utility>
#include <algorithm>
#include <sys/resource.h>

#include "errors.hh"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Epetra_Comm.h"
#include "Epetra_MpiComm.h"

#include "chemistry_pk.hh"
#ifdef ALQUIMIA_ENABLED
#include "alquimia_chemistry_pk.hh"
#endif
#include "MPC.hh"
#include "State.hh"
#include "Darcy_PK.hh"
#include "Richards_PK.hh"
#include "Transport_PK.hh"

#define BOOST_FILESYSTEM_NO_DEPRECATED
#include "boost/filesystem/operations.hpp"
#include "boost/filesystem/path.hpp"
#include "TimeStepManager.hh"

// make sure that we use default parameters
// that are consistent with the input translator
#include "InputParserIS_Defs.hh"

#include "TimerManager.hh"

#include "DataDebug.hh"

namespace Amanzi {

bool reset_info_compfunc(std::pair<double,double> x, std::pair<double,double> y) {
  return (x.first < y.first);
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
  mpc_parameter_list = parameter_list.sublist("MPC");

  vo_ = Teuchos::rcp(new VerboseObject("MPC", mpc_parameter_list));
  Teuchos::OSTab tab = vo_->getOSTab();

  read_parameter_list();

  // let users selectively disable individual process kernels
  // to allow for testing of the process kernels separately
  transport_enabled =
      (mpc_parameter_list.get<std::string>("disable Transport_PK","no") == "no");

  std::string chemistry_model = mpc_parameter_list.get<std::string>("Chemistry Model","Off");
  if (chemistry_model != "Off") {
    chemistry_enabled = true;
  }

  flow_enabled =
      (mpc_parameter_list.get<std::string>("disable Flow_PK","no") == "no");

  if (vo_->os_OK(Teuchos::VERB_LOW)) {
    *vo_->os() << "The following process kernels are enabled: ";
    if (flow_enabled) *vo_->os() << "Flow ";
    if (transport_enabled) *vo_->os() << "Transport ";
    if (chemistry_enabled) *vo_->os() << "Chemistry ";
    *vo_->os() << std::endl;
  }

  if (transport_enabled || flow_enabled || chemistry_enabled) {
    Teuchos::ParameterList state_parameter_list = parameter_list.sublist("State");
    S = Teuchos::rcp(new State(state_parameter_list));
    S->RegisterMesh("domain",mesh_maps);
  }


  // create auxilary state objects for the process models

  // Get the transport solute component names from the MPC parameter lists.
  std::vector<std::string> component_names;
  {
    Teuchos::Array<std::string> comp_names;
    if (mpc_parameter_list.isParameter("component names")) {
      comp_names = mpc_parameter_list.get<Teuchos::Array<std::string> >("component names");
    }
    for (int i = 0; i < comp_names.length(); i++) {
      component_names.push_back(comp_names[i]);
    }
  }

  // chemistry...
  if (chemistry_enabled) {
    Teuchos::ParameterList chemistry_parameter_list = parameter_list.sublist("Chemistry");

#ifdef ALQUIMIA_ENABLED
    if (chemistry_model == "Alquimia") {
      // Start up the chemistry engine.
      if (!chemistry_parameter_list.isParameter("Engine")) 
      {
        Errors::Message msg;
        msg << "MPC::mpc_init(): \n";
        msg << "  No 'Engine' parameter found in parameter list for 'Chemistry'.\n";
        Exceptions::amanzi_throw(msg);
      }
      if (!chemistry_parameter_list.isParameter("Engine Input File")) 
      {
        Errors::Message msg;
        msg << "MPC::mpc_init(): \n";
        msg << "  No 'Engine Input File' parameter found in parameter list for 'Chemistry'.\n";
        Exceptions::amanzi_throw(msg);
      }
      std::string chemEngineName = chemistry_parameter_list.get<std::string>("Engine");
      std::string chemEngineInputFile = chemistry_parameter_list.get<std::string>("Engine Input File");
      chem_engine = Teuchos::rcp(new AmanziChemistry::ChemistryEngine(chemEngineName, chemEngineInputFile));

      // Overwrite the component names with the primary species names from the engine.
      chem_engine->GetPrimarySpeciesNames(component_names);
    }
#endif
    // Initialize the chemistry state.
    CS = Teuchos::rcp(new AmanziChemistry::Chemistry_State(chemistry_parameter_list, component_names, S));

#ifdef ALQUIMIA_ENABLED
    // Set up auxiliary chemistry data using the ChemistryEngine.
    if (chemistry_model == "Alquimia") {
      std::vector<std::string> auxNames;
      chem_engine->GetAuxiliaryOutputNames(auxNames);
      CS->SetAuxDataNames(auxNames); 
    }
#endif
  }

  // flow...
  if (flow_enabled) {
    flow_model = mpc_parameter_list.get<std::string>("Flow model", "Darcy");
    if (flow_model == "Darcy") {
      FPK = Teuchos::rcp(new Flow::Darcy_PK(parameter_list, S));
    } else if (flow_model == "Steady State Saturated") {
      FPK = Teuchos::rcp(new Flow::Darcy_PK(parameter_list, S));
    } else if (flow_model == "Richards") {
      FPK = Teuchos::rcp(new Flow::Richards_PK(parameter_list, S));
    } else if (flow_model == "Steady State Richards") {
      FPK = Teuchos::rcp(new Flow::Richards_PK(parameter_list, S));
    } else {
      std::cout << "MPC: unknown flow model: " << flow_model << std::endl;
      throw std::exception();
    }
  }

  if (flow_model == "Steady State Richards") {
    if (vo_->os_OK(Teuchos::VERB_LOW)) {
      *vo_->os() << "Flow will be off during the transient phase" << std::endl;
    }
  }

  if (transport_enabled) {
#ifdef ALQUIMIA_ENABLED
    if (chemistry_model == "Alquimia") {
      // When Alquimia is used, the Transport PK must interact with the 
      // chemistry engine to obtain boundary values for the components.
      // The component names are fetched from the chemistry engine.
      TPK = Teuchos::rcp(new Transport::Transport_PK(parameter_list, S, CS, chem_engine));
    }
    else {
#endif
      TPK = Teuchos::rcp(new Transport::Transport_PK(parameter_list, S, component_names));
#ifdef ALQUIMIA_ENABLED
    }
#endif
  }

  S->Setup();
  S->InitializeFields();
  S->InitializeEvaluators();
  S->GetMeshPartition("materials");

  if (chemistry_enabled) {
    CS->Initialize();
  }
  if (flow_enabled) {
    FPK->InitializeFields();
  }
  if (transport_enabled) {
    TPK->InitializeFields();
  }

  S->CheckAllFieldsInitialized();

  if (flow_enabled) {
    FPK->Initialize(S.ptr());
  }

  if (transport_enabled) {
    bool subcycling = parameter_list.sublist("MPC").get<bool>("transport subcycling", false);
    transport_subcycling = (subcycling) ? 1 : 0;
    TPK->Initialize(S.ptr());
  }

  if (chemistry_enabled) {
    try {
      if (chemistry_model == "Alquimia") {
#ifdef ALQUIMIA_ENABLED
        CPK = Teuchos::rcp(new AmanziChemistry::Alquimia_Chemistry_PK(parameter_list, CS, chem_engine));
#else
        std::cout << "MPC: Alquimia chemistry model is not enabled for this build.\n";
        throw std::exception();
#endif
      } else if (chemistry_model == "Amanzi") {
        Teuchos::ParameterList chemistry_parameter_list = parameter_list.sublist("Chemistry");
        CPK = Teuchos::rcp(new AmanziChemistry::Chemistry_PK(chemistry_parameter_list, CS));
      } else {
        std::cout << "MPC: unknown chemistry model: " << chemistry_model << std::endl;
        throw std::exception();
      }
      CPK->InitializeChemistry();
    } catch (const Errors::Message& chem_error) {
      std::ostringstream error_message;
      error_message << "MPC:mpc_init(): Chemistry_PK.InitializeChemistry returned an error status:\n";
      error_message << chem_error.what();

      Errors::Message message(error_message.str());
      Exceptions::amanzi_throw(message);
    }   
  } 
  // done creating auxilary state objects and  process models

  // create the observations
  if (parameter_list.isSublist("Observation Data")) {
    Teuchos::ParameterList observation_plist = parameter_list.sublist("Observation Data");
    observations = Teuchos::rcp(new Amanzi::Unstructured_observations(observation_plist, output_observations, comm));

    if (mpc_parameter_list.isParameter("component names")) {
      Teuchos::Array<std::string> comp_names;
      comp_names = mpc_parameter_list.get<Teuchos::Array<std::string> >("component names");
      observations->RegisterComponentNames(comp_names.toVector());
    }
  } else {
    observations = Teuchos::null;
  }

  // create the visualization object
  if (parameter_list.isSublist("Visualization Data"))  {
    Teuchos::ParameterList vis_parameter_list = parameter_list.sublist("Visualization Data");
    visualization = Teuchos::ptr(new Amanzi::Visualization(vis_parameter_list, comm));
    visualization->set_mesh(mesh_maps);
    visualization->CreateFiles();
  } else {  // create a dummy vis object
    visualization = Teuchos::ptr(new Amanzi::Visualization());
  }


  // create the restart object
  if (parameter_list.isSublist("Checkpoint Data")) {
    Teuchos::ParameterList checkpoint_parameter_list = parameter_list.sublist("Checkpoint Data");
    restart = Teuchos::ptr(new Amanzi::Checkpoint(checkpoint_parameter_list, comm));
  } else {
    restart = Teuchos::ptr(new Amanzi::Checkpoint());
  }

  // are we restarting from a file?
  // first assume we're not
  restart_requested = false;


  if (flow_enabled) {
    if (parameter_list.isSublist("Walkabout Data")) {
      Teuchos::ParameterList walkabout_parameter_list = parameter_list.sublist("Walkabout Data");
      walkabout = Teuchos::ptr(new Amanzi::Checkpoint(walkabout_parameter_list, comm));
    } else {
      walkabout = Teuchos::ptr(new Amanzi::Checkpoint());
    }    
  }
  
  // then check if indeed we are
  if (mpc_parameter_list.isSublist("Restart from Checkpoint Data File")) {

    restart_requested = ! mpc_parameter_list.sublist("Restart from Checkpoint Data File").get<bool>("initialize from checkpoint data file and do not restart",false);

    Teuchos::ParameterList& restart_parameter_list =
        mpc_parameter_list.sublist("Restart from Checkpoint Data File");

    restart_from_filename = restart_parameter_list.get<std::string>("Checkpoint Data File Name");

    if (restart_requested) {
      if (vo_->os_OK(Teuchos::VERB_LOW)) {
        *vo_->os() << "Restarting from checkpoint file: " << restart_from_filename << std::endl;
      }
    } else {
      if (vo_->os_OK(Teuchos::VERB_LOW)) {
        *vo_->os() << "Initializing data from checkpoint file: " << restart_from_filename << std::endl
                   << "    (Ignoring all initial conditions.)" << std::endl;
      }
    }

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
    dTtransient = 1e+99;

    do_picard_ = steady_list.get<bool>("Use Picard",false);
  } else if ( ti_list.isSublist("Transient") ) {
    ti_mode = TRANSIENT;

    Teuchos::ParameterList& transient_list = ti_list.sublist("Transient");

    T0 = transient_list.get<double>("Start");
    T1 = transient_list.get<double>("End");
    dTtransient =  transient_list.get<double>("Initial Time Step");
    dTsteady = 1e+99;

    do_picard_ = false;
  } else if ( ti_list.isSublist("Transient with Static Flow") ) {
    ti_mode = TRANSIENT_STATIC_FLOW;

    Teuchos::ParameterList& transient_list = ti_list.sublist("Transient with Static Flow");

    T0 = transient_list.get<double>("Start");
    T1 = transient_list.get<double>("End");
    dTtransient =  transient_list.get<double>("Initial Time Step");
    dTsteady = 1e+99;

    do_picard_ = false;
  } else {
    Errors::Message message("MPC: no valid Time Integration Mode was specified, you must specify exactly one of Initialize To Steady, Steady, or Transient.");
    Exceptions::amanzi_throw(message);
  }

  if (mpc_parameter_list.isSublist("Time Period Control")) {
    Teuchos::ParameterList& tpc_list =  mpc_parameter_list.sublist("Time Period Control");
    Teuchos::Array<double> reset_times    = tpc_list.get<Teuchos::Array<double> >("Start Times");
    Teuchos::Array<double> reset_times_dt = tpc_list.get<Teuchos::Array<double> >("Initial Time Step");
    if (reset_times.size() != reset_times_dt.size()) {
      Errors::Message message("You must specify the same number of Reset Times and Initial Time Steps under Time Period Control");
      Exceptions::amanzi_throw(message);
    }
    Teuchos::Array<double>::const_iterator it_tim;
    Teuchos::Array<double>::const_iterator it_dt;
    for (it_tim = reset_times.begin(), it_dt = reset_times_dt.begin();
         it_tim != reset_times.end();
         ++it_tim, ++it_dt) {
      reset_info_.push_back(std::make_pair(*it_tim, *it_dt));
    }
    // now we sort in ascending order by time
    std::sort(reset_info_.begin(), reset_info_.end(), reset_info_compfunc);
  }

  ti_rescue_factor_ = mpc_parameter_list.get<double>("time integration rescue reduction factor",0.5);
}


/* *******************************************************************/
void MPC::cycle_driver() {
  Teuchos::OSTab tab = vo_->getOSTab();

  // Amanzi::timer_manager.add("AnalyticJacobian", Amanzi::Timer::ACCUMULATE);
  // Amanzi::timer_manager.add("Function", Amanzi::Timer::ACCUMULATE);
  // Amanzi::timer_manager.add("Update precon", Amanzi::Timer::ACCUMULATE);
  // Amanzi::timer_manager.add("Apply precon", Amanzi::Timer::ACCUMULATE);

  // start timers
  Amanzi::timer_manager.add("Chemistry PK", Amanzi::Timer::ACCUMULATE);
  Amanzi::timer_manager.add("Flow PK", Amanzi::Timer::ACCUMULATE);
  Amanzi::timer_manager.add("Transport PK", Amanzi::Timer::ACCUMULATE);
  Amanzi::timer_manager.add("I/O", Amanzi::Timer::ACCUMULATE);

  // create the time step manager
  Teuchos::Ptr<Amanzi::TimeStepManager> TSM = Teuchos::ptr(new TimeStepManager());
  // register visualization times with the time step manager
  visualization->RegisterWithTimeStepManager(TSM);
  restart->RegisterWithTimeStepManager(TSM);
  // register observation times with the time step manager
  if (observations != Teuchos::null) observations->RegisterWithTimeStepManager(TSM);
  // register reset_times
  for(std::vector<std::pair<double,double> >::const_iterator it = reset_info_.begin();
      it != reset_info_.end(); ++it) TSM->RegisterTimeEvent(it->first);
  // if this is an init to steady run, register the switchover time
  if (ti_mode == INIT_TO_STEADY) TSM->RegisterTimeEvent(Tswitch);
  // register the final time
  TSM->RegisterTimeEvent(T1);

  enum time_step_limiter_type {FLOW_LIMITS, TRANSPORT_LIMITS, CHEMISTRY_LIMITS, MPC_LIMITS};
  time_step_limiter_type tslimiter;

  if (transport_enabled || flow_enabled || chemistry_enabled) {
    S->set_time(T0);  // start at time T=T0;
    S->set_initial_time(T0);
    S->set_intermediate_time(Tswitch);
  }

  int iter = 0;  // set the iteration counter to zero
  S->set_cycle(iter);

  double restart_dT(1.0e99);
  // read the checkpoint file as requested
  if (restart_requested == true) {
    // re-initialize the state object
    restart_dT = ReadCheckpoint(comm, Teuchos::ptr(&*S), restart_from_filename);
    iter = S->cycle();
    for (std::vector<std::pair<double,double> >::iterator it = reset_info_.begin();
         it != reset_info_.end(); ++it) {
      if (it->first < S->time()) it = reset_info_.erase(it);
      if (it == reset_info_.end() ) break;
    }
    S->set_initial_time(S->time());
  } else { // no restart, we will call the PKs to allow them to init their auxilary data and massage initial conditions
    Amanzi::timer_manager.start("Flow PK");
    if (flow_enabled) FPK->InitializeAuxiliaryData();
    if (do_picard_) {
      FPK->InitPicard(S->time());
      FPK->CommitState(0.0, S.ptr());
    }
    Amanzi::timer_manager.stop("Flow PK");
  }

  Amanzi::timer_manager.start("Flow PK");
  if (flow_enabled) {
    if (ti_mode == STEADY  && flow_model != std::string("Steady State Saturated")) {
      // this is the case Richards time stepped to steady state
      // we simply initialize the Flow PK accordingly and the
      // time stepping loop below takes care of the rest
      FPK->InitSteadyState(S->time(), dTsteady);
    } else if ( ti_mode == TRANSIENT && flow_model !=std::string("Steady State Richards")) {
      FPK->InitTransient(S->time(), dTtransient);
    } else if ( ti_mode == TRANSIENT_STATIC_FLOW ) {
      FPK->InitTransient(S->time(), dTtransient);
    } else if ( (ti_mode == INIT_TO_STEADY || ti_mode == STEADY) && flow_model == std::string("Steady State Saturated")) {
      // this is the case where we need to solve the Darcy problem first
      // and then either stop (in the STEADY case), or move to the switch time
      // to get the transient part going (INIT_TO_STEADY).
      // note that if a restart was requested, we get the flow field from
      // the checkpoint file we skip the linear Darcy solve and proceed with
      // the initialization of the transient problem
      if (!restart_requested) {
        FPK->InitSteadyState(S->time(), dTsteady);
        FPK->InitializeSteadySaturated();
        FPK->CommitState(0.0, S.ptr());
        if (ti_mode == INIT_TO_STEADY) S->advance_time(Tswitch-T0);
        if (ti_mode == STEADY)         S->advance_time(T1-T0);
      } else {
        FPK->InitTransient(S->time(), dTtransient);
      }
    } else if (ti_mode == INIT_TO_STEADY) {
      if (!restart_requested) {
        if (S->time() < Tswitch) {
          FPK->InitSteadyState(S->time(), dTsteady);
        } else {
          if (flow_model !=std::string("Steady State Richards")) {
            FPK->InitTransient(S->time(), dTtransient);
          }
        }
      } else {
        if (S->time() < Tswitch) {
          FPK->InitSteadyState(S->time(), restart_dT);
        } else {
          if (flow_model != std::string("Steady State Richards")) {
            // only initialize here if we're not on a time integrator restart time
            // in that case we initialize later
            std::vector<std::pair<double,double> >::const_iterator it;
            for (it = reset_info_.begin(); it != reset_info_.end(); ++it) {
              if (it->first == S->time()) break;
            }
            if (it == reset_info_.end() && fabs(S->time()-Tswitch)>1e-7) {
              FPK->InitTransient(S->time(), restart_dT);
            }
          }
        }
      }
    } 
  }
  Amanzi::timer_manager.stop("Flow PK");


  // write visualization output as requested
  Amanzi::timer_manager.start("I/O");
  //visualization->set_mesh(mesh_maps);
  visualization->CreateFiles();


  if (!restart_requested) {
    if (flow_enabled) FPK->UpdateAuxilliaryData();
  }

  if (!restart_requested) {
    if (vo_->os_OK(Teuchos::VERB_MEDIUM)) {
      *vo_->os() << "Cycle " << S->cycle() << ": writing visualization file" << std::endl;
    }
    if (chemistry_enabled) {
      // get the auxillary data from chemistry
      Teuchos::RCP<Epetra_MultiVector> aux = CPK->get_extra_chemistry_output_data();
      // write visualization data for timestep
      WriteVis(visualization,S.ptr()); // TODO(berndt) make sure that aux names are used for vis
    } else {
      //always write the initial visualization dump
      WriteVis(visualization,S.ptr());
    }

    // write a restart dump if requested (determined in dump_state)
    if (restart->DumpRequested(S->cycle(),S->time())) {
      if (vo_->os_OK(Teuchos::VERB_MEDIUM)) {
        *vo_->os() << "Cycle " << S->cycle() << ": writing checkpoint" << std::endl;
      }
      if (ti_mode == STEADY || ti_mode == INIT_TO_STEADY) {
        WriteCheckpoint(restart,S.ptr(),dTsteady);
      } else if ( (ti_mode == TRANSIENT) || (ti_mode == TRANSIENT_STATIC_FLOW) ) {
        WriteCheckpoint(restart,S.ptr(),dTtransient);
      }
    }

    // write walkabout data if requested
    if (flow_enabled) {
      if (walkabout->DumpRequested(S->cycle(), S->time())) {
        if (vo_->os_OK(Teuchos::VERB_MEDIUM)) {
          *vo_->os() << "Cycle " << S->cycle() << ": writing walkabout data" << std::endl;
        }
        FPK->WriteWalkabout(walkabout.ptr());
      }
    }
  }

  Amanzi::timer_manager.stop("I/O");

  if (flow_enabled || transport_enabled || chemistry_enabled) {
    if (observations != Teuchos::null) {
      if (observations->DumpRequested(S->cycle(), S->time())) {
        if (flow_enabled) FPK->UpdateAuxilliaryData();
        observations->MakeObservations(*S);
      }
    }
    // we need to create an EpetraMultiVector that will store the
    // intermediate value for the total component concentration
    if (chemistry_enabled || transport_enabled) {
      total_component_concentration_star =
          Teuchos::rcp(new Epetra_MultiVector(*S->GetFieldData("total_component_concentration")->ViewComponent("cell", true)));
    }
    // then start time stepping
    while ((S->time() < T1) && ((end_cycle == -1) || (iter <= end_cycle))) {

      // log that we are starting a time step
      if (vo_->os_OK(Teuchos::VERB_LOW)) {
        *vo_->os() << std::setprecision(5) << std::endl
                   << "Cycle " << iter
                   << ": starting time step at time(y) = "<< std::scientific << S->time() / (365.25*60*60*24)
                   << std::endl;
      }

      // determine the time step we are now going to take
      double chemistry_dT = 1e+99, transport_dT = 1e+99, flow_dT = 1e+99;
      double mpc_dT = 1e+99, limiter_dT = 1e+99, observation_dT = 1e+99;

      // Update our reset times (delete the next one if we just did it)
      if (!restart_requested) {
        if (!reset_info_.empty()) {
          if (S->last_time() >= reset_info_.front().first) {
            reset_info_.erase(reset_info_.begin());
          }
        }
      }

      // catch the switchover time to transient
      Amanzi::timer_manager.start("Flow PK");
      if (flow_enabled) {
        //if (ti_mode == INIT_TO_STEADY && S->last_time() < Tswitch && S->time() >= Tswitch) {
        if (ti_mode == INIT_TO_STEADY && fabs(S->time() - Tswitch) < 1e-7) {
          if (vo_->os_OK(Teuchos::VERB_LOW)) {
            *vo_->os() << "Steady state computation complete... now running in transient mode." << std::endl
                       << "Tswitch = " << Tswitch 
                       << " S->time() = " << S->time() << " S->last_time() = " << S->last_time() << std::endl;
          }
          // only init the transient problem if we need to
          if (flow_model != "Steady State Richards") { //  && flow_model != "Steady State Saturated" )  {
            FPK->InitTransient(S->time(), dTtransient);
          }
        }
      }

      // find the flow time step
      if (flow_enabled && (ti_mode != TRANSIENT_STATIC_FLOW)) {
        // only if we are actually running with flow

        if (!restart_requested) {
          if ((ti_mode == STEADY) ||
              (ti_mode == TRANSIENT && flow_model != std::string("Steady State Richards")) ||
              (ti_mode == INIT_TO_STEADY &&
               ( (flow_model == std::string("Steady State Richards") && S->time() >= Tswitch) ||
                 (flow_model == std::string("Steady State Saturated") && S->time() >= Tswitch) ||
                 (flow_model == std::string("Richards") ) ) ) ) {
            flow_dT = FPK->get_dt();
          }
        } else {
          flow_dT = restart_dT;
        }
      }
      Amanzi::timer_manager.stop("Flow PK");

      if ( (ti_mode == TRANSIENT) || 
	   (ti_mode == INIT_TO_STEADY && S->time() >= Tswitch) || 
	   (ti_mode == TRANSIENT_STATIC_FLOW) ) {
        if (transport_enabled) {
          Amanzi::timer_manager.start("Transport PK");
          double transport_dT_tmp = TPK->get_dt();
          if (transport_subcycling == 0) transport_dT = transport_dT_tmp;
          Amanzi::timer_manager.stop("Transport PK");
        }
        if (chemistry_enabled) {
          Amanzi::timer_manager.start("Chemistry PK");
          chemistry_dT = CPK->time_step();
          Amanzi::timer_manager.stop("Chemistry PK");
        }
      }

      // take the mpc time step as the min of all suggested time steps
      mpc_dT = std::min(flow_dT, transport_dT);

      // take the mpc time step as the min of the last limiter and itself
      mpc_dT = TSM->TimeStep(S->time(), mpc_dT);

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
      if (ti_mode == INIT_TO_STEADY && S->time() >= Tswitch && S->last_time() < Tswitch) {
        mpc_dT = std::min( mpc_dT, dTtransient );
        tslimiter = MPC_LIMITS;
      }

      // make sure that if we are currently on a reset time, to reset the time step
      if (! (ti_mode == STEADY)) {
        if (!reset_info_.empty()) {
          // this is probably iffy...
          if (S->time() == reset_info_.front().first) {
            if (vo_->os_OK(Teuchos::VERB_LOW)) {
              *vo_->os() << std::setprecision(5) << "Resetting the time integrator at time(y) = "
                         << std::fixed << S->time()/(365.25*24*60*60) << std::endl;
            }
            mpc_dT = reset_info_.front().second;
            mpc_dT = TSM->TimeStep(S->time(), mpc_dT);
            tslimiter = MPC_LIMITS;
            // now reset the flow time integrator..
            if (S->time() >= Tswitch) {
              // we are now in transient mode
              Amanzi::timer_manager.start("Flow PK");
              if (flow_enabled) FPK->InitTransient(S->time(), mpc_dT);
              Amanzi::timer_manager.stop("Flow PK");
            } else {
              // we are still in steady state mode
              // and cannot handle this case, currently
              // so we throw an exception
              Errors::Message message("MPC: Restarting the time integrator before switch time is not supported.");
              Exceptions::amanzi_throw(message);	      
            }
          }
        }
      }

      if (restart_requested) {
        //mpc_dT = restart_dT;
        restart_requested = false;
      }

      if (vo_->os_OK(Teuchos::VERB_MEDIUM)) {
        *vo_->os() << std::setprecision(5)
                   << "Cycle " << iter
                   << ": proposed time step before flow step dT(y) = " 
                   << std::scientific << mpc_dT / (365.25*60*60*24) << std::endl;
      }


      // steady flow is special, it might redo a time step, so we print
      // time step info after we've advanced steady flow
      // first advance flow
      Amanzi::timer_manager.start("Flow PK");
      if (flow_enabled) {
        if ((ti_mode == STEADY) ||
            (ti_mode == TRANSIENT && (flow_model != std::string("Steady State Richards"))) || 
            (ti_mode == INIT_TO_STEADY &&
             ( (flow_model == std::string("Steady State Richards") && S->time() >= Tswitch) ||
               (flow_model == std::string("Steady State Saturated") && S->time() >= Tswitch) ||
               (flow_model == std::string("Richards")))) ) {
          double actual_dT(mpc_dT);
          FPK->Advance(mpc_dT, actual_dT);
          // adjust the mpc timestep to the one that the flow
          // pk actually took
          if (actual_dT != mpc_dT) {
            mpc_dT = actual_dT;
            tslimiter = FLOW_LIMITS;
          }
          FPK->CommitState(mpc_dT, S.ptr());
        }
      }
      S->set_final_time(S->initial_time() + mpc_dT);
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

      // then advance transport and chemistry
      if ( (ti_mode == TRANSIENT) || 
           (ti_mode == INIT_TO_STEADY && S->time() >= Tswitch) || 
           (ti_mode == TRANSIENT_STATIC_FLOW) ) {
        double tc_dT(mpc_dT);
        double c_dT(chemistry_dT);

        S->set_intermediate_time(S->initial_time());

        double t_dT(transport_dT);
        if (chemistry_enabled) {
          Amanzi::timer_manager.start("Chemistry PK");

          // The chemistry and transport PKs will both proceed with the 
          // lesser of their respective timesteps.
          if (transport_enabled) {
            t_dT = TPK->get_dt();
          }
          tc_dT = std::min(std::min(t_dT, c_dT), mpc_dT);

          if ((vo_->os_OK(Teuchos::VERB_LOW)) && (tc_dT < mpc_dT)) {
            *vo_->os() << "Subcycling info: MPC is allowing chemistry to subcycle" << std::endl;
            *vo_->os() << "  chemistry subcycling timestep = " << tc_dT << std::endl;
          }
          Amanzi::timer_manager.stop("Chemistry PK");
        }

        // Now advance one step at a time till we reach the end of the flow timestep.
        while (S->intermediate_time() < S->final_time()) {
          // Try to subcycle with tc_dT. If a timestep fails, we rely on the Chemistry PK 
          // to select a new one and try again.
          bool chem_step_succeeded = false;
          try {
            // First we do a transport step. If transport is off, we simply prepare
            // total_component_concentration_star for the chemistry step
            if (transport_enabled) {
              Amanzi::timer_manager.start("Transport PK");
              double tmp_dT;
              int ok = TPK->Advance(tc_dT, tmp_dT);
              if (ok == 0) {
                *total_component_concentration_star = *TPK->total_component_concentration()->ViewComponent("cell", true);
              } else {
                Errors::Message message("MPC: error... Transport_PK.advance returned an error status");
                Exceptions::amanzi_throw(message);
              }
            } else if (chemistry_enabled) { // if we're not advancing transport we still need to prepare for chemistry
              *total_component_concentration_star = *S->GetFieldData("total_component_concentration")->ViewComponent("cell", true);
            }

            // Second, we do a chemistry step. If chemistry is off, we simply update
            // total_component_concentration in state
            if (chemistry_enabled) {
              Amanzi::timer_manager.start("Chemistry PK");
              CPK->Advance(tc_dT, total_component_concentration_star);
              Amanzi::timer_manager.stop("Chemistry PK");

              chem_step_succeeded = true;
              *S->GetFieldData("total_component_concentration","state")->ViewComponent("cell", true)
                = * CPK->get_total_component_concentration();
            } else if (transport_enabled) {
              *S->GetFieldData("total_component_concentration","state")->ViewComponent("cell", true)
                  = *total_component_concentration_star;
            }

            // all went well, so we advance the intermediate time and call 
            // CommitState() on the chemistry PK.
            S->set_intermediate_time(S->intermediate_time() + tc_dT);
            Amanzi::timer_manager.start("Chemistry PK");
            if (chemistry_enabled) CPK->CommitState(CS, S->intermediate_time());
            Amanzi::timer_manager.stop("Chemistry PK");
          }
          catch (const Errors::Message& chem_error) {
            // If the chemistry step failed, we have to try again.
            // Let the user know that the chemistry step failed.
            if (vo_->os_OK(Teuchos::VERB_LOW)) {
              *vo_->os() << chem_error.what() << std::endl
                << "Chemistry step failed, reducing chemistry time step." << std::endl;
            } // end if
          } // end catch

          if (chemistry_enabled) {
            // Compute the next chemistry time step.
            c_dT = CPK->time_step();

            // If the new time step is essentially equal to the old one and the last one 
            // failed, chop it in half.
	    if (!chem_step_succeeded && (c_dT/tc_dT > 0.99)) {
              c_dT *= 0.5;
    	    } 
          } // end if

          // Factor in the time remaining till the end of the MPC time step.
          double remaining_dT = S->final_time() - S->intermediate_time();
          tc_dT = std::min(std::min(t_dT, c_dT), remaining_dT);
        } // end while
      } // end if

      // update the times in the state object
      S->advance_time(mpc_dT);
      S->set_initial_time(S->time());

      // We're done with this time step, so commit the state in the process kernels

      if ( (ti_mode == TRANSIENT) || 
           (ti_mode == INIT_TO_STEADY && S->time() >= Tswitch) ||
      	   (ti_mode == TRANSIENT_STATIC_FLOW) ) {
        Amanzi::timer_manager.start("Transport PK");
        if (transport_enabled && !chemistry_enabled) TPK->CommitState(0.0, S.ptr());
        Amanzi::timer_manager.stop("Transport PK");

        Amanzi::timer_manager.start("Chemistry PK");
        if (chemistry_enabled) CPK->CommitState(CS, S->time());
        Amanzi::timer_manager.stop("Chemistry PK");
      }

      if (vo_->os_OK(Teuchos::VERB_LOW)) {
        *vo_->os() << std::setprecision(5)
                   << "Cycle " << iter
                   << ": complete, new time = " << S->time() / (365.25*60*60*24)
                   << " " << limitstring << std::endl;
      }

      // advance the iteration count
      iter++;
      S->set_cycle(iter);


      // make observations
      if (observations != Teuchos::null) {
        if (observations->DumpRequested(S->cycle(), S->time())) {
          if (flow_enabled) FPK->UpdateAuxilliaryData();
          observations->MakeObservations(*S);
        }
      }

      // write visualization if requested
      // force a vis dump and checkpoint in certain cases,
      // such as at the end of the simulation, and
      // at switch-over time
      bool force(false);
      if (abs(S->time() - T1) < 1e-7) {
        force = true;
      }

      if (ti_mode == INIT_TO_STEADY)
        if (abs(S->time() - Tswitch) < 1e-7) {
          force = true;
        }

      // figure out if in the next iteration, we
      // will reset the time integrator, if so we
      // force a checkpoint
      bool force_checkpoint(false);
      if (ti_mode != STEADY)
        if (!reset_info_.empty())
          if (S->time() == reset_info_.front().first)
            force_checkpoint = true;

      if (force || force_checkpoint || visualization->DumpRequested(S->cycle(), S->time()) ||
          restart->DumpRequested(S->cycle(), S->time())) {
        if (flow_enabled) FPK->UpdateAuxilliaryData();
      }

      Amanzi::timer_manager.start("I/O");
      if (chemistry_enabled) {
        // get the auxillary data
        Teuchos::RCP<Epetra_MultiVector> aux = CPK->get_extra_chemistry_output_data();
        if (force || visualization->DumpRequested(S->cycle(), S->time())) {
	  if (!visualization->is_disabled()) {
            if (vo_->os_OK(Teuchos::VERB_MEDIUM)) {
              *vo_->os() << "Cycle " << S->cycle() << ": writing visualization file" << std::endl;
	    }
	    WriteVis(visualization,S.ptr());
	  }
        }
      } else {
        if (force || visualization->DumpRequested(S->cycle(), S->time())) {
	  if (!visualization->is_disabled()) {
            if (vo_->os_OK(Teuchos::VERB_MEDIUM)) {
              *vo_->os() << "Cycle " << S->cycle() << ": writing visualization file" << std::endl;
	    }
	    WriteVis(visualization,S.ptr());
	  }
        }
      }

      // write restart dump if requested
      if (force || force_checkpoint || restart->DumpRequested(S->cycle(), S->time())) {
	if (!restart->is_disabled()) {
          if (vo_->os_OK(Teuchos::VERB_MEDIUM)) {
            *vo_->os() << "Cycle " << S->cycle() << ": writing checkpoint file" << std::endl;
	  }
	  WriteCheckpoint(restart,S.ptr(),mpc_dT);
	}
      }

      // write walkabout data if requested
      if (flow_enabled) {
        if (force || force_checkpoint || walkabout->DumpRequested(S->cycle(), S->time())) {
	  if (!walkabout->is_disabled()) {
            if (vo_->os_OK(Teuchos::VERB_MEDIUM)) {
              *vo_->os() << "Cycle " << S->cycle() << ": writing walkabout data" << std::endl;
	    }
	    FPK->WriteWalkabout(walkabout.ptr());
	  }
        }
      }
      
      Amanzi::timer_manager.stop("I/O");
    }
  }
  if (flow_enabled) FPK->UpdateAuxilliaryData();

  // write final visualization dump and checkpoint
  // if no time stepping was done
  if (iter == 0) {
    ++iter;
    S->set_cycle(iter);
    Amanzi::timer_manager.start("I/O");
    if (!visualization->is_disabled()) {
      if (vo_->os_OK(Teuchos::VERB_MEDIUM)) {
        *vo_->os() << "Cycle " << S->cycle() << ": writing visualization file" << std::endl;
      }
      if (flow_enabled) FPK->UpdateAuxilliaryData();
      WriteVis(visualization,S.ptr());
    }
      
    if (!restart->is_disabled()) {
      if (vo_->os_OK(Teuchos::VERB_MEDIUM)) {
        *vo_->os() << "Cycle " << S->cycle() << ": writing checkpoint file" << std::endl;
      }
      WriteCheckpoint(restart,S.ptr(),0.0);
    }

    if (flow_enabled) {
      if (!walkabout->is_disabled()) {
        if (vo_->os_OK(Teuchos::VERB_MEDIUM)) {
          *vo_->os() << "Cycle " << S->cycle() << ": writing walkabout file" << std::endl;
	}
	FPK->WriteWalkabout(walkabout.ptr());
      }
    }

    Amanzi::timer_manager.stop("I/O");
  }


  // some final output
  if (vo_->os_OK(Teuchos::VERB_LOW)) {
    *vo_->os() << "Simulation complete at cycle " << iter
               << " and Time(y) = "<< S->time() / (365.25*60*60*24) << std::endl;
  }
  
  if (vo_->os_OK(Teuchos::VERB_MEDIUM)) {
    Epetra_Map cell_map = mesh_maps->cell_map(false);
    double mem = rss_usage();
    
    double percell(mem);
    if (cell_map.NumMyElements() > 0) {
      percell = mem/cell_map.NumMyElements();
    }

    double max_percell(0.0);
    double min_percell(0.0);
    comm->MinAll(&percell,&min_percell,1);
    comm->MaxAll(&percell,&max_percell,1);

    double total_mem(0.0);
    double max_mem(0.0);
    double min_mem(0.0);
    comm->SumAll(&mem,&total_mem,1);
    comm->MinAll(&mem,&min_mem,1);
    comm->MaxAll(&mem,&max_mem,1);
    
    *vo_->os() << std::endl
               << "Memory usage (high water mark):" << std::endl
               << std::fixed << std::setprecision(1);
    *vo_->os() << "  Maximum per core:   " << std::setw(7) << max_mem 
               << " MBytes,  maximum per cell: " << std::setw(7) << max_percell*1024*1024 
               << " Bytes" << std::endl;
    *vo_->os() << "  Minumum per core:   " << std::setw(7) << min_mem 
               << " MBytes,  minimum per cell: " << std::setw(7) << min_percell*1024*1024 
               << " Bytes" << std::endl;
    *vo_->os() << "  Total:              " << std::setw(7) << total_mem 
               << " MBytes,  total per cell:   " << std::setw(7) << total_mem/cell_map.NumGlobalElements()*1024*1024 
               << " Bytes" << std::endl;
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
