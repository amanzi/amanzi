#include <sstream>
#include <string>
#include <algorithm>

#include "errors.hh"
#include "exceptions.hh"
#include "dbc.hh"

#include "InputParserIS.hh"
#include "InputParserIS_Defs.hh"

namespace Amanzi {
namespace AmanziInput {

/* ******************************************************************
* Empty
****************************************************************** */
Teuchos::ParameterList InputParserIS::CreateTimePeriodControlList_(Teuchos::ParameterList* plist)
{
  Teuchos::ParameterList tpc_list;

  Teuchos::Array<double> start_times;
  Teuchos::Array<double> initial_time_step;

  std::map<double, double> time_map;
  double default_initial_time_step(RESTART_TIME_STEP);
  
  /* EIB: proposed v1.2.2 update - Add Maximum Cycle Number for Transient modes */
  double max_cycle_number;

  if (plist->isSublist("Execution Control")) {
    Teuchos::ParameterList exe_sublist = plist->sublist("Execution Control");

    // get the default initial time step
    if (exe_sublist.isSublist("Time Period Control")) {
      default_initial_time_step = exe_sublist.sublist("Time Period Control").get<double>("Default Initial Time Step",RESTART_TIME_STEP);
    }

    if (exe_sublist.isParameter("Flow Model")) {
      std::string flow_model = exe_sublist.get<std::string>("Flow Model");
      if (flow_model != "Off") { // we need to process boudary conditions to possibly add additional time periods

        Teuchos::ParameterList& bc_sublist = plist->sublist("Boundary Conditions");

        for (Teuchos::ParameterList::ConstIterator i = bc_sublist.begin(); i != bc_sublist.end(); i++) {
          // look at sublists
          if (bc_sublist.isSublist(bc_sublist.name(i))) {
            Teuchos::ParameterList& bc = bc_sublist.sublist(bc_sublist.name(i));

            Teuchos::ParameterList bc_list;

            if (bc.isSublist("BC: Flux")) {
              bc_list = bc.sublist("BC: Flux");
            } else if (bc.isSublist("BC: Uniform Pressure")) {
              bc_list = bc.sublist("BC: Uniform Pressure");
            } else if (bc.isSublist("BC: Seepage")) {
              bc_list = bc.sublist("BC: Seepage");
            } else if (bc.isSublist("BC: Zero Flow")) {
              bc_list = bc.sublist("BC: Zero Flow");
            }

            if (bc_list.isParameter("Times")) {
              Teuchos::Array<double> times = bc_list.get<Teuchos::Array<double> >("Times");

              for (Teuchos::Array<double>::const_iterator times_it = times.begin();
                   times_it != times.end(); ++times_it) {
                // skip the first one, there is no jump
                if (times_it != times.begin()) {
                  time_map[*times_it] = default_initial_time_step;
                }
              }
            }
          }
        }

        Teuchos::ParameterList& src_sublist = plist->sublist("Sources");

        for (Teuchos::ParameterList::ConstIterator i = src_sublist.begin(); i != src_sublist.end(); i++) {
          // look at sublists
          if (src_sublist.isSublist(src_sublist.name(i))) {
            Teuchos::ParameterList& src = src_sublist.sublist(src_sublist.name(i));

            Teuchos::ParameterList src_list;

            if (src.isSublist("Source: Uniform")) {
              src_list = src.sublist("Source: Uniform");
            } else if (src.isSublist("Source: Volume Weighted")) {
              src_list = src.sublist("Source: Volume Weighted");
            } else if (src.isSublist("Source: Permeability Weighted")) {
              src_list = src.sublist("Source: Permeability Weighted");
            }

            if (src_list.isParameter("Times")) {
              Teuchos::Array<double> times = src_list.get<Teuchos::Array<double> >("Times");

              for (Teuchos::Array<double>::const_iterator times_it = times.begin();
                   times_it != times.end(); ++times_it) {
                // skip the first one, there is no jump
                if (times_it != times.begin()) {
                  time_map[*times_it] = default_initial_time_step;
                }
              }
            }
          }
        }
      }
    }

    // add the these last so that the default initial time steps get overwritten
    if (exe_sublist.isSublist("Time Period Control")) {
      start_times = exe_sublist.sublist("Time Period Control").get<Teuchos::Array<double> >("Start Times");
      initial_time_step = exe_sublist.sublist("Time Period Control").get<Teuchos::Array<double> >("Initial Time Step");

      Teuchos::Array<double>::const_iterator initial_time_step_it = initial_time_step.begin();
      for (Teuchos::Array<double>::const_iterator start_times_it = start_times.begin();
           start_times_it != start_times.end(); ++start_times_it) {
        time_map[*start_times_it] = *initial_time_step_it;
        ++initial_time_step_it;
      }
    }
  for (std::map<double,double>::const_iterator map_it = time_map.begin();
       map_it != time_map.end(); ++map_it) {
  }

    // delete the start, switch, and end times, since the user must specify initial time steps for those seperately
    if (exe_sublist.isSublist("Time Integration Mode")) {
      if (exe_sublist.sublist("Time Integration Mode").isSublist("Initialize To Steady")) {
        double start_time = exe_sublist.sublist("Time Integration Mode").sublist("Initialize To Steady").get<double>("Start");
        double switch_time = exe_sublist.sublist("Time Integration Mode").sublist("Initialize To Steady").get<double>("Switch");
        double end_time = exe_sublist.sublist("Time Integration Mode").sublist("Initialize To Steady").get<double>("End");

        time_map.erase(start_time);
        time_map.erase(switch_time);
        time_map.erase(end_time);
      }
      if (exe_sublist.sublist("Time Integration Mode").isSublist("Steady")) {
        double start_time = exe_sublist.sublist("Time Integration Mode").sublist("Steady").get<double>("Start");
        double end_time = exe_sublist.sublist("Time Integration Mode").sublist("Steady").get<double>("End");

        time_map.erase(start_time);
        time_map.erase(end_time);
      }
      if (exe_sublist.sublist("Time Integration Mode").isSublist("Transient")) {
        double start_time = exe_sublist.sublist("Time Integration Mode").sublist("Transient").get<double>("Start");
        double end_time = exe_sublist.sublist("Time Integration Mode").sublist("Transient").get<double>("End");

        time_map.erase(start_time);
        time_map.erase(end_time);
        
        /* EIB: proposed v1.2.2 update - Add Maximum Cycle Number for Transient modes */
        if (exe_sublist.sublist("Time Integration Mode").sublist("Transient").isParameter("Maximum Cycle Number"))
          max_cycle_number =exe_sublist.sublist("Time Integration Mode").sublist("Transient").get<double>("Maximum Cycle Number");
      }
      if (exe_sublist.sublist("Time Integration Mode").isSublist("Transient with Static Flow")) {
        double start_time = exe_sublist.sublist("Time Integration Mode").sublist("Transient with Static Flow").get<double>("Start");
        double end_time = exe_sublist.sublist("Time Integration Mode").sublist("Transient with Static Flow").get<double>("End");

        time_map.erase(start_time);
        time_map.erase(end_time);
        
        /* EIB: proposed v1.2.2 update - Add Maximum Cycle Number for Transient modes */
        if (exe_sublist.sublist("Time Integration Mode").sublist("Transient").isParameter("Maximum Cycle Number"))
            max_cycle_number =exe_sublist.sublist("Time Integration Mode").sublist("Transient").get<double>("Maximum Cycle Number");
      }	
    }
  }

  /* TODO: need to do something with the new Maximum Cycle Number parameter in unstructured */
  
  start_times.clear();
  initial_time_step.clear();

  for (std::map<double,double>::const_iterator map_it = time_map.begin();
       map_it != time_map.end(); ++map_it) {
    start_times.push_back(map_it->first);
    initial_time_step.push_back(map_it->second);
  }

  ASSERT(start_times.size() == initial_time_step.size());

  tpc_list.set<Teuchos::Array<double> >("Start Times", start_times);
  tpc_list.set<Teuchos::Array<double> >("Initial Time Step", initial_time_step);

  return tpc_list;
}


/* ******************************************************************
* Empty
****************************************************************** */
Teuchos::ParameterList InputParserIS::CreateMPC_List_(Teuchos::ParameterList* plist)
{
  Teuchos::ParameterList mpc_list;

  bool transport_on(false), chemistry_on(false), flow_on(false);

  if (plist->isSublist("Execution Control")) {
    Teuchos::ParameterList exe_sublist = plist->sublist("Execution Control");

    mpc_list.sublist("Time Integration Mode") = exe_sublist.sublist("Time Integration Mode");

    //if (exe_sublist.isSublist("Time Period Control")) {
    mpc_list.sublist("Time Period Control") = CreateTimePeriodControlList_(plist);
    //}

    // now interpret the modes
    if (exe_sublist.isParameter("Transport Model")) {
      if (exe_sublist.get<std::string>("Transport Model") == "Off" || exe_sublist.get<std::string>("Transport Model") == "off") {
        mpc_list.set<std::string>("disable Transport_PK","yes");
      } else if (exe_sublist.get<std::string>("Transport Model") == "On" || exe_sublist.get<std::string>("Transport Model") == "on") {
        transport_on = true;
        mpc_list.set<std::string>("disable Transport_PK","no");

      } else {
        Exceptions::amanzi_throw(Errors::Message("Transport Model must either be On or Off"));
      }
    }
    else {
      Exceptions::amanzi_throw(Errors::Message("The parameter Transport Model must be specified."));
    }

    if (exe_sublist.isParameter("Chemistry Model")) {
      if (exe_sublist.get<std::string>("Chemistry Model") == "Off") {
        mpc_list.set<std::string>("Chemistry Model","Off");
      } else {
        chemistry_on = true;
        std::string chem_model = exe_sublist.get<std::string>("Chemistry Model");
        mpc_list.set<std::string>("Chemistry Model", chem_model);
      }
    }
    else {
      Exceptions::amanzi_throw(Errors::Message("The parameter \'Chemistry Model\' must be specified."));
    }

    // set defaults
    if (transport_on) {
      mpc_list.set<bool>("transport subcycling", TRANSPORT_SUBCYCLING);
    }
    if (transport_on && chemistry_on) {
      mpc_list.set<double>("max chemistry to transport timestep ratio",CHEM_TRANS_DT_RATIO);
    }

    if (transport_on) {
      if (plist->sublist("Execution Control").isSublist("Numerical Control Parameters")) {
        if (plist->sublist("Execution Control").sublist("Numerical Control Parameters").isSublist("Unstructured Algorithm")) {
          Teuchos::ParameterList& ncpu_list = plist->sublist("Execution Control").sublist("Numerical Control Parameters").sublist("Unstructured Algorithm");

          if (ncpu_list.isSublist("Transport Process Kernel")) {
            Teuchos::ParameterList& ncp_list = ncpu_list.sublist("Transport Process Kernel");

            if (ncp_list.isParameter("transport subcycling")) {
              mpc_list.set<bool>("transport subcycling",ncp_list.get<bool>("transport subcycling"));
            }
          }
        }
      }
    }

    if (transport_on && chemistry_on) {
      if (plist->sublist("Execution Control").isSublist("Numerical Control Parameters")) {
        if (plist->sublist("Execution Control").sublist("Numerical Control Parameters").isSublist("Unstructured Algorithm")) {
          Teuchos::ParameterList& ncpu_list = plist->sublist("Execution Control").sublist("Numerical Control Parameters").sublist("Unstructured Algorithm");
          if (ncpu_list.isSublist("Chemistry Process Kernel")) {
            Teuchos::ParameterList& ncp_list = ncpu_list.sublist("Chemistry Process Kernel");

            if (ncp_list.isParameter("max chemistry to transport timestep ratio")) {
              mpc_list.set<double>("max chemistry to transport timestep ratio",ncp_list.get<double>("max chemistry to transport timestep ratio"));
            }
          }
        }
      }
    }

    if (transport_on || chemistry_on) {
      mpc_list.set<Teuchos::Array<std::string> >("component names", comp_names_all_);
    }

    if (exe_sublist.isParameter("Flow Model")) {
      if (exe_sublist.get<std::string>("Flow Model") == "Off" || exe_sublist.get<std::string>("Flow Model") == "off") {
        mpc_list.set<std::string>("disable Flow_PK", "yes");
        flow_on = false;
      } else if (exe_sublist.get<std::string>("Flow Model") == "Richards") {
        mpc_list.set<std::string>("disable Flow_PK", "no");
        mpc_list.set<std::string>("Flow model","Richards");
        flow_on = true;
      } else if (exe_sublist.get<std::string>("Flow Model") == "Single Phase") {
        mpc_list.set<std::string>("disable Flow_PK", "no");
        mpc_list.set<std::string>("Flow model","Steady State Saturated");
        flow_on = true;
      } else {
        Exceptions::amanzi_throw(Errors::Message("Flow Model must either be Richards, Single Phase, or Off"));
      }
    }
    else {
      Exceptions::amanzi_throw(Errors::Message("The parameter Flow Model must be specified."));
    }

    if (flow_on) {
      double ti_rescue(TI_RESCUE_REDUCTION_FACTOR);
      Teuchos::ParameterList& exe_list = plist->sublist("Execution Control");

      if (exe_list.isSublist("Numerical Control Parameters")) {
        if (exe_list.sublist("Numerical Control Parameters").isSublist("Unstructured Algorithm")) {
          Teuchos::ParameterList& ncpu_list = exe_list.sublist("Numerical Control Parameters").sublist("Unstructured Algorithm");
          if (ncpu_list.isSublist("MPC")) {
            ti_rescue = ncpu_list.sublist("MPC").get<double>("time integration rescue reduction factor", TI_RESCUE_REDUCTION_FACTOR);
          }
        }
      }
      mpc_list.set<double>("time integration rescue reduction factor", TI_RESCUE_REDUCTION_FACTOR);

      // Add Picard option
      use_picard_ = USE_PICARD;
      if (exe_list.sublist("Numerical Control Parameters").sublist("Unstructured Algorithm")
          .isSublist("Flow Process Kernel")) {
        Teuchos::ParameterList fpk_params = exe_list.sublist("Numerical Control Parameters")
                                                    .sublist("Unstructured Algorithm")
                                                    .sublist("Flow Process Kernel");
        if (fpk_params.isParameter("Use Picard")) {
          use_picard_ = fpk_params.get<bool>("Use Picard", USE_PICARD);
        }
      }

      if (use_picard_) {
        Teuchos::ParameterList& tim_list = mpc_list.sublist("Time Integration Mode");
        if (tim_list.isSublist("Initialize To Steady")) {
          tim_list.sublist("Initialize To Steady").set<bool>("Use Picard", use_picard_);
        } else if (tim_list.isSublist("Steady")) {
          tim_list.sublist("Steady").set<bool>("Use Picard", use_picard_);
        }
      }
    }

    /* EIB: proposed v1.2.2 update - Change Restart name */
    //if (plist->sublist("Execution Control").isSublist("Restart from Checkpoint Data File") &&
    if (plist->sublist("Execution Control").isSublist("Restart") &&
        plist->sublist("Initial Conditions").isParameter("Init from Checkpoint File")) {
      // this is an error, you can either restart or re-init, but not both
      Exceptions::amanzi_throw(Errors::Message("You can either restart from a checkpoint or initialize from a checkpoint, but not both."));
    }
    
    /* EIB: proposed v1.2.2 update - Change Restart name */
    /*
    if (plist->sublist("Execution Control").isSublist("Restart from Checkpoint Data File")) {
      mpc_list.sublist("Restart from Checkpoint Data File") =
      plist->sublist("Execution Control").sublist("Restart from Checkpoint Data File");
    }
     */
    if (plist->sublist("Execution Control").isSublist("Restart")) {
      mpc_list.sublist("Restart from Checkpoint Data File") =
          plist->sublist("Execution Control").sublist("Restart");
    }
    
    /* EIB: proposed v1.2.2 update - Change Restart name */
    if (plist->sublist("Initial Conditions").isParameter("Init from Checkpoint File")) {
      //Teuchos::ParameterList& rest_list = mpc_list.sublist("Restart from Checkpoint Data File");
      Teuchos::ParameterList& rest_list = mpc_list.sublist("Restart");
      std::string file = plist->sublist("Initial Conditions").get<std::string>("Init from Checkpoint File");
      //rest_list.set<std::string>("Checkpoint Data File Name",file);
      rest_list.set<std::string>("File Name",file);
      rest_list.set<bool>("initialize from checkpoint data file and do not restart",true);
    }
  }

  mpc_list.sublist("VerboseObject") = CreateVerbosityList_(verbosity_level);

  return mpc_list;
}

}  // namespace AmanziInput
}  // namespace Amanzi
