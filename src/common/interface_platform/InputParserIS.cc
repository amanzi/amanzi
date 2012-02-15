#include "InputParserIS.hh"

#include "Teuchos_XMLParameterListHelpers.hpp"

#include <string>

#include "errors.hh"
#include "exceptions.hh"
#include "dbc.hh"


namespace Amanzi {
namespace AmanziInput {

Teuchos::ParameterList translate(Teuchos::ParameterList* plist, int numproc) 
{
  numproc_ = numproc;

  Teuchos::ParameterList new_list, tmp_list;

  init_global_info(plist);

  // unstructured header
  new_list.set<bool>("Native Unstructured Input", true);
  new_list.set<std::string>("grid_option", "Unstructured");

  // checkpoint list is optional
  tmp_list = create_Checkpoint_Data_List (plist);
  if (tmp_list.begin() != tmp_list.end()) {
    new_list.sublist("Checkpoint Data") = tmp_list;
  }

  tmp_list = create_Visualization_Data_List (plist);
  if (tmp_list.begin() != tmp_list.end()) {
    new_list.sublist("Visualization Data") = tmp_list;
  }

  if (plist->sublist("Output").isSublist("Observation Data")) {
    Teuchos::ParameterList& od_list = plist->sublist("Output").sublist("Observation Data");
    if (od_list.begin() != od_list.end()) {
      new_list.sublist("Observation Data") = create_Observation_Data_List (plist);
    }
  }

  new_list.sublist("Regions") = get_Regions_List(plist);
  new_list.sublist("Mesh") = translate_Mesh_List(plist);
  new_list.sublist("Domain") = get_Domain_List(plist);
  new_list.sublist("MPC") = create_MPC_List(plist);
  new_list.sublist("Transport") = create_Transport_List(plist);
  new_list.sublist("State") = create_State_List(plist);
  new_list.sublist("Flow") = create_Flow_List(plist);

  return new_list;
}



Teuchos::ParameterList get_Time_Macro (const std::string& macro_name, Teuchos::ParameterList* plist) 
{
  Teuchos::ParameterList time_macro;

  if ( plist->sublist("Output").sublist("Time Macros").isSublist(macro_name) ) {
    if (plist->sublist("Output").sublist("Time Macros").sublist(macro_name).isParameter("Start_Period_Stop")) {
      Teuchos::Array<double> time_range;
      time_range = plist->sublist("Output").sublist("Time Macros").sublist(macro_name)
          .get<Teuchos::Array<double> >("Start_Period_Stop");

      time_macro.set<Teuchos::Array<double> >("Start_Period_Stop",time_range);

    }
    if (plist->sublist("Output").sublist("Time Macros").sublist(macro_name).isParameter("Values")) {
      Teuchos::Array<double> values;
      values = plist->sublist("Output").sublist("Time Macros").sublist(macro_name)
          .get<Teuchos::Array<double> >("Values");
      time_macro.set<Teuchos::Array<double> >("Values",values);
    }
  } else {
    std::stringstream ss;
    ss << "The time macro " << macro_name << " does not exist in the input file";
    Exceptions::amanzi_throw(Errors::Message(ss.str().c_str()));
  }

  return time_macro;
}


Teuchos::Array<int> get_Cycle_Macro ( const std::string& macro_name, Teuchos::ParameterList* plist ) 
{
  Teuchos::Array<int> cycle_range;

  if ( plist->sublist("Output").sublist("Cycle Macros").isSublist(macro_name) ) {
    cycle_range = plist->sublist("Output").sublist("Cycle Macros").sublist(macro_name)
        .get<Teuchos::Array<int> >("Start_Period_Stop");
  } else {
    std::stringstream ss;
    ss << "The cycle macro " << macro_name << " does not exist in the input file";
    Exceptions::amanzi_throw(Errors::Message(ss.str().c_str()));
  }

  return cycle_range;
}


Teuchos::Array<std::string> get_Variable_Macro ( Teuchos::Array<std::string>& macro_name, Teuchos::ParameterList* plist ) 
{
  std::vector<std::string> vars;

  for (int i=0; i<macro_name.size(); i++) {
    if ( plist->sublist("Output").sublist("Variable Macros").isSublist(macro_name[i]) ) {
      Teuchos::ParameterList& macro_list = plist->sublist("Output").sublist("Variable Macros").sublist(macro_name[i]);

      if ( macro_list.isParameter("Phase") ) {
        std::string macro_phase = macro_list.get<std::string>("Phase");
        if (macro_phase == "All") {
          vars.push_back(phase_comp_name);
        } else {  // not All, must equal phase_comp_name
          if ( macro_list.isParameter("Component") ) {
            std::string macro_comp = macro_list.get<std::string>("Component");
            if (macro_comp == "All") {
              vars.push_back(phase_comp_name);
            } else { // not All, must equal
              if ( macro_comp != phase_comp_name ) {
                std::stringstream ss;
                ss << "The phase component name " << macro_comp << " is refered to in a variable macro but is not defined";
                Exceptions::amanzi_throw(Errors::Message(ss.str().c_str()));
              }
              vars.push_back(macro_comp);
            }

          }
        }
      }

      if ( macro_list.isParameter("Solute") ) {
        std::string macro_solute = macro_list.get<std::string>("Solute");
        if ( macro_solute == "All" ) {
          for ( int i=0; i<comp_names.size(); i++) vars.push_back(comp_names[i]);
        } else {
          vars.push_back(macro_solute);
        }
      }
    } else {
      std::stringstream ss;
      ss << "The variable macro " << macro_name[i] << " does not exist in the input file";
      Exceptions::amanzi_throw(Errors::Message(ss.str().c_str()));
    }
  }

  Teuchos::Array<std::string> ret_vars(vars.size());

  for (int i=0; i<vars.size(); i++) {
    ret_vars[i] = vars[i];
  }

  return ret_vars;
}


void init_global_info( Teuchos::ParameterList* plist ) 
{
  Teuchos::ParameterList& phase_list = plist->sublist("Phase Definitions");
  if ( (++ phase_list.begin()) == phase_list.end() ) {
    phase_name = phase_list.name(phase_list.begin());
    phase_comp_name = phase_list.sublist(phase_name).sublist("Phase Components").name(phase_list.sublist(phase_name).sublist("Phase Components").begin());

    comp_names = phase_list.sublist(phase_name).sublist("Phase Components").sublist(phase_comp_name).get<Teuchos::Array<std::string> >("Component Solutes");

    // create a map for the components
    for (int i=0; i<comp_names.size(); i++) {
      comp_names_map[comp_names[i]] = i;
    }
  } else {
    // we can only do single phase
  }

  if ( plist->isSublist("Execution Control") ) {

    if ( plist->sublist("Execution Control").isParameter("Verbosity") ) {
      std::string verbosity = plist->sublist("Execution Control").get<std::string>("Verbosity");

      if ( verbosity == "None" ) {
        verbosity_level = "none";
      } else if ( verbosity == "Low" ) {
        verbosity_level = "low";
      } else if ( verbosity == "Medium" ) {
        verbosity_level = "medium";
      } else if ( verbosity == "High" ) {
        verbosity_level = "high";
      } else if ( verbosity == "Extreme" ) {
        verbosity_level = "high";
      } else {
        Exceptions::amanzi_throw(Errors::Message("Verbosity must be one of None, Low, Medium, High, or Extreme."));
      }

    } else {
      verbosity_level = "low";
    }
  }
}


Teuchos::ParameterList create_Checkpoint_Data_List ( Teuchos::ParameterList* plist ) 
{
  Teuchos::ParameterList restart_list;

  if ( plist->isSublist("Output") ) {

    if ( plist->sublist("Output").isSublist("Checkpoint Data") ) {
      restart_list = plist->sublist("Output").sublist("Checkpoint Data");

      // check if the cycle range is defined via a macro
      if ( restart_list.isParameter("Cycle Macro") ) {
        std::string cycle_macro = restart_list.get<std::string>("Cycle Macro");

        Teuchos::Array<int> range = get_Cycle_Macro(cycle_macro, plist);
        Teuchos::ParameterList& c_restart_list = restart_list.sublist("Cycle Data");

        c_restart_list.set<int>("Start", range[0]);
        c_restart_list.set<int>("End", range[2]);
        c_restart_list.set<int>("Interval", range[1]);
        // now delete the Cycle Macro paramter

        restart_list.remove("Cycle Macro");
      }
    }
  }

  return restart_list;
}



Teuchos::ParameterList create_Visualization_Data_List ( Teuchos::ParameterList* plist ) 
{
  Teuchos::ParameterList vis_list;

  if ( plist->isSublist("Output") ) {

    if ( plist->sublist("Output").isSublist("Visualization Data") ) {
      vis_list = plist->sublist("Output").sublist("Visualization Data");

      // check if the cycle range is defined via a macro
      if ( vis_list.isParameter("Cycle Macro") ) {
        std::string cycle_macro = vis_list.get<std::string>("Cycle Macro");
        Teuchos::ParameterList &cdata = vis_list.sublist("Cycle Data");
        Teuchos::Array<int> cm = get_Cycle_Macro(cycle_macro,plist);

        cdata.set("Start",cm[0]);
        cdata.set("End",cm[2]);
        cdata.set("Interval",cm[1]);

        // delete the cycle macro
        vis_list.remove("Cycle Macro");
      }
    }
  }
  return vis_list;
}


Teuchos::ParameterList create_Observation_Data_List ( Teuchos::ParameterList* plist ) 
{
  Teuchos::ParameterList obs_list;

  if ( plist->isSublist("Output") ) {
    if ( plist->sublist("Output").isSublist("Observation Data") ) {

      Teuchos::ParameterList olist = plist->sublist("Output").sublist("Observation Data");

      if (olist.isParameter("Observation Output Filename")) {
        obs_list.set<std::string>("Observation Output Filename",olist.get<std::string>("Observation Output Filename"));
      } else {
        Exceptions::amanzi_throw(Errors::Message("The required parameter Observation Output Filename was not specified."));
      }

      for ( Teuchos::ParameterList::ConstIterator i = olist.begin();
            i != olist.end(); i++ ) {
        if (  olist.isSublist( i->first ) ) {
          // copy the observation data sublist
          obs_list.sublist(i->first) = olist.sublist(i->first);

          if ( obs_list.sublist(i->first).isParameter("Time Macro") ) {
            std::string time_macro = obs_list.sublist(i->first).get<std::string>("Time Macro");
            Teuchos::ParameterList time_macro_list = get_Time_Macro(time_macro, plist);
            if (time_macro_list.isParameter("Start_Period_Stop")) {
              obs_list.sublist(i->first).set("Start_Period_Stop",time_macro_list.get<Teuchos::Array<double> >("Start_Period_Stop"));
            }
            if (time_macro_list.isParameter("Values")) {
              obs_list.sublist(i->first).set("Values",time_macro_list.get<Teuchos::Array<double> >("Values"));
            }
            obs_list.sublist(i->first).remove("Time Macro");
          }

          if ( obs_list.sublist(i->first).isParameter("Cycle Macro") ) {
            std::string cycle_macro = obs_list.sublist(i->first).get<std::string>("Cycle Macro");
            obs_list.sublist(i->first).set("Start_Period_Stop",get_Cycle_Macro(cycle_macro, plist));
            obs_list.sublist(i->first).remove("Cycle Macro");
          }

          // if ( obs_list.sublist(i->first).isParameter("Variable Macro") ) {
          //   Teuchos::Array<std::string> var_macro = obs_list.sublist(i->first).get<Teuchos::Array<std::string> >("Variable Macro");
          //   obs_list.sublist(i->first).set("Variables",  get_Variable_Macro(var_macro, plist));
          //   obs_list.sublist(i->first).remove("Variable Macro");
          // }
        }
      }
    }
  }
  return obs_list;
}



Teuchos::ParameterList get_Regions_List ( Teuchos::ParameterList* plist ) 
{
  Teuchos::ParameterList reg_list;

  if ( plist->isSublist("Regions") ) {
    reg_list = plist->sublist("Regions");
  }

  return reg_list;
}


Teuchos::ParameterList get_Mesh_List ( Teuchos::ParameterList* plist ) 
{
  Teuchos::ParameterList msh_list;

  if ( plist->isSublist("Mesh") ) {
    msh_list = plist->sublist("Mesh");
  }

  return msh_list;
}


Teuchos::ParameterList translate_Mesh_List ( Teuchos::ParameterList* plist ) 
{
  Teuchos::ParameterList msh_list;

  if ( plist->isSublist("Mesh") ) {
    if (plist->sublist("Mesh").isSublist("Unstructured") ) {
      if (plist->sublist("Mesh").sublist("Unstructured").isSublist("Generate Mesh")) {
        Teuchos::ParameterList& generate = plist->sublist("Mesh").sublist("Unstructured").sublist("Generate Mesh").sublist("Uniform Structured");
        Teuchos::Array<int> ncells = generate.get<Teuchos::Array<int> >("Number of Cells");
        Teuchos::Array<double> low = generate.get<Teuchos::Array<double> >("Domain Low Corner");
        Teuchos::Array<double> high = generate.get<Teuchos::Array<double> >("Domain High Corner");

        // msh_list.set<std::string>("Framework","MSTK");
        Teuchos::ParameterList& msh_gen = msh_list.sublist("Unstructured").sublist("Generate Mesh");

        msh_gen.set<int>("Number of Cells in X",ncells[0]);
        msh_gen.set<int>("Number of Cells in Y",ncells[1]);
        msh_gen.set<int>("Number of Cells in Z",ncells[2]);

        msh_gen.set<double>("X_Min",low[0]);
        msh_gen.set<double>("X_Max",high[0]);

        msh_gen.set<double>("Y_Min",low[1]);
        msh_gen.set<double>("Y_Max",high[1]);

        msh_gen.set<double>("Z_Min",low[2]);
        msh_gen.set<double>("Z_Max",high[2]);

      } else if (plist->sublist("Mesh").sublist("Unstructured").isSublist("Read Mesh File")) {
        std::string format = plist->sublist("Mesh").sublist("Unstructured").sublist("Read Mesh File").get<std::string>("Format");
        if ( format == "Exodus II") {
          // process the file name to replace .exo with .par in the case of a parallel run
          Teuchos::ParameterList& fn_list =  msh_list.sublist("Unstructured").sublist("Read Mesh File");
          fn_list.set<std::string>("Format","Exodus II");
          std::string file =  plist->sublist("Mesh").sublist("Unstructured").sublist("Read Mesh File").get<std::string>("File");
          std::string suffix(file.substr(file.size()-4,4));

          if ( suffix != ".exo" ) {
            // exodus files must have the .exo suffix
            Exceptions::amanzi_throw(Errors::Message("Exodus II was specified as a mesh file format but the suffix of the file that was specified is not .exo"));
          }

          // figure out if this is a parallel run
          if (numproc_ > 1) {
            std::string par_file = file.replace(file.size()-4,4,std::string(".par"));

            fn_list.set<std::string>("File",par_file);
          } else {
            // don't translate the suffix if this is a serial run
            fn_list.set<std::string>("File",file);
          }

          msh_list.sublist("Unstructured").sublist("Read Mesh File") = fn_list;

        } else {
          msh_list.sublist("Unstructured").sublist("Read Mesh File") = plist->sublist("Mesh").sublist("Unstructured").sublist("Read Mesh File");
        }
      }

      if (plist->sublist("Mesh").sublist("Unstructured").isSublist("Expert")) {
        msh_list.sublist("Unstructured").sublist("Expert") = plist->sublist("Mesh").sublist("Unstructured").sublist("Expert");
      }
    }
  }

  return msh_list;
}



Teuchos::ParameterList get_Domain_List ( Teuchos::ParameterList* plist ) 
{
  Teuchos::ParameterList dom_list;

  if ( plist->isSublist("Domain") ) {
    dom_list = plist->sublist("Domain");
  }

  return dom_list;
}



Teuchos::ParameterList create_MPC_List ( Teuchos::ParameterList* plist ) 
{
  Teuchos::ParameterList mpc_list;

  if ( plist->isSublist("Execution Control") ) {
    Teuchos::ParameterList exe_sublist = plist->sublist("Execution Control");

    mpc_list.sublist("Time Integration Mode") = exe_sublist.sublist("Time Integration Mode");

    // now interpret the modes
    if ( exe_sublist.isParameter("Transport Model") ) {
      if ( exe_sublist.get<std::string>("Transport Model") == "Off" ) {
        mpc_list.set<std::string>("disable Transport_PK","yes");
      } else if ( exe_sublist.get<std::string>("Transport Model") == "On" ) {
        mpc_list.set<std::string>("disable Transport_PK","no");
      } else {
        Exceptions::amanzi_throw(Errors::Message("Transport Model must either be On or Off"));
      }
    } else {
      Exceptions::amanzi_throw(Errors::Message("The parameter Transport Model must be specified."));
    }

    if ( exe_sublist.isParameter("Flow Model") ) {
      if ( exe_sublist.get<std::string>("Flow Model") == "Off" ) {
        mpc_list.set<std::string>("disable Flow_PK", "yes");
      } else if ( exe_sublist.get<std::string>("Flow Model") == "Richards" ) {
        mpc_list.set<std::string>("disable Flow_PK", "no");
        mpc_list.set<std::string>("Flow model","Richards");
      } else {
        Exceptions::amanzi_throw(Errors::Message("Flow Model must either be Richards or Off"));
      }
    } else {
      Exceptions::amanzi_throw(Errors::Message("The parameter Flow Model must be specified."));
    }

    if ( exe_sublist.isParameter("Chemistry Model") ) {
      if ( exe_sublist.get<std::string>("Chemistry Model") == "Off" ) {
        mpc_list.set<std::string>("disable Chemistry_PK","yes");
      } else {
        Exceptions::amanzi_throw(Errors::Message("Chemistry Model must be Off, we currently do not support Chemistry through the inpur spec."));
      }
    } else {
      Exceptions::amanzi_throw(Errors::Message("The parameter Chemistry Model must be specified."));
    }


    // if ( plist->sublist("Execution control").isSublist("Restart from Checkpoint Data File") ) {
    //   mpc_list.sublist("Restart from Checkpoint Data File") =
    //       plist->sublist("Execution control").sublist("Restart from Checkpoint Data File");
    // }
  }

  mpc_list.sublist("VerboseObject") = create_Verbosity_List(verbosity_level);

  return mpc_list;
}


Teuchos::ParameterList create_Transport_List ( Teuchos::ParameterList* plist ) 
{
  Teuchos::ParameterList trp_list;

  if ( plist->isSublist("Execution Control") ) {
    if ( plist->sublist("Execution Control").isParameter("Transport Model") ) {
      if ( plist->sublist("Execution Control").get<std::string>("Transport Model") == "On" ) {
        if (plist->sublist("Execution Control").isSublist("Numerical Control Parameters")) {
          Teuchos::ParameterList& ncp_list = plist->sublist("Execution Control").sublist("Numerical Control Parameters");
          if (ncp_list.isParameter("Transport Integration Algorithm")) {
            std::string tia = ncp_list.get<std::string>("Transport Integration Algorithm");
            if ( tia == "Explicit First-Order" ) {
              trp_list.set<int>("discretization order",1);
            } else if ( tia == "Explicit Second-Order" ) {
              trp_list.set<int>("discretization order",2);
            }
          }
        } else {
          trp_list.set<int>("discretization order",2);
        }
      }

      // continue to set some reasonable defaults
      trp_list.set<std::string>("enable internal tests","no");
      trp_list.set<double>("CFL",1.0);
    }
  }

  // now generate the boundary conditions
  // loop over the boundary condition sublists and extract the relevant data

  int n_transport_bcs = 0;

  Teuchos::ParameterList& bc_sublist = plist->sublist("Boundary Conditions");

  for (Teuchos::ParameterList::ConstIterator i = bc_sublist.begin(); i != bc_sublist.end(); i++) {
    // only count sublists
    if (bc_sublist.isSublist(bc_sublist.name(i))) {
      if ( bc_sublist.sublist((bc_sublist.name(i))).isSublist("Solute BC"))
        n_transport_bcs++;
    }
  }

  if (n_transport_bcs >= 0) {
    Teuchos::ParameterList& tbc_list = trp_list.sublist("Transport BCs");

    Teuchos::ParameterList& phase_list = plist->sublist("Phase Definitions");

    int bc_counter = 0;
    if ( (++ phase_list.begin()) == phase_list.end() ) {
      Teuchos::ParameterList& bc_sublist = plist->sublist("Boundary Conditions");

      for (Teuchos::ParameterList::ConstIterator i = bc_sublist.begin(); i != bc_sublist.end(); i++) {
        // read the assigned regions
        Teuchos::Array<std::string> regs = bc_sublist.sublist(bc_sublist.name(i)).get<Teuchos::Array<std::string> >("Assigned Regions");

        // only count sublists
        if (bc_sublist.isSublist(bc_sublist.name(i))) {
          if ( bc_sublist.sublist((bc_sublist.name(i))).isSublist("Solute BC")) {
            // read the solute bc stuff
            Teuchos::ParameterList& solbc = bc_sublist.sublist((bc_sublist.name(i))).sublist("Solute BC");

            Teuchos::ParameterList& comps = bc_sublist.sublist((bc_sublist.name(i))).sublist("Solute BC").sublist(phase_name).sublist(phase_comp_name);

            for (Teuchos::Array<std::string>::const_iterator i = comp_names.begin();
                 i != comp_names.end(); i++) {
              if (  comps.isSublist(*i) ) {
                std::stringstream compss;
                compss << "Component " << comp_names_map[*i];

                // for now just read the first value from the
                if ( comps.sublist(*i).isSublist("BC: Uniform Concentration") ) {
                  std::stringstream ss;
                  ss << "BC " << bc_counter;
                  Teuchos::ParameterList& bc = tbc_list.sublist(ss.str());

                  Teuchos::ParameterList& bcsub = comps.sublist(*i).sublist("BC: Uniform Concentration");

                  Teuchos::Array<double> values = bcsub.get<Teuchos::Array<double> >("Values");
                  Teuchos::Array<double> times = bcsub.get<Teuchos::Array<double> >("Times");
                  Teuchos::Array<std::string> time_fns = bcsub.get<Teuchos::Array<std::string> >("Time Functions");
                  bc.set<Teuchos::Array<double> >(compss.str(), values );
                  bc.set<Teuchos::Array<double> >("Times", times);
                  bc.set<Teuchos::Array<std::string> >("Time Functions", time_fns);
                  bc.set<Teuchos::Array<std::string> >("Regions", regs);

                  bc_counter++;
                }
              }
            }
          }
        }
      }
      tbc_list.set<int>("number of BCs", bc_counter);
    } else {
      Exceptions::amanzi_throw(Errors::Message( "Unstructured Amanzi can only have one phase, but the input file specifies more than one."));
    }
  }

  return trp_list;
}


Teuchos::ParameterList create_Flow_List ( Teuchos::ParameterList* plist ) 
{
  Teuchos::ParameterList flw_list;

  if ( plist->isSublist("Execution Control") ) {
    if ( plist->sublist("Execution Control").isParameter("Flow Model") ) {
      if ( plist->sublist("Execution Control").get<std::string>("Flow Model") == "Steady" ) {

      } else if ( plist->sublist("Execution Control").get<std::string>("Flow Model") == "Richards" ) {
        Teuchos::ParameterList& richards_problem = flw_list.sublist("Richards Problem");
        richards_problem.set<bool>("Upwind relative permeability", true);
        // this one should come from the input file...
        richards_problem.set<double>("Atmospheric pressure", 101325.0);

        Teuchos::ParameterList& richards_model_evaluator = richards_problem.sublist("Richards model evaluator");
        // set some reasonable defaults...
        richards_model_evaluator.set<double>("Absolute error tolerance",1.0);
        richards_model_evaluator.set<double>("Relative error tolerance",0.0);
        richards_model_evaluator.sublist("VerboseObject") = create_Verbosity_List(verbosity_level);

        Teuchos::ParameterList& time_integrator = richards_problem.sublist("Time integrator");
        // set some reasonable defaults...
        time_integrator.set<int>("Nonlinear solver max iterations", 6);
        time_integrator.set<int>("NKA max vectors", 5);
        time_integrator.set<int>("Maximum number of BDF tries", 10);
        time_integrator.set<double>("Nonlinear solver tolerance", 0.01);
        time_integrator.set<double>("NKA drop tolerance", 5.0e-2);
        time_integrator.sublist("VerboseObject") = create_Verbosity_List(verbosity_level);


        // set paramters for the steady state time integrator
        Teuchos::ParameterList& steady_time_integrator = richards_problem.sublist("steady time integrator");
        if (plist->sublist("Execution Control").isSublist("Numerical Control Parameters")) {
          if (plist->sublist("Execution Control").sublist("Numerical Control Parameters").isSublist("Unstructured Algorithm")) {
            Teuchos::ParameterList& num_list = plist->sublist("Execution Control").sublist("Numerical Control Parameters").sublist("Unstructured Algorithm");
	    steady_time_integrator.set<int>("steady max iterations", num_list.get<int>("steady max iterations",10));
	    steady_time_integrator.set<int>("steady min iterations", num_list.get<int>("steady min iterations",5));
            steady_time_integrator.set<int>("steady limit iterations", num_list.get<int>("steady limit iterations",20));
            steady_time_integrator.set<double>("steady nonlinear tolerance", num_list.get<double>("steady nonlinear tolerance",1.0));
            steady_time_integrator.set<double>("steady time step reduction factor", num_list.get<double>("steady time step reduction factor",0.8));
            steady_time_integrator.set<double>("steady time step increase factor", num_list.get<double>("steady time step increase factor",1.2));
	    steady_time_integrator.set<double>("steady max time step", num_list.get<double>("steady max time step",1.0e+8));
	    steady_time_integrator.set<int>("steady max preconditioner lag iterations", num_list.get<int>("steady max preconditioner lag iterations",5));
          } else {
            // set some probably not so good defaults for the steady computation
            steady_time_integrator.set<int>("steady max iterations",10);
            steady_time_integrator.set<int>("steady min iterations",5);
            steady_time_integrator.set<int>("steady limit iterations",20);
            steady_time_integrator.set<double>("steady nonlinear tolerance",1.0);
            steady_time_integrator.set<double>("steady time step reduction factor",0.8);
            steady_time_integrator.set<double>("steady time step increase factor",1.2);
	    steady_time_integrator.set<double>("steady max time step", 1.0e+8);
	    steady_time_integrator.set<int>("steady max preconditioner lag iterations", 5);
          }
        } else {
          // set some probably not so good defaults for the steady computation
          steady_time_integrator.set<int>("steady max iterations",10);
          steady_time_integrator.set<int>("steady min iterations",5);
          steady_time_integrator.set<int>("steady limit iterations",20);
          steady_time_integrator.set<double>("steady nonlinear tolerance",1.0);
          steady_time_integrator.set<double>("steady time step reduction factor",0.8);
          steady_time_integrator.set<double>("steady time step increase factor",1.2);
	  steady_time_integrator.set<double>("steady max time step", 1.0e+8);
	  steady_time_integrator.set<int>("steady max preconditioner lag iterations", 5);
        }
        steady_time_integrator.sublist("VerboseObject") = create_Verbosity_List(verbosity_level);

        // insert the water retention models sublist
        Teuchos::ParameterList &water_retention_models = richards_problem.sublist("Water retention models");
        water_retention_models = create_WRM_List(plist);

        // insert the flow BC sublist
        Teuchos::ParameterList& flow_bc = richards_problem.sublist("boundary conditions");
        flow_bc = create_SS_FlowBC_List(plist);

        // insert the diffusion preconditioner sublist
        Teuchos::ParameterList &diffprecon = richards_problem.sublist("Diffusion Preconditioner");
        diffprecon = create_DPC_List(plist);
      } else {
        // something's wrong
      }
    }
  }

  return flw_list;
}


Teuchos::ParameterList create_WRM_List ( Teuchos::ParameterList* plist ) 
{
  Teuchos::ParameterList wrm_list;

  // loop through the material properties list and extract the water retention model info

  Teuchos::ParameterList& matprop_list = plist->sublist("Material Properties");

  int counter = 0;
  for (Teuchos::ParameterList::ConstIterator i = matprop_list.begin(); i != matprop_list.end(); i++) {
    // get the wrm parameters

    std::string rel_perm = matprop_list.sublist(i->first).sublist("Capillary Pressure: van Genuchten").get<std::string>("Relative Permeability");
    if (rel_perm != "Mualem") {
      std::stringstream ss;
      ss << "Currently we can only deal with Mualem as the relative permeability model";
      Exceptions::amanzi_throw(Errors::Message(ss.str().c_str()));
    }

    double alpha = matprop_list.sublist(i->first).sublist("Capillary Pressure: van Genuchten").get<double>("alpha");
    double Sr = matprop_list.sublist(i->first).sublist("Capillary Pressure: van Genuchten").get<double>("Sr");
    double m = matprop_list.sublist(i->first).sublist("Capillary Pressure: van Genuchten").get<double>("m");

    // now get the assigned regions
    Teuchos::Array<std::string> regions = matprop_list.sublist(i->first).get<Teuchos::Array<std::string> >("Assigned Regions");

    for (Teuchos::Array<std::string>::const_iterator i = regions.begin();
         i != regions.end(); i++) {
      std::stringstream ss;
      ss << "Water Retention Model for " << *i;

      Teuchos::ParameterList& wrm_sublist = wrm_list.sublist(ss.str());

      wrm_sublist.set<std::string>("Water retention model", "van Genuchten");
      wrm_sublist.set<std::string>("Region",*i);
      wrm_sublist.set<double>("van Genuchten m", m);
      wrm_sublist.set<double>("van Genuchten alpha",alpha);
      wrm_sublist.set<double>("van Genuchten residual saturation", Sr);
    }
  }

  return wrm_list;
}


Teuchos::ParameterList create_DPC_List ( Teuchos::ParameterList* plist ) 
{
  Teuchos::ParameterList dpc_list;

  Teuchos::ParameterList& ml_list = dpc_list.sublist("ML Parameters");
  ml_list.set<int>("ML output", 0);
  ml_list.set<int>("max levels", 40);
  ml_list.set<std::string>("prec type","MGV");
  ml_list.set<int>("cycle applications", 2);
  ml_list.set<std::string>("aggregation: type", "Uncoupled-MIS");
  ml_list.set<double>("aggregation: damping factor", 1.33333);
  ml_list.set<double>("aggregation: threshold", 0.0);
  ml_list.set<std::string>("eigen-analysis: type","cg");
  ml_list.set<int>("eigen-analysis: iterations", 10);
  ml_list.set<int>("smoother: sweeps", 3);
  ml_list.set<double>("smoother: damping factor", 1.0);
  ml_list.set<std::string>("smoother: pre or post", "both");
  ml_list.set<std::string>("smoother: type", "Jacobi");
  ml_list.set<double>("smoother: damping factor", 1.0);
  ml_list.set<std::string>("coarse: type", "Amesos-KLU");
  ml_list.set<int>("coarse: max size", 256);

  return dpc_list;
}


Teuchos::ParameterList create_SS_FlowBC_List ( Teuchos::ParameterList* plist ) 
{
  Teuchos::ParameterList ssf_list;

  Teuchos::ParameterList& bc_sublist = plist->sublist("Boundary Conditions");

  int bc_counter = 0;

  for (Teuchos::ParameterList::ConstIterator i = bc_sublist.begin(); i != bc_sublist.end(); i++) {
    // look at sublists
    if (bc_sublist.isSublist(bc_sublist.name(i))) {
      Teuchos::ParameterList& bc = bc_sublist.sublist(bc_sublist.name(i));

      // get the regions
      Teuchos::Array<std::string> regions = bc.get<Teuchos::Array<std::string> >("Assigned Regions");

      if ( bc.isSublist("BC:Zero Flow") ) {
        // this is the natural BC for flow and we need not list it explicitly
      }

      else if ( bc.isSublist("BC: Flux") ) {
        Teuchos::ParameterList& bc_flux = bc.sublist("BC: Flux");

        Teuchos::Array<double> times = bc_flux.get<Teuchos::Array<double> >("Times");
        Teuchos::Array<std::string> time_fns = bc_flux.get<Teuchos::Array<std::string> >("Time Functions");

        if (! (bc_flux.isParameter("Inward Mass Flux") || bc_flux.isParameter("Outward Mass Flux"))  )  {
          // we can only handle mass fluxes right now
          Exceptions::amanzi_throw(Errors::Message("In BC: Flux we can only handle Mass Flux"));
        }

        Teuchos::Array<double> flux;

        if (bc_flux.isParameter("Inward Mass Flux")) {
          flux = bc_flux.get<Teuchos::Array<double> >("Inward Mass Flux");
        } else if (bc_flux.isParameter("Outward Mass Flux")) {
          flux = bc_flux.get<Teuchos::Array<double> >("Outward Mass Flux");
        }

        if (bc_flux.isParameter("Inward Mass Flux")) {
          for (int i=0; i<flux.size(); i++) flux[i] = - flux[i];
        }

        std::stringstream ss;
        ss << "BC " << bc_counter++;


        Teuchos::ParameterList& tbc = ssf_list.sublist("mass flux").sublist(ss.str());
        tbc.set<Teuchos::Array<std::string> >("regions", regions );


        if ( times.size() == 1 ) {
          Teuchos::ParameterList& tbcs = tbc.sublist("outward mass flux").sublist("function-constant");
          tbcs.set<double>("value",flux[0]);
        } else {
          Teuchos::ParameterList& tbcs = tbc.sublist("outward mass flux").sublist("function-tabular");

          tbcs.set<Teuchos::Array<double> >("x values", times);
          tbcs.set<Teuchos::Array<double> >("y values", flux);

          std::vector<std::string> forms_(time_fns.size());

          for (int i=0; i<time_fns.size(); i++)
            if (time_fns[i] == "Linear") {
              forms_[i] = "linear";
            } else if (time_fns[i] == "Constant") {
              forms_[i] = "constant";
            } else {
              Exceptions::amanzi_throw(Errors::Message("In the definition of BCs: tabular function can only be Linear or Constant"));
            }

          Teuchos::Array<std::string> forms = forms_;
          tbcs.set<Teuchos::Array<std::string> >("forms", forms);
        }

      } else if ( bc.isSublist("BC: Uniform Pressure") ) {
        Teuchos::ParameterList& bc_dir = bc.sublist("BC: Uniform Pressure");

        Teuchos::Array<double>      times = bc_dir.get<Teuchos::Array<double> >("Times");
        Teuchos::Array<std::string> time_fns = bc_dir.get<Teuchos::Array<std::string> >("Time Functions");
        Teuchos::Array<double>      values = bc_dir.get<Teuchos::Array<double> >("Values");

        std::stringstream ss;
        ss << "BC " << bc_counter++;

        Teuchos::ParameterList& tbc = ssf_list.sublist("pressure").sublist(ss.str());
        tbc.set<Teuchos::Array<std::string> >("regions", regions );

        if ( times.size() == 1 ) {
          Teuchos::ParameterList& tbcs = tbc.sublist("boundary pressure").sublist("function-constant");
          tbcs.set<double>("value",values[0]);
        } else {
          Teuchos::ParameterList& tbcs = tbc.sublist("boundary pressure").sublist("function-tabular");

          tbcs.set<Teuchos::Array<double> >("x values", times);
          tbcs.set<Teuchos::Array<double> >("y values", values);

          std::vector<std::string> forms_(time_fns.size());

          for (int i=0; i<time_fns.size(); i++)
            if (time_fns[i] == "Linear") {
              forms_[i] = "linear";
            } else if (time_fns[i] == "Constant") {
              forms_[i] = "constant";
            } else {
              Exceptions::amanzi_throw(Errors::Message("Tabular function can only be Linear or Constant"));
            }
          Teuchos::Array<std::string> forms = forms_;
          tbcs.set<Teuchos::Array<std::string> >("forms", forms);
        }

      } else if (  bc.isSublist("BC: Hydrostatic") ) {
        Teuchos::ParameterList& bc_dir = bc.sublist("BC: Hydrostatic");

        Teuchos::Array<double>      times = bc_dir.get<Teuchos::Array<double> >("Times");
        Teuchos::Array<std::string> time_fns = bc_dir.get<Teuchos::Array<std::string> >("Time Functions");
        Teuchos::Array<double>      values = bc_dir.get<Teuchos::Array<double> >("Water Table Height");

        std::stringstream ss;
        ss << "BC " << bc_counter++;

        Teuchos::ParameterList& tbc = ssf_list.sublist("static head").sublist(ss.str());
        tbc.set<Teuchos::Array<std::string> >("regions", regions );


        if ( times.size() == 1 ) {
          Teuchos::ParameterList& tbcs = tbc.sublist("water table elevation").sublist("function-constant");
          tbcs.set<double>("value",values[0]);
        } else {
          Teuchos::ParameterList& tbcs = tbc.sublist("water table elevation").sublist("function-tabular");

          tbcs.set<Teuchos::Array<double> >("x values", times);
          tbcs.set<Teuchos::Array<double> >("y values", values);

          std::vector<std::string> forms_(time_fns.size());

          for (int i=0; i<time_fns.size(); i++)
            if (time_fns[i] == "Linear") {
              forms_[i] = "linear";
            } else if (time_fns[i] == "Constant") {
              forms_[i] = "constant";
            } else  {
              Exceptions::amanzi_throw(Errors::Message("Tabular function can only be Linear or Constant"));
            }
          Teuchos::Array<std::string> forms = forms_;
          tbcs.set<Teuchos::Array<std::string> >("forms", forms);
        }
      }

      // TODO...
      // add the rest of the boundary conditions
    }
  }

  return ssf_list;
}


Teuchos::ParameterList create_State_List ( Teuchos::ParameterList* plist ) 
{
  Teuchos::ParameterList stt_list;

  stt_list.set<double>("Gravity x", 0.0);
  stt_list.set<double>("Gravity y", 0.0);
  stt_list.set<double>("Gravity z", -9.81);

  // find the viscosity
  Teuchos::ParameterList& phase_list = plist->sublist("Phase Definitions");

  // make sure there is only one phase
  if ( (++ phase_list.begin()) == phase_list.end() ) {
    // write the array of component solutes
    stt_list.set<Teuchos::Array<std::string> >("Component Solutes", comp_names);
    stt_list.set<int>("Number of component concentrations", comp_names.size());

    double viscosity = phase_list.sublist(phase_name).sublist("Phase Properties").sublist("Viscosity: Uniform").get<double>("Viscosity");
    double density = phase_list.sublist(phase_name).sublist("Phase Properties").sublist("Density: Uniform").get<double>("Density");

    stt_list.set<double>("Constant viscosity", viscosity);
    stt_list.set<double>("Constant water density", density);

    std::map<std::string,int> region_to_matid;
    std::map<int,std::string> matid_to_material;

    int matid_ctr = 0;
    // loop over the material properties
    Teuchos::ParameterList& matprop_list = plist->sublist("Material Properties");
    for (Teuchos::ParameterList::ConstIterator i = matprop_list.begin(); i != matprop_list.end(); i++) {
      // get the regions
      Teuchos::Array<std::string> regions = matprop_list.sublist(matprop_list.name(i)).get<Teuchos::Array<std::string> >("Assigned Regions");

      // record the material ID for each region that this material occupies
      matid_ctr++;
      for (int ii=0; ii<regions.size(); ii++) {
        if (region_to_matid.find(regions[ii]) == region_to_matid.end()) {
          region_to_matid[regions[ii]] = matid_ctr;
          matid_to_material[matid_ctr] = matprop_list.name(i);
        } else {
          std::stringstream ss;
          ss << "There is more than one material assinged to region " << regions[ii] << ".";
          Exceptions::amanzi_throw(Errors::Message(ss.str().c_str()));
        }
      }

      double porosity = matprop_list.sublist(matprop_list.name(i)).sublist("Porosity: Uniform").get<double>("Value");
      double perm_vert, perm_horiz;

      if (matprop_list.sublist(matprop_list.name(i)).isSublist("Intrinsic Permeability: Uniform")) {
        perm_vert = matprop_list.sublist(matprop_list.name(i)).sublist("Intrinsic Permeability: Uniform").get<double>("Value");
        perm_horiz = perm_vert;
      } else if (matprop_list.sublist(matprop_list.name(i)).isSublist("Intrinsic Permeability: Anisotropic Uniform")) {
        perm_vert = matprop_list.sublist(matprop_list.name(i)).sublist("Intrinsic Permeability: Anisotropic Uniform").get<double>("Vertical");
        perm_horiz = matprop_list.sublist(matprop_list.name(i)).sublist("Intrinsic Permeability: Anisotropic Uniform").get<double>("Horizontal");
      } else {
        Exceptions::amanzi_throw(Errors::Message("Permeability can only be specified as Intrinsic Permeability: Uniform, or Intrinsic Permeability: Anisotropic Uniform."));
      }

      for (Teuchos::Array<std::string>::const_iterator i=regions.begin(); i!=regions.end(); i++) {
        std::stringstream sss;
        sss << "Mesh block " << *i;

        Teuchos::ParameterList& stt_mat = stt_list.sublist(sss.str());

        stt_mat.set<double>("Constant porosity", porosity);
        stt_mat.set<double>("Constant vertical permeability", perm_vert);
        stt_mat.set<double>("Constant horizontal permeability", perm_horiz);
        stt_mat.set<std::string>("Region", *i);

        // find the initial conditions for the current region
        Teuchos::ParameterList& ic_list = plist->sublist("Initial Conditions");
        Teuchos::ParameterList* ic_for_region = NULL;
        for (Teuchos::ParameterList::ConstIterator it = ic_list.begin(); it != ic_list.end(); it++) {
          if (ic_list.isSublist(it->first)) {
            Teuchos::Array<std::string> ass_regions = ic_list.sublist(it->first).get<Teuchos::Array<std::string> >("Assigned Regions");
            if (ass_regions.size() == 1 && ass_regions[0] == "All") {
              ic_for_region = &(ic_list.sublist(it->first));
            } else {
              // check if the current region is part of the current initial condition's assigned regions
              for (int ii=0; ii<ass_regions.size(); ii++) {
                if (ass_regions[ii] == *i) {
                  ic_for_region = &(ic_list.sublist(it->first));
                }
              }
            }
          } else {
            Exceptions::amanzi_throw(Errors::Message("The list Initial Conditions contains an entry that is not a ParamterList itself."));
          }
        }
        // make sure that we actually have found an IC list that defines initial conditions for the current region
        if (ic_for_region == NULL) {
          std::stringstream ss;
          ss << "There is no sublist of the Initial Conditions list for region " << *i << ".";
          Exceptions::amanzi_throw(Errors::Message(ss.str().c_str()));
        }
        // at this point ic_for_region is the list that defines the inital conditions for the current region

        // write the initial conditions for pressure
        if ( ic_for_region->isSublist("IC: Uniform Pressure")) {
          Teuchos::ParameterList& sublist = stt_mat.sublist("uniform pressure");
          sublist.set<double>("value",ic_for_region->sublist("IC: Uniform Pressure").get<double>("Value"));
        } else if (ic_for_region->isSublist("IC: Linear Pressure")) {
          Teuchos::ParameterList& sublist = stt_mat.sublist("linear pressure");
          sublist.set<Teuchos::Array<double> >("gradient", ic_for_region->sublist("IC: Linear Pressure").get<Teuchos::Array<double> >("Gradient Value"));
          sublist.set<Teuchos::Array<double> >("reference coordinate", ic_for_region->sublist("IC: Linear Pressure").get<Teuchos::Array<double> >("Reference Coordinate"));
          sublist.set<double>("reference value", ic_for_region->sublist("IC: Linear Pressure").get<double>("Reference Value"));
        } else {
          Exceptions::amanzi_throw(Errors::Message("An initial condition for pressure must be specified. It must either be IC: Uniform Pressure, or IC: Linear Pressure."));
        }

        // write the initial conditions for saturation, since this is not a primary variable, this is not required
        // so we don't throw if these initial conditions are missing
        if ( ic_for_region->isSublist("IC: Uniform Saturation") ) {
          Teuchos::ParameterList& sublist = stt_mat.sublist("uniform saturation");
          sublist.set<double>("value",ic_for_region->sublist("IC: Uniform Saturation").get<double>("Value"));
        } else if (ic_for_region->isSublist("IC: Linear Saturation")) {
          Teuchos::ParameterList& sublist = stt_mat.sublist("linear saturation");
        }

        // write the initial conditions for the solutes, note that we hardcode for there only being one phase, with one phase component
        for (int ii=0; ii<comp_names.size(); ii++) {
          if (! ic_for_region->sublist("Solute IC").sublist(phase_name).sublist(phase_comp_name).isSublist(comp_names[ii])) {
            std::stringstream ss;
            ss << "Initial condition for solute " << comp_names[ii] << " in region " << *i << " is missing.";
            Exceptions::amanzi_throw(Errors::Message(ss.str().c_str()));
          }

          double conc = ic_for_region->sublist("Solute IC").sublist(phase_name).sublist(phase_comp_name).sublist(comp_names[ii]).sublist("IC: Uniform Concentration").get<double>("Value");

          std::stringstream ss;
          ss << "Constant component concentration " << comp_names_map[ comp_names[ii] ];

          stt_mat.set<double>(ss.str(), conc);
        }
      }
    }

    // write the mapping between region name and material id
    // (here material ID is an atificial integer that is only used for visualization)
    Teuchos::Array<int> matids(region_to_matid.size());
    Teuchos::Array<std::string> regnames(region_to_matid.size());
    Teuchos::Array<std::string> matnames(matid_to_material.size());

    int ii=0;
    for (std::map<std::string,int>::const_iterator it = region_to_matid.begin(); it != region_to_matid.end(); it++) {
      matids[ii] = it->second;
      regnames[ii] = it->first;
      ii++;
    }

    for (int k=0; k<matnames.size(); k++) {
      matnames[k] = matid_to_material[k+1];
    }

    stt_list.set<Teuchos::Array<int> >("Region Name to Material ID Map (Material IDs)",matids);
    stt_list.set<Teuchos::Array<std::string> >("Region Name to Material ID Map (Region Names)",regnames);
    stt_list.set<Teuchos::Array<std::string> >("Material Names",matnames);

  } else {
    Exceptions::amanzi_throw("There is more than one phase, however, amanzi-u only supports one phase");
  }

  return stt_list;
}


Teuchos::ParameterList create_Verbosity_List ( const std::string& vlevel ) 
{
  Teuchos::ParameterList vlist;

  if (vlevel == "low") {
    vlist.set<std::string>("Verbosity Level","low");
  } else if (vlevel == "medium") {
    vlist.set<std::string>("Verbosity Level","medium");
  } else if (vlevel == "high") {
    vlist.set<std::string>("Verbosity Level","high");
  } else if (vlevel == "none") {
    vlist.set<std::string>("Verbosity Level","none");
  }

  return vlist;
}

}  // namespace AmanziInput
}  // namespace Amanzi
