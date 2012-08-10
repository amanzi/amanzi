#include "InputParserIS.hh"
#include "InputParserIS-defaults.hh"
#include "Teuchos_XMLParameterListHelpers.hpp"

#include <string>
#include <algorithm>
#include <boost/lambda/lambda.hpp>
#include <boost/bind.hpp>

#include "errors.hh"
#include "exceptions.hh"
#include "dbc.hh"

#define BOOST_FILESYSTEM_VERSION 2
#include "boost/filesystem/operations.hpp"
#include "boost/filesystem/path.hpp"


namespace Amanzi {
namespace AmanziInput {


/**
 * /author   Nathan Barnett
 * /fn       compareEpsilon
 * /brief    When handed an element in an array and an epsilon value,
 *           compares the absolute difference between the two values
 * /returns  boolean - true if out of epsilon value
 */
template <typename T>
bool compareEpsilon(T& first, T eps) {
  return fabs(first-*(&first-1))<eps;
}


/* ******************************************************************
 * Empty
 ****************************************************************** */
Teuchos::ParameterList translate(Teuchos::ParameterList* plist, int numproc) {
  numproc_ = numproc;

  Teuchos::ParameterList new_list, tmp_list;

  init_global_info(plist);

  // unstructured header
  new_list.set<bool>("Native Unstructured Input", true);
  new_list.set<std::string>("grid_option", "Unstructured");
  new_list.set<std::string>("input file name", 
      plist->get<std::string>("input file name", "unit_test.xml"));

  // checkpoint list is optional
  tmp_list = create_Checkpoint_Data_List(plist);
  if (tmp_list.begin() != tmp_list.end()) {
    new_list.sublist("Checkpoint Data") = tmp_list;
  }

  tmp_list = create_Visualization_Data_List(plist);
  if (tmp_list.begin() != tmp_list.end()) {
    new_list.sublist("Visualization Data") = tmp_list;
  }

  if (plist->sublist("Output").isSublist("Observation Data")) {
    Teuchos::ParameterList& od_list = plist->sublist("Output").sublist("Observation Data");
    if (od_list.begin() != od_list.end()) {
      new_list.sublist("Observation Data") = create_Observation_Data_List(plist);
    }
  }

  new_list.sublist("Regions") = get_Regions_List(plist);
  new_list.sublist("Mesh") = translate_Mesh_List(plist);
  new_list.sublist("Domain") = get_Domain_List(plist);
  new_list.sublist("MPC") = create_MPC_List(plist);
  new_list.sublist("Transport") = create_Transport_List(plist);
  new_list.sublist("State") = create_State_List(plist);
  new_list.sublist("Flow") = create_Flow_List(plist);
  new_list.sublist("Preconditioners") = create_Preconditioners_List(plist);
  new_list.sublist("Solvers") = create_Solvers_List(plist);

  if (new_list.sublist("MPC").get<std::string>("Chemistry Model") != "Off") {
    new_list.sublist("Chemistry") = CreateChemistryList(plist);
  }

  return new_list;
}


/* ******************************************************************
 * Empty
 ****************************************************************** */
Teuchos::ParameterList get_Time_Macro(const std::string& macro_name, Teuchos::ParameterList* plist) {
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


/* ******************************************************************
 * Empty
 ****************************************************************** */
Teuchos::Array<int> get_Cycle_Macro(const std::string& macro_name, Teuchos::ParameterList* plist) {
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


/* ******************************************************************
 * Empty
 ****************************************************************** */
Teuchos::Array<std::string> get_Variable_Macro(Teuchos::Array<std::string>& macro_name, Teuchos::ParameterList* plist) {
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


/* ******************************************************************
 * Empty
 ****************************************************************** */
void init_global_info(Teuchos::ParameterList* plist) {

  // spatial dimension
  if ( plist->isSublist("Domain") ) {
    spatial_dimension_ = plist->sublist("Domain").get<int>("Spatial Dimension",0);
  } else {
    spatial_dimension_ = 0;
  }

  phase_name = "Aqueous";
  phase_comp_name = "Water";
  // don't know the history of these variables, clear them just to be safe.
  comp_names.clear();
  mineral_names_.clear();
  sorption_site_names_.clear();

  Teuchos::ParameterList& phase_list = plist->sublist("Phase Definitions");
  Teuchos::ParameterList::ConstIterator item;
  for (item = phase_list.begin(); item != phase_list.end(); ++item) {
    if (phase_list.name(item) == phase_name) {
      Teuchos::ParameterList aqueous_list = phase_list.sublist("Aqueous");
      if (aqueous_list.isSublist("Phase Components")) {
        Teuchos::ParameterList phase_components =
            aqueous_list.sublist("Phase Components");
        if (phase_components.isSublist("Water")) {
          Teuchos::ParameterList water_components =
              phase_components.sublist(phase_comp_name);
          comp_names =
              water_components.get<Teuchos::Array<std::string> >("Component Solutes");
        }  // end water
      }  // end phase components
    }  // end Aqueous phase
    else if (phase_list.name(item) == "Solid") {
      Teuchos::ParameterList solid_list = phase_list.sublist("Solid");
      // this is the order that the chemistry expects
      if (solid_list.isParameter("Minerals")) {
        mineral_names_ = solid_list.get<Teuchos::Array<std::string> >("Minerals");
      }
      if (solid_list.isParameter("Sorption Sites")) {
        sorption_site_names_ = solid_list.get<Teuchos::Array<std::string> >("Sorption Sites");
      }
    }  // end Solid phase
    else {
      std::stringstream message;
      message << "Error: InputParserIS::init_global_info(): "
              << "The only phases supported on unstructured meshes at this time are '"
              << phase_name << "' and 'Solid'!\n"
              << phase_list << std::endl;
      Exceptions::amanzi_throw(Errors::Message(message.str()));
    }
  }

  if (comp_names.size() > 0) {
    // create a map for the components
    for (int i=0; i<comp_names.size(); i++) {
      comp_names_map[comp_names[i]] = i;
    }
  } else {
    std::stringstream message;
    message << "Error: InputParserIS::init_global_info(): "
            << "component names must be defined in the phase definitions block!\n";
    Exceptions::amanzi_throw(Errors::Message(message.str()));
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


/* ******************************************************************
 * Empty
 ****************************************************************** */
Teuchos::ParameterList create_Checkpoint_Data_List(Teuchos::ParameterList* plist) {

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


/* ******************************************************************
 * Populate visualization list.
 ****************************************************************** */
Teuchos::ParameterList create_Visualization_Data_List(Teuchos::ParameterList* plist) {

  Teuchos::ParameterList  vis_list;
  Teuchos::Array<double>  visualizationPoints;

  if ( plist->isSublist("Output") ) {
    if ( plist->sublist("Output").isSublist("Visualization Data") ) {
      vis_list = plist->sublist("Output").sublist("Visualization Data");

      // Cycle Macro
      if ( vis_list.isParameter("Cycle Macro") ) {
        std::string cycle_macro = vis_list.get<std::string>("Cycle Macro");
        Teuchos::Array<int> cm = get_Cycle_Macro(cycle_macro,plist);

        vis_list.sublist("cycle start period stop").sublist(cycle_macro).set("start period stop",cm);

        // delete the cycle macro
        vis_list.remove("Cycle Macro");
      }

      // Time Macro
      if ( vis_list.isParameter("Time Macro") ) {
        std::vector<std::string> time_macros;
        time_macros = vis_list.get<Teuchos::Array<std::string> >("Time Macro").toVector();

        for (int i=0; i < time_macros.size(); i++) {
          // Create a local parameter list and store the time macro (3 doubles)
          Teuchos::ParameterList time_macro_list = get_Time_Macro(time_macros[i], plist);
          if (time_macro_list.isParameter("Start_Period_Stop")) {
            vis_list.sublist("time start period stop").sublist(time_macros[i]).set("start period stop",time_macro_list.get<Teuchos::Array<double> >("Start_Period_Stop"));
          }
          if (time_macro_list.isParameter("Values")) {
            vis_list.set("times",time_macro_list.get<Teuchos::Array<double> >("Values"));
          }
        }
        vis_list.remove("Time Macro");
      }
    }
  }
  return vis_list;
}




/* ******************************************************************
 * Empty
 ****************************************************************** */
Teuchos::ParameterList create_Observation_Data_List(Teuchos::ParameterList* plist) {
  using namespace boost;
  using boost::bind;

  // Create a parameter list for holding data
  Teuchos::ParameterList obs_list;
  Teuchos::Array<double> observationPoints;

  // Check if there is an "Output" XML node
  if (plist->isSublist("Output")) {
    // If "Output" exists, check if there is an "Observation Data" subnode
    if (plist->sublist("Output").isSublist("Observation Data")) {
      // If both exist, initialize a structure with the XML data
      Teuchos::ParameterList olist = plist->sublist("Output").sublist("Observation Data");
      // If the node has value refering to the name of the output file, grab it
      if (olist.isParameter("Observation Output Filename")) {
        obs_list.set<std::string>("Observation Output Filename", olist.get<std::string>("Observation Output Filename"));
      } else {
        Exceptions::amanzi_throw(Errors::Message("The required parameter Observation Output Filename was not specified."));
      }
      // Iterate through the array
      for (Teuchos::ParameterList::ConstIterator i = olist.begin(); i != olist.end(); i++) {
        // If the current iteration node is a "tree"
        if (olist.isSublist(i->first)) {
          // copy the observation data sublist into the local list
          obs_list.sublist(i->first) = olist.sublist(i->first);

          if (obs_list.sublist(i->first).isParameter("Time Macro")) {
            std::string time_macro = obs_list.sublist(i->first).get<std::string>("Time Macro");
            // Create a local parameter list and store the time macro (3 doubles)
            Teuchos::ParameterList time_macro_list = get_Time_Macro(time_macro, plist);
            if (time_macro_list.isParameter("Start_Period_Stop")) {
              obs_list.sublist(i->first).sublist("time start period stop").sublist(time_macro).set("start period stop", time_macro_list.get<Teuchos::Array<double> >("Start_Period_Stop"));
            }
            if (time_macro_list.isParameter("Values")) {
              obs_list.sublist(i->first).set("Values",time_macro_list.get<Teuchos::Array<double> >("Values"));
              Teuchos::Array<double> values = time_macro_list.get<Teuchos::Array<double> >("Values");
              obs_list.sublist(i->first).set<Teuchos::Array<double> >("times",values);
            }
            obs_list.sublist(i->first).remove("Time Macro");
          }
          if (obs_list.sublist(i->first).isParameter("Cycle Macro")) {
            std::string cycle_macro = obs_list.sublist(i->first).get<std::string>("Cycle Macro");
            obs_list.sublist(i->first).sublist("cycle start period stop").sublist(cycle_macro).set("start period stop", get_Cycle_Macro(cycle_macro, plist));
            obs_list.sublist(i->first).remove("Cycle Macro");
          }
        }
      }
    }
  }
  return obs_list;
}


/* ******************************************************************
 * Empty
 ****************************************************************** */
Teuchos::ParameterList get_Regions_List(Teuchos::ParameterList* plist) {
  Teuchos::ParameterList reg_list;
  if (plist->isSublist("Regions")) {
    // check labeled set consistency, file specified here must match file specified
    // in the Mesh list
    // first find the mesh file name specified in the Mesh list
    std::string meshfile;
    if (plist->sublist("Mesh").sublist("Unstructured").isSublist("Read Mesh File")) {
      meshfile = plist->sublist("Mesh").sublist("Unstructured").sublist("Read Mesh File").get<std::string>("File");
    }
    Teuchos::ParameterList& rlist = plist->sublist("Regions");
    // loop over all regions and find possible labeled set definitions
    for (Teuchos::ParameterList::ConstIterator i = rlist.begin(); i != rlist.end(); i++) {
      // only count sublists
      if (rlist.isSublist(rlist.name(i))) {
	if (rlist.sublist((rlist.name(i))).isSublist("Region: Labeled Set")) {
	  std::string file = rlist.sublist((rlist.name(i))).sublist("Region: Labeled Set").get<std::string>("File");
	  boost::filesystem::path meshfile_path(meshfile);
	  boost::filesystem::path labeled_set_meshfile_path(file);
	  if (!boost::filesystem::equivalent(meshfile_path, labeled_set_meshfile_path)) {
	    Errors::Message message("There is a labeled set region that refers to a mesh file that is different from the mesh file that is defined in the Mesh list: " + file);
	    Exceptions::amanzi_throw(message);
	  }
	}
      }   
    }
    // all is well, return the Regions list for insertion into the
    // native list
    reg_list = plist->sublist("Regions");
  }

  return reg_list;
}


/* ******************************************************************
 * Empty
 ****************************************************************** */
Teuchos::ParameterList get_Mesh_List(Teuchos::ParameterList* plist) {
  Teuchos::ParameterList msh_list;

  if (plist->isSublist("Mesh")) {
    msh_list = plist->sublist("Mesh");
  }

  return msh_list;
}


/* ******************************************************************
 * Empty
 ****************************************************************** */
Teuchos::ParameterList translate_Mesh_List(Teuchos::ParameterList* plist) {
  Teuchos::ParameterList msh_list;

  if ( plist->isSublist("Mesh") ) {
    if (plist->sublist("Mesh").isSublist("Unstructured") ) {
      if (plist->sublist("Mesh").sublist("Unstructured").isSublist("Generate Mesh")) {
        Teuchos::ParameterList& generate = plist->sublist("Mesh").sublist("Unstructured").sublist("Generate Mesh").sublist("Uniform Structured");
        Teuchos::Array<int> ncells = generate.get<Teuchos::Array<int> >("Number of Cells");
        Teuchos::Array<double> low = generate.get<Teuchos::Array<double> >("Domain Low Corner");
        Teuchos::Array<double> high = generate.get<Teuchos::Array<double> >("Domain High Corner");

        Teuchos::ParameterList& msh_gen = msh_list.sublist("Unstructured").sublist("Generate Mesh");

        msh_gen.set< Teuchos::Array<int> >("Number of Cells",ncells);
        msh_gen.set< Teuchos::Array<double> >("Domain Low Corner",low);
        msh_gen.set< Teuchos::Array<double> >("Domain High Corner",high);

      } else if (plist->sublist("Mesh").sublist("Unstructured").isSublist("Read Mesh File")) {
        std::string format = plist->sublist("Mesh").sublist("Unstructured").sublist("Read Mesh File").get<std::string>("Format");
        if (format == "Exodus II") {
          // process the file name to replace .exo with .par in the case of a parallel run
          Teuchos::ParameterList& fn_list =  msh_list.sublist("Unstructured").sublist("Read Mesh File");
          fn_list.set<std::string>("Format", "Exodus II");
          std::string file =  plist->sublist("Mesh").sublist("Unstructured").sublist("Read Mesh File").get<std::string>("File");
          std::string suffix(file.substr(file.size()-4, 4));

          if (suffix != ".exo" ) {
            Exceptions::amanzi_throw(Errors::Message("Exodus II was specified as a mesh file format but the suffix of the file that was specified is not .exo"));
          }

          // figure out if this is a parallel run
          std::string framework("Unspecified");
          if (plist->sublist("Mesh").sublist("Unstructured").isSublist("Expert")) {
            framework = plist->sublist("Mesh").sublist("Unstructured").sublist("Expert").get<std::string>("Framework");
          }          

          // Assume that if the framework is unspecified then stk::mesh is used
          // This is obviously a kludge but I don't know how to get around it
          if (numproc_ > 1 && (framework == "Unspecified" || framework == "stk::mesh")) {
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


/* ******************************************************************
 * Empty
 ****************************************************************** */
Teuchos::ParameterList get_Domain_List(Teuchos::ParameterList* plist) {
  Teuchos::ParameterList dom_list;

  if ( plist->isSublist("Domain") ) {
    dom_list = plist->sublist("Domain");
  }

  return dom_list;
}


/* ******************************************************************
 * Empty
 ****************************************************************** */
Teuchos::ParameterList create_MPC_List(Teuchos::ParameterList* plist) {
  Teuchos::ParameterList mpc_list;

  if ( plist->isSublist("Execution Control") ) {
    Teuchos::ParameterList exe_sublist = plist->sublist("Execution Control");

    mpc_list.sublist("Time Integration Mode") = exe_sublist.sublist("Time Integration Mode");

    if (exe_sublist.isSublist("Time Period Control")) {
      mpc_list.sublist("Time Period Control") = exe_sublist.sublist("Time Period Control");
    }

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

    // detect whether transport subcycling is on
    if (plist->sublist("Execution Control").isSublist("Numerical Control Parameters")) {
      if (plist->sublist("Execution Control").sublist("Numerical Control Parameters").isSublist("Unstructured Algorithm")) {
        Teuchos::ParameterList& ncp_list = plist->sublist("Execution Control").sublist("Numerical Control Parameters").sublist("Unstructured Algorithm");
        if (ncp_list.isParameter("transport subcycling")) {
          mpc_list.set<bool>("transport subcycling",ncp_list.get<bool>("transport subcycling"));
        } else {
          mpc_list.set<bool>("transport subcycling", false);
        }
	if (ncp_list.isParameter("max chemistry to transport timestep ratio")) {
	  mpc_list.set<double>("max chemistry to transport timestep ratio",ncp_list.get<double>("max chemistry to transport timestep ratio"));
	} else {
	  mpc_list.set<double>("max chemistry to transport timestep ratio",CHEM_TRANS_DT_RATIO);
	}
      } else {
        mpc_list.set<bool>("transport subcycling", false);
	mpc_list.set<double>("max chemistry to transport timestep ratio",CHEM_TRANS_DT_RATIO);
      }
    } else {
      mpc_list.set<bool>("transport subcycling", false);
      mpc_list.set<double>("max chemistry to transport timestep ratio",CHEM_TRANS_DT_RATIO);
    }


    if ( exe_sublist.isParameter("Flow Model") ) {
      if ( exe_sublist.get<std::string>("Flow Model") == "Off" ) {
        mpc_list.set<std::string>("disable Flow_PK", "yes");
      } else if ( exe_sublist.get<std::string>("Flow Model") == "Richards" ) {
        mpc_list.set<std::string>("disable Flow_PK", "no");
        mpc_list.set<std::string>("Flow model","Richards");
      } else if (exe_sublist.get<std::string>("Flow Model") == "Steady State Richards") {
        mpc_list.set<std::string>("disable Flow_PK", "no");
        mpc_list.set<std::string>("Flow model","Steady State Richards");
      } else if (exe_sublist.get<std::string>("Flow Model") == "Steady State Saturated") {
	mpc_list.set<std::string>("disable Flow_PK", "no");
	mpc_list.set<std::string>("Flow model","Steady State Saturated");      
      } else {
        Exceptions::amanzi_throw(Errors::Message("Flow Model must either be Richards, Steady State Richards, Steady State Saturated, or Off"));
      }
    } else {
      Exceptions::amanzi_throw(Errors::Message("The parameter Flow Model must be specified."));
    }

    if ( exe_sublist.isParameter("Chemistry Model") ) {
      if ( exe_sublist.get<std::string>("Chemistry Model") == "Off" ) {
        mpc_list.set<std::string>("Chemistry Model","Off");
      } else {
        std::string chem_model = exe_sublist.get<std::string>("Chemistry Model");
        mpc_list.set<std::string>("Chemistry Model", chem_model);
      }
    } else {
      Exceptions::amanzi_throw(Errors::Message("The parameter \'Chemistry Model\' must be specified."));
    }


    if ( plist->sublist("Execution Control").isSublist("Restart from Checkpoint Data File") ) {
      mpc_list.sublist("Restart from Checkpoint Data File") =
          plist->sublist("Execution Control").sublist("Restart from Checkpoint Data File");
    }
  }

  mpc_list.sublist("VerboseObject") = create_Verbosity_List(verbosity_level);

  return mpc_list;
}


/* ******************************************************************
 * Populate Transport parameters.
 ****************************************************************** */
Teuchos::ParameterList create_Transport_List(Teuchos::ParameterList* plist) {
  Teuchos::ParameterList trp_list;

  if ( plist->isSublist("Execution Control") ) {
    if ( plist->sublist("Execution Control").isParameter("Transport Model") ) {
      if ( plist->sublist("Execution Control").get<std::string>("Transport Model") == "On" ) {
        if (plist->sublist("Execution Control").isSublist("Numerical Control Parameters")) {
          if (plist->sublist("Execution Control").sublist("Numerical Control Parameters").isSublist("Unstructured Algorithm")) {
            Teuchos::ParameterList& ncp_list = plist->sublist("Execution Control").sublist("Numerical Control Parameters").sublist("Unstructured Algorithm");
            if (ncp_list.isParameter("Transport Integration Algorithm")) {
              std::string tia = ncp_list.get<std::string>("Transport Integration Algorithm");
              if ( tia == "Explicit First-Order" ) {
                trp_list.set<int>("discretization order",1);
              } else if ( tia == "Explicit Second-Order" ) {
                trp_list.set<int>("discretization order",2);
              }
            } else {
              trp_list.set<int>("discretization order",1);
            }
          } else {
            trp_list.set<int>("discretization order",1);
          }
        } else {
          trp_list.set<int>("discretization order",1);
        }
      }

      // continue to set some reasonable defaults
      trp_list.sublist("VerboseObject") = create_Verbosity_List(verbosity_level);
      trp_list.set<std::string>("enable internal tests", "no");
      trp_list.set<double>("CFL", 1.0);
      trp_list.set<std::string>("flow mode", "transient");
      trp_list.set<std::string>("advection limiter", "Tensorial");
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

    // TODO: these simple checks for one transported phase will not
    // work with the addition of the solid phase

    //if ( (++ phase_list.begin()) == phase_list.end() ) {
    if (true) {
      Teuchos::ParameterList& bc_sublist = plist->sublist("Boundary Conditions");

      for (Teuchos::ParameterList::ConstIterator i = bc_sublist.begin(); i != bc_sublist.end(); i++) {
        // read the assigned regions
        Teuchos::Array<std::string> regs = bc_sublist.sublist(bc_sublist.name(i)).get<Teuchos::Array<std::string> >("Assigned Regions");

        // only count sublists
	std::string bc_root_str(bc_sublist.name(i));
        if (bc_sublist.isSublist(bc_sublist.name(i))) {
          if ( bc_sublist.sublist((bc_sublist.name(i))).isSublist("Solute BC")) {
            // read the solute bc stuff
            Teuchos::ParameterList& solbc = bc_sublist.sublist((bc_sublist.name(i))).sublist("Solute BC");

            Teuchos::ParameterList& comps = bc_sublist.sublist((bc_sublist.name(i))).sublist("Solute BC").sublist(phase_name).sublist(phase_comp_name);

            for (Teuchos::Array<std::string>::const_iterator i = comp_names.begin();
                 i != comp_names.end(); i++) {
              if (  comps.isSublist(*i) ) {
                std::stringstream compss;
		compss << *i;
                // for now just read the first value from the
                if ( comps.sublist(*i).isSublist("BC: Uniform Concentration") ) {
		  // create the unique name for this boundary condition
                  std::stringstream ss;
		  ss << bc_root_str << " " << *i;
		  // and create the native boundary condition list
                  Teuchos::ParameterList& bc = tbc_list.sublist(ss.str());

                  Teuchos::ParameterList& bcsub = comps.sublist(*i).sublist("BC: Uniform Concentration");

                  Teuchos::Array<double> values = bcsub.get<Teuchos::Array<double> >("Values");
                  Teuchos::Array<double> times = bcsub.get<Teuchos::Array<double> >("Times");
                  Teuchos::Array<std::string> time_fns = bcsub.get<Teuchos::Array<std::string> >("Time Functions");
                  bc.set<Teuchos::Array<double> >(compss.str(), values );
                  bc.set<Teuchos::Array<double> >("Times", times);
                  bc.set<Teuchos::Array<std::string> >("Time Functions", time_fns);
                  bc.set<Teuchos::Array<std::string> >("Regions", regs);
		}
              }
            }
          }
        }
      }
    } else {
      Exceptions::amanzi_throw(Errors::Message( "Unstructured Amanzi can only have one phase, but the input file specifies more than one."));
    }
  }

  return trp_list;
}


/* ******************************************************************
 * Collects default preconditioners
 ****************************************************************** */
Teuchos::ParameterList create_Preconditioners_List(Teuchos::ParameterList* plist) {
  Teuchos::ParameterList prec_list;
  prec_list.sublist("Trilinos ML") = create_DPC_List(plist);
  prec_list.sublist("Hypre AMG") = create_HypreAMG_List(plist);
  prec_list.sublist("Block ILU") = create_BILU_List(plist);
  return prec_list;
}


/* ******************************************************************
 * Collects linear solvers
 ****************************************************************** */
Teuchos::ParameterList create_Solvers_List(Teuchos::ParameterList* plist) {
  Teuchos::ParameterList solver_list;
  Teuchos::ParameterList& aztecoo_list = solver_list.sublist("AztecOO");

  // define defaults...
  double tol = 1e-14;
  int maxiter = 400;
  std::string method = "GMRES";
  // get values from Execution control list if they exist
  if ( plist->isSublist("Execution Control") ) {  
    if (plist->sublist("Execution Control").isSublist("Numerical Control Parameters")) {
      if (plist->sublist("Execution Control").sublist("Numerical Control Parameters").isSublist("Unstructured Algorithm")) {
	Teuchos::ParameterList& num_list = plist->sublist("Execution Control").sublist("Numerical Control Parameters").sublist("Unstructured Algorithm");
	if (num_list.isParameter("linear solver tolerance"))
	  tol = num_list.get<double>("linear solver tolerance");
	if (num_list.isParameter("linear solver maximum iterations")) 
	  maxiter = num_list.get<int>("linear solver maximum iterations");
	if (num_list.isParameter("linear solver method"))
	  method = num_list.get<std::string>("linear solver method");
      }
    }
  }
  aztecoo_list.set<double>("error tolerance", tol);
  aztecoo_list.set<std::string>("iterative method", method);
  aztecoo_list.set<int>("maximum number of iterations", maxiter);
  return solver_list;
}


/* ******************************************************************
 * Empty
 ****************************************************************** */
Teuchos::ParameterList create_Flow_List(Teuchos::ParameterList* plist) {
  Teuchos::ParameterList flw_list;

  bool use_hypre(false);
  bool use_block_ilu(false);
  if ( plist->isSublist("Execution Control") ) 
    if (plist->sublist("Execution Control").isSublist("Numerical Control Parameters")) 
      if (plist->sublist("Execution Control").sublist("Numerical Control Parameters").isSublist("Unstructured Algorithm")) {    
	Teuchos::ParameterList& num_list = plist->sublist("Execution Control").sublist("Numerical Control Parameters").sublist("Unstructured Algorithm");
	use_hypre =  num_list.get<bool>("use Hypre AMG", false);
	use_block_ilu = num_list.get<bool>("use Block ILU", false);
	if (use_hypre && use_block_ilu) {
	  Exceptions::amanzi_throw(Errors::Message("You can only specify either use HYPRE AMG, or use Block ILU, but not both."));	  
	}
      }
  
  if ( plist->isSublist("Execution Control") ) {
    if ( plist->sublist("Execution Control").isParameter("Flow Model") ) {
      std::string flow_model = plist->sublist("Execution Control").get<std::string>("Flow Model");
      if (flow_model == "Steady State Saturated") {
	Teuchos::ParameterList& darcy_problem = flw_list.sublist("Darcy Problem"); 
	darcy_problem.sublist("VerboseObject") = create_Verbosity_List(verbosity_level); 
	darcy_problem.set<double>("atmospheric pressure", ATMOSPHERIC_PRESSURE); 
	Teuchos::ParameterList& steady_time_integrator = darcy_problem.sublist("steady state time integrator"); 
	steady_time_integrator.set<std::string>("linear solver","AztecOO");
	if (!use_hypre) {
	  steady_time_integrator.set<std::string>("preconditioner", "Trilinos ML");
	} else {
	  steady_time_integrator.set<std::string>("preconditioner", "Hypre AMG");
	}
	// insert the flow BC sublist 
	Teuchos::ParameterList& flow_bc = darcy_problem.sublist("boundary conditions"); 
	flow_bc = create_SS_FlowBC_List(plist); 
      } else if (flow_model == "Richards" || flow_model == "Steady State Richards") {
	Teuchos::ParameterList& richards_problem = flw_list.sublist("Richards Problem"); 
        richards_problem.set<std::string>("relative permeability", "upwind with Darcy flux");
        // this one should come from the input file...
        richards_problem.sublist("VerboseObject") = create_Verbosity_List(verbosity_level);
        richards_problem.set<double>("atmospheric pressure", ATMOSPHERIC_PRESSURE);
	// see if we need to generate a Picard list
	bool use_picard(false);
	Teuchos::ParameterList& ti_mode_list = plist->sublist("Execution Control").sublist("Time Integration Mode");
	if (ti_mode_list.isSublist("Steady")) {
	  use_picard = ti_mode_list.sublist("Steady").get<bool>("Use Picard",false);
	} else if (ti_mode_list.isSublist("Initialize To Steady")) {
	  use_picard = ti_mode_list.sublist("Initialize To Steady").get<bool>("Use Picard",false);
	}
	if (use_picard) {
	  bool have_num_params_list(false);
	  Teuchos::ParameterList num_params_list;
	  if (plist->sublist("Execution Control").isSublist("Numerical Control Parameters"))
	    if (plist->sublist("Execution Control").sublist("Numerical Control Parameters").isSublist("Unstructured Algorithm")) {
	      have_num_params_list = true;
	      num_params_list = plist->sublist("Execution Control").sublist("Numerical Control Parameters").sublist("Unstructured Algorithm");
	    }

	  Teuchos::ParameterList& picard_list = richards_problem.sublist("initial guess pseudo time integrator");
	  
	  if (have_num_params_list) {
	    picard_list.set<bool>("initialize with darcy",num_params_list.get<bool>("pseudo time integrator initialize with darcy",true));
	    picard_list.set<double>("clipping saturation value",num_params_list.get<double>("pseudo time integrator clipping saturation value",0.9));
	    picard_list.set<std::string>("time integration method",num_params_list.get<std::string>("pseudo time integrator time integration method","Picard"));
	    if (!use_hypre) {
	      picard_list.set<std::string>("preconditioner",num_params_list.get<std::string>("pseudo time integrator preconditioner","Trilinos ML"));
	    } else {
	      picard_list.set<std::string>("preconditioner",num_params_list.get<std::string>("pseudo time integrator preconditioner","Hypre AMG"));
	    }
	    picard_list.set<std::string>("linear solver",num_params_list.get<std::string>("pseudo time integrator linear solver","AztecOO"));
	    Teuchos::Array<std::string> error_ctrl(1);
	    error_ctrl[0] = std::string("pressure");
	    picard_list.set<Teuchos::Array<std::string> >("error control options",num_params_list.get<Teuchos::Array<std::string> >("pseudo time integrator error control options",error_ctrl));
	    picard_list.sublist("Picard").set<double>("convergence tolerance",num_params_list.get<double>("pseudo time integrator picard convergence tolerance",PICARD_TOLERANCE));
	    picard_list.sublist("Picard").set<int>("maximum number of iterations",num_params_list.get<int>("pseudo time integrator picard maximum number of iterations",400));
	  } else {
	    picard_list.set<bool>("initialize with darcy",true);
	    picard_list.set<double>("clipping saturation value",0.9);
	    picard_list.set<std::string>("time integration method","Picard");
	    if (!use_hypre) {
	      picard_list.set<std::string>("preconditioner","Trilinos ML");
	    } else {
	      picard_list.set<std::string>("preconditioner","Hypre AMG");
	    }
	    picard_list.set<std::string>("linear solver","AztecOO");
	    Teuchos::Array<std::string> error_ctrl(1);
	    error_ctrl[0] = std::string("pressure");
	    picard_list.set<Teuchos::Array<std::string> >("error control options",error_ctrl);
	    picard_list.sublist("Picard").set<double>("convergence tolerance",PICARD_TOLERANCE);
	    picard_list.sublist("Picard").set<int>("maximum number of iterations",400);
	  }
	}
	  

	bool have_unstructured_algorithm_sublist(false);
	// create sublists for the steady state time integrator
	Teuchos::ParameterList& steady_time_integrator = richards_problem.sublist("steady state time integrator");
	steady_time_integrator.set<std::string>("time integration method","BDF1");
	Teuchos::ParameterList& sti_bdf1 = steady_time_integrator.sublist("BDF1");
	Teuchos::ParameterList& sti_bdf1_param = sti_bdf1.sublist("BDF1 parameters");
	
	// link to preconditioner and linear solver for the steady state time integration
	if (!use_hypre && !use_block_ilu) {
	  steady_time_integrator.set<std::string>("preconditioner", "Trilinos ML");
	} else if (use_hypre) {
	  steady_time_integrator.set<std::string>("preconditioner", "Hypre AMG");
	} else if (use_block_ilu) {
	  steady_time_integrator.set<std::string>("preconditioner", "Block ILU");
	}	  
	steady_time_integrator.set<std::string>("linear solver", "AztecOO");
	
	if (plist->sublist("Execution Control").isSublist("Numerical Control Parameters")) {
	  if (plist->sublist("Execution Control").sublist("Numerical Control Parameters").isSublist("Unstructured Algorithm")) {
	    have_unstructured_algorithm_sublist = true;
	    Teuchos::ParameterList& num_list = plist->sublist("Execution Control").sublist("Numerical Control Parameters").sublist("Unstructured Algorithm");
	    sti_bdf1_param.set<int>("max iterations", num_list.get<int>("steady max iterations",10));
	    sti_bdf1_param.set<int>("min iterations", num_list.get<int>("steady min iterations",5));
	    sti_bdf1_param.set<int>("limit iterations", num_list.get<int>("steady limit iterations",20));
	    sti_bdf1_param.set<double>("nonlinear tolerance", num_list.get<double>("steady nonlinear tolerance",STEADY_NONLINEAR_TOLERANCE));
	    sti_bdf1_param.set<double>("time step reduction factor", num_list.get<double>("steady time step reduction factor",0.8));
	    sti_bdf1_param.set<double>("time step increase factor", num_list.get<double>("steady time step increase factor",1.2));
	    sti_bdf1_param.set<double>("max time step", num_list.get<double>("steady max time step",1.0e+8));
	    sti_bdf1_param.set<int>("max preconditioner lag iterations", num_list.get<int>("steady max preconditioner lag iterations",5));
	    sti_bdf1_param.set<double>("error abs tol", num_list.get<double>("steady error abs tol",1.0));
	    sti_bdf1_param.set<double>("error rel tol", num_list.get<double>("steady error rel tol",0.0));
	    sti_bdf1_param.set<int>("max divergent iterations", num_list.get<int>("steady max divergent iterations",MAX_DIVERGENT_ITERATIONS));
	  }
	}
	if (have_unstructured_algorithm_sublist == false) {
	  // set some probably not so good defaults for the steady computation
	  sti_bdf1_param.set<int>("max iterations",10);
	  sti_bdf1_param.set<int>("min iterations",5);
	  sti_bdf1_param.set<int>("limit iterations",20);
	  sti_bdf1_param.set<double>("nonlinear tolerance",STEADY_NONLINEAR_TOLERANCE);
	  sti_bdf1_param.set<double>("time step reduction factor",0.8);
	  sti_bdf1_param.set<double>("time step increase factor",1.2);
	  sti_bdf1_param.set<double>("max time step", 1.0e+8);
	  sti_bdf1_param.set<int>("max preconditioner lag iterations", 5);
	  sti_bdf1_param.set<double>("error abs tol", 1.0);
	  sti_bdf1_param.set<double>("error rel tol", 0.0);
	  sti_bdf1_param.set<int>("max divergent iterations",MAX_DIVERGENT_ITERATIONS);
	}
	// crerate sublists for the transient time integrator
        Teuchos::ParameterList& transient_time_integrator = richards_problem.sublist("transient time integrator");
        transient_time_integrator.set<std::string>("time integration method", "BDF1");
        Teuchos::ParameterList& tti_bdf1 = transient_time_integrator.sublist("BDF1");
        Teuchos::ParameterList& tti_bdf1_param = tti_bdf1.sublist("BDF1 parameters");

        // link to preconditioner and linear solver for the transient time integrator
	if (!use_hypre && !use_block_ilu) {
	  transient_time_integrator.set<std::string>("preconditioner", "Trilinos ML");
	} else if (use_hypre) {
	  transient_time_integrator.set<std::string>("preconditioner", "Hypre AMG");
	} else if (use_block_ilu) {
	  transient_time_integrator.set<std::string>("preconditioner", "Block ILU");
	}
        transient_time_integrator.set<std::string>("linear solver", "AztecOO");

        have_unstructured_algorithm_sublist = false;
        if (plist->sublist("Execution Control").isSublist("Numerical Control Parameters")) {
          if (plist->sublist("Execution Control").sublist("Numerical Control Parameters").isSublist("Unstructured Algorithm")) {
            have_unstructured_algorithm_sublist = true;
            Teuchos::ParameterList& num_list = plist->sublist("Execution Control").sublist("Numerical Control Parameters").sublist("Unstructured Algorithm");
            tti_bdf1_param.set<int>("max iterations", num_list.get<int>("transient max iterations",10));
            tti_bdf1_param.set<int>("min iterations", num_list.get<int>("transient min iterations",5));
            tti_bdf1_param.set<int>("limit iterations", num_list.get<int>("transient limit iterations",20));
            tti_bdf1_param.set<double>("nonlinear tolerance", num_list.get<double>("transient nonlinear tolerance",TRANSIENT_NONLINEAR_TOLERANCE));
            tti_bdf1_param.set<double>("time step reduction factor", num_list.get<double>("transient time step reduction factor",0.8));
            tti_bdf1_param.set<double>("time step increase factor", num_list.get<double>("transient time step increase factor",1.2));
            tti_bdf1_param.set<double>("max time step", num_list.get<double>("transient max time step",1.0e+8));
            tti_bdf1_param.set<int>("max preconditioner lag iterations", num_list.get<int>("transient max preconditioner lag iterations",5));
            tti_bdf1_param.set<double>("error abs tol", num_list.get<double>("transient error abs tol",1.0));
            tti_bdf1_param.set<double>("error rel tol", num_list.get<double>("transient error rel tol",0.0));
	    tti_bdf1_param.set<int>("max divergent iterations", num_list.get<int>("transient max divergent iterations",MAX_DIVERGENT_ITERATIONS));     
	  }
        }
        if (have_unstructured_algorithm_sublist == false) {
          // set some probably not so good defaults for the steady computation
          tti_bdf1_param.set<int>("max iterations",10);
          tti_bdf1_param.set<int>("min iterations",5);
          tti_bdf1_param.set<int>("limit iterations",20);
          tti_bdf1_param.set<double>("nonlinear tolerance",TRANSIENT_NONLINEAR_TOLERANCE);
          tti_bdf1_param.set<double>("time step reduction factor",0.8);
          tti_bdf1_param.set<double>("time step increase factor",1.2);
          tti_bdf1_param.set<double>("max time step", 1.0e+8);
          tti_bdf1_param.set<int>("max preconditioner lag iterations", 5);
          tti_bdf1_param.set<double>("error abs tol", 1.0);
          tti_bdf1_param.set<double>("error rel tol", 0.0);
	  tti_bdf1_param.set<int>("max divergent iterations",MAX_DIVERGENT_ITERATIONS);
        }


        transient_time_integrator.sublist("VerboseObject") = create_Verbosity_List(verbosity_level);

        // insert the water retention models sublist
        Teuchos::ParameterList &water_retention_models = richards_problem.sublist("Water retention models");
        water_retention_models = create_WRM_List(plist);

        // insert the flow BC sublist
        Teuchos::ParameterList& flow_bc = richards_problem.sublist("boundary conditions");
        flow_bc = create_SS_FlowBC_List(plist);

      } else {
        // something's wrong
      }
    }
  }

  flw_list.sublist("VerboseObject") = create_Verbosity_List(verbosity_level);

  return flw_list;
}


/* ******************************************************************
 * WRM sublist
 ****************************************************************** */
Teuchos::ParameterList create_WRM_List(Teuchos::ParameterList* plist)
{
  Teuchos::ParameterList wrm_list;

  // loop through the material properties list and extract the water retention model info

  Teuchos::ParameterList& matprop_list = plist->sublist("Material Properties");

  int counter = 0;
  for (Teuchos::ParameterList::ConstIterator i = matprop_list.begin(); i != matprop_list.end(); i++) {
    // get the wrm parameters
    Teuchos::ParameterList& cp_list =  matprop_list.sublist(i->first);
    // we can have either van Genuchten or Brooks Corey
    if (cp_list.isSublist("Capillary Pressure: van Genuchten")) {
      Teuchos::ParameterList vG_list = cp_list.sublist("Capillary Pressure: van Genuchten");
      std::string rel_perm = vG_list.get<std::string>("Relative Permeability");
      if (rel_perm != "Mualem" && rel_perm != "Burdine") {
        std::stringstream ss;
        ss << "Currently we only have Mualem or Burdine as the relative permeability models";
        Exceptions::amanzi_throw(Errors::Message(ss.str().c_str()));
      }

      double alpha = vG_list.get<double>("alpha");
      double Sr = vG_list.get<double>("Sr");
      double m = vG_list.get<double>("m");
      double ell = vG_list.get<double>("ell",0.5);
      double krel_smooth = vG_list.get<double>("krel smoothing interval", 0.0);
      if (krel_smooth < 0.0) {
        Exceptions::amanzi_throw(Errors::Message("If krel smoothing interval is specified it must be positive."));
      }


      // now get the assigned regions
      Teuchos::Array<std::string> regions = cp_list.get<Teuchos::Array<std::string> >("Assigned Regions");

      for (Teuchos::Array<std::string>::const_iterator i = regions.begin();
           i != regions.end(); i++) {
        std::stringstream ss;
        ss << "Water Retention Model for " << *i;

        Teuchos::ParameterList& wrm_sublist = wrm_list.sublist(ss.str());

        wrm_sublist.set<std::string>("Water retention model", "van Genuchten");
        wrm_sublist.set<std::string>("Region",*i);
        wrm_sublist.set<double>("van Genuchten m", m);
        wrm_sublist.set<double>("van Genuchten l", ell);
        wrm_sublist.set<double>("van Genuchten alpha",alpha);
        wrm_sublist.set<double>("residual saturation", Sr);
        wrm_sublist.set<double>("regularization interval", krel_smooth);
        wrm_sublist.set<std::string>("relative permeability model", rel_perm);
      }
    } else if (cp_list.isSublist("Capillary Pressure: Brooks Corey")) {
      Teuchos::ParameterList& BC_list = cp_list.sublist("Capillary Pressure: Brooks Corey");

      std::string rel_perm = BC_list.get<std::string>("Relative Permeability");
      if (rel_perm != "Mualem" && rel_perm != "Burdine") {
        std::stringstream ss;
        ss << "Currently we only have Mualem or Burdine as the relative permeability models";
        Exceptions::amanzi_throw(Errors::Message(ss.str().c_str()));
      }

      double lambda      = BC_list.get<double>("lambda");
      double alpha       = BC_list.get<double>("alpha");
      double ell         = BC_list.get<double>("ell",0.5);
      double Sr          = BC_list.get<double>("Sr");
      double krel_smooth = BC_list.get<double>("krel smoothing interval",0.0);

      if (krel_smooth < 0.0) {
        Exceptions::amanzi_throw(Errors::Message("If krel smoothing interval is specified it must be positive."));
      }

      // now get the assigned regions
      Teuchos::Array<std::string> regions = cp_list.get<Teuchos::Array<std::string> >("Assigned Regions");

      for (Teuchos::Array<std::string>::const_iterator i = regions.begin();
           i != regions.end(); i++) {
        std::stringstream ss;
        ss << "Water Retention Model for " << *i;

        Teuchos::ParameterList& wrm_sublist = wrm_list.sublist(ss.str());

        wrm_sublist.set<std::string>("Water retention model", "Brooks Corey");
        wrm_sublist.set<std::string>("Region",*i);
        wrm_sublist.set<double>("Brooks Corey lambda", lambda);
        wrm_sublist.set<double>("Brooks Corey alpha",alpha);
        wrm_sublist.set<double>("Brooks Corey l",ell);
        wrm_sublist.set<double>("residual saturation", Sr);
        wrm_sublist.set<double>("regularization interval", krel_smooth);
        wrm_sublist.set<std::string>("relative permeability model", rel_perm);
      }
    } else {
      // not implemented error
      Exceptions::amanzi_throw(Errors::Message("An unknown capillary pressure model was specified, must specify either van Genuchten or Brooks Corey"));
    }
  }
  return wrm_list;
}


/* ******************************************************************
 * ML preconditioner sublist
 ****************************************************************** */
Teuchos::ParameterList create_DPC_List(Teuchos::ParameterList* plist)
{
  Teuchos::ParameterList dpc_list;

  double aggthr(0.0);
  std::string smthtyp("Jacobi");
  int ncycles(2);
  int nsmooth(3);

  if (plist->sublist("Execution Control").isSublist("Numerical Control Parameters")) {
    if (plist->sublist("Execution Control").sublist("Numerical Control Parameters").isSublist("Unstructured Algorithm")) {
      Teuchos::ParameterList& ncp_list = plist->sublist("Execution Control").sublist("Numerical Control Parameters").sublist("Unstructured Algorithm");
      if (ncp_list.isParameter("ML aggregation threshold")) {
        aggthr = ncp_list.get<double>("ML aggregation threshold");
      }
      if (ncp_list.isParameter("ML smoother type")) {
        smthtyp = ncp_list.get<std::string>("ML smoother type");
      }
      if (ncp_list.isParameter("ML cycle applications")) {
        ncycles = ncp_list.get<int>("ML cycle applications");
      }
      if (ncp_list.isParameter("ML smoother sweeps")) {
        nsmooth = ncp_list.get<int>("ML smoother sweeps");
      }
    }
  }

  Teuchos::ParameterList& ml_list = dpc_list.sublist("ML Parameters");
  ml_list.set<int>("ML output", 0);
  ml_list.set<int>("max levels", 40);
  ml_list.set<std::string>("prec type","MGV");
  ml_list.set<int>("cycle applications", ncycles);
  ml_list.set<std::string>("aggregation: type", "Uncoupled-MIS");
  ml_list.set<double>("aggregation: damping factor", 1.33333);
  ml_list.set<double>("aggregation: threshold", aggthr);
  ml_list.set<std::string>("eigen-analysis: type","cg");
  ml_list.set<int>("eigen-analysis: iterations", 10);
  ml_list.set<int>("smoother: sweeps", nsmooth);
  ml_list.set<double>("smoother: damping factor", 1.0);
  ml_list.set<std::string>("smoother: pre or post", "both");
  ml_list.set<std::string>("smoother: type", smthtyp);
  ml_list.set<double>("smoother: damping factor", 1.0);
  ml_list.set<std::string>("coarse: type", "Amesos-KLU");
  ml_list.set<int>("coarse: max size", 256);

  return dpc_list;
}

/* ******************************************************************
 * Block ILU preconditioner sublist
 ****************************************************************** */
Teuchos::ParameterList create_BILU_List(Teuchos::ParameterList* plist)
{
  Teuchos::ParameterList bilu_list;
  
  double bilu_relax_value(1.0);
  double bilu_abs_thresh(0.0);
  double bilu_rel_thresh(1.0);
  int bilu_level_of_fill(0);
  int bilu_overlap(0);

  if (plist->sublist("Execution Control").isSublist("Numerical Control Parameters")) {
    if (plist->sublist("Execution Control").sublist("Numerical Control Parameters").isSublist("Unstructured Algorithm")) {
      Teuchos::ParameterList& ncp_list = plist->sublist("Execution Control").sublist("Numerical Control Parameters").sublist("Unstructured Algorithm");
      if (ncp_list.isParameter("Block ILU relax value")) {
        bilu_relax_value = ncp_list.get<double>("Block ILU relax value");
      }
      if (ncp_list.isParameter("Block ILU relative threshold")) {
        bilu_rel_thresh = ncp_list.get<double>("Block ILU relative threshold");
      }      
      if (ncp_list.isParameter("Block ILU absolute threshold")) {
        bilu_abs_thresh = ncp_list.get<double>("Block ILU absolute threshold");
      }
      if (ncp_list.isParameter("Block ILU level of fill")) {
        bilu_level_of_fill = ncp_list.get<int>("Block ILU level of fill");
      }
      if (ncp_list.isParameter("Block ILU overlap")) {
        bilu_overlap = ncp_list.get<int>("Block ILU overlap");
      }      
    }
  }
      
  Teuchos::ParameterList& p_list = bilu_list.sublist("Block ILU Parameters");
  p_list.set<double>("fact: relax value",bilu_relax_value);
  p_list.set<double>("fact: absolute threshold",bilu_abs_thresh);
  p_list.set<double>("fact: relative threshold",bilu_rel_thresh);
  p_list.set<int>("fact: level-of-fill",bilu_level_of_fill);
  p_list.set<int>("overlap",bilu_overlap);
  p_list.set<std::string>("schwarz: combine mode","Add");

  return bilu_list;
}



/* ******************************************************************
 * Hypre BoomerAMG preconditioner sublist
 ****************************************************************** */
Teuchos::ParameterList create_HypreAMG_List(Teuchos::ParameterList* plist)
{
  Teuchos::ParameterList dpc_list;

  dpc_list.set<std::string>("discretization method", "optimized mfd");

  double tol(0.0);
  int ncycles(5);
  int nsmooth(3);
  double strong_threshold(0.5); 

  if (plist->sublist("Execution Control").isSublist("Numerical Control Parameters")) {
    if (plist->sublist("Execution Control").sublist("Numerical Control Parameters").isSublist("Unstructured Algorithm")) {
      Teuchos::ParameterList& ncp_list = plist->sublist("Execution Control").sublist("Numerical Control Parameters").sublist("Unstructured Algorithm");
      if (ncp_list.isParameter("Hypre AMG cycle applications")) {
        ncycles = ncp_list.get<int>("Hypre AMG cycle applications");
      }
      if (ncp_list.isParameter("Hypre AMG smoother sweeps")) {
        nsmooth = ncp_list.get<int>("Hypre AMG smoother sweeps");
      }
      if (ncp_list.isParameter("Hypre AMG tolerance")) {
	tol = ncp_list.get<double>("Hypre AMG tolerance");
      }
      if (ncp_list.isParameter("Hypre AMG strong threshold")) {
	strong_threshold = ncp_list.get<double>("Hypre AMG strong threshold");
      }
    }
  }      

  Teuchos::ParameterList& amg_list = dpc_list.sublist("BoomerAMG Parameters");
  amg_list.set<double>("tolerance", tol);
  amg_list.set<int>("smoother sweeps", nsmooth);
  amg_list.set<int>("cycle applications", ncycles);
  amg_list.set<double>("strong threshold", strong_threshold);

  return dpc_list;
}


/* ******************************************************************
 * DPC sublist
 ****************************************************************** */
Teuchos::ParameterList create_SS_FlowBC_List(Teuchos::ParameterList* plist) {
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
      } else if ( bc.isSublist("BC: Seepage") ) {
        Teuchos::ParameterList& bc_flux = bc.sublist("BC: Seepage");

        Teuchos::Array<double> times = bc_flux.get<Teuchos::Array<double> >("Times");
        Teuchos::Array<std::string> time_fns = bc_flux.get<Teuchos::Array<std::string> >("Time Functions");

        if (! bc_flux.isParameter("Inward Mass Flux") )  {
          // we can only handle mass fluxes right now
          Exceptions::amanzi_throw(Errors::Message("In BC: Seepage we can only handle Inward Mass Flux"));
        }

        Teuchos::Array<double> flux;

        flux = bc_flux.get<Teuchos::Array<double> >("Inward Mass Flux");
        for (int i=0; i<flux.size(); i++) flux[i] = - flux[i];

        std::stringstream ss;
        ss << "BC " << bc_counter++;

        Teuchos::ParameterList& tbc = ssf_list.sublist("seepage face").sublist(ss.str());
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

      }

      // TODO...
      // add the rest of the boundary conditions
    }
  }

  return ssf_list;
}


/* ******************************************************************
 * populates parameters in the State list.
 ****************************************************************** */
Teuchos::ParameterList create_State_List(Teuchos::ParameterList* plist) {
  Teuchos::ParameterList stt_list;

  stt_list.set<double>("Gravity x", 0.0);
  if (spatial_dimension_ == 2) {
    stt_list.set<double>("Gravity y", - GRAVITY_MAGNITUDE);
  } else {
    stt_list.set<double>("Gravity y", 0.0);
    stt_list.set<double>("Gravity z", - GRAVITY_MAGNITUDE);
  }

  // find the viscosity
  Teuchos::ParameterList& phase_list = plist->sublist("Phase Definitions");

  // make sure there is only one phase
  //if ( (++ phase_list.begin()) == phase_list.end() ) {
  // TODO: these simple checks for one transported phase will not work
  // with the addition of the solid phase
  if (true) {
    // write the array of component solutes
    stt_list.set<Teuchos::Array<std::string> >("Component Solutes", comp_names);
    stt_list.set<int>("Number of component concentrations", comp_names.size());

    // write the array of mineral names. hopefully order is preserved...
    if (mineral_names_.size()) {
      stt_list.set<Teuchos::Array<std::string> >("Minerals", mineral_names_);
    }
    if (sorption_site_names_.size()) {
      stt_list.set<Teuchos::Array<std::string> >("Sorption Sites", sorption_site_names_);
    }

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
      
      // particle density, for now we make the default 1.0 since we do not want to require this input parameter in the input spec, yet
      double particle_density = matprop_list.sublist(matprop_list.name(i)).sublist("Particle Density: Uniform").get<double>("Value",PARTICLE_DENSITY);

      // extract the mineralogy for this material
      Teuchos::ParameterList mineralogy;
      if (matprop_list.sublist(matprop_list.name(i)).isSublist("Mineralogy")) {
        mineralogy = matprop_list.sublist(matprop_list.name(i)).sublist("Mineralogy");
      }
      // extract the isotherms for this material
      Teuchos::ParameterList isotherms;
      if (matprop_list.sublist(matprop_list.name(i)).isSublist("Sorption Isotherms")) {
        isotherms = matprop_list.sublist(matprop_list.name(i)).sublist("Sorption Isotherms");
      }

      // extract the surface complexation sites for this material
      Teuchos::ParameterList surface_sites;
      if (matprop_list.sublist(matprop_list.name(i)).isSublist("Surface Complexation Sites")) {
        surface_sites = matprop_list.sublist(matprop_list.name(i)).sublist("Surface Complexation Sites");
      }

      // extract the ion exchange information, may be a list in the future.
      // real cation exchange capacities are > 0, so we can use cec < 0 as a flag.
      double cec = -1.0;
      if (matprop_list.sublist(matprop_list.name(i)).isParameter("Cation Exchange Capacity")) {
        cec = matprop_list.sublist(matprop_list.name(i)).get<double>("Cation Exchange Capacity");
      }

      for (Teuchos::Array<std::string>::const_iterator i=regions.begin(); i!=regions.end(); i++) {
        std::stringstream sss;
        sss << "Mesh block " << *i;

        Teuchos::ParameterList& stt_mat = stt_list.sublist(sss.str());

        stt_mat.set<double>("Constant porosity", porosity);
        stt_mat.set<double>("Constant vertical permeability", perm_vert);
        stt_mat.set<double>("Constant horizontal permeability", perm_horiz);
	stt_mat.set<double>("Constant particle density", particle_density);
        stt_mat.set<std::string>("Region", *i);

        if (  mineralogy.begin() != mineralogy.end() ) { // this is to avoid creating an empty Mineralogy list
          Teuchos::ParameterList& region_mineralogy = stt_mat.sublist("Mineralogy");
          region_mineralogy = mineralogy;
        }

        if ( isotherms.begin() != isotherms.end() ) { // this is to avoid creating an empty Sorption Isotherms list
          Teuchos::ParameterList& region_isotherms = stt_mat.sublist("Sorption Isotherms");
          region_isotherms = isotherms;
        }

        if ( surface_sites.begin() != surface_sites.end() ) { // this is to avoid creating an empty Surface Complexation Sites list
          Teuchos::ParameterList& region_surface_sites = stt_mat.sublist("Surface Complexation Sites");
          region_surface_sites = surface_sites;
        }

        if (cec > 0.0) {
          stt_mat.set<double>("Cation Exchange Capacity", cec);
        }
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
        } else if (ic_for_region->isSublist("IC: File Pressure")) {
          Teuchos::ParameterList& sublist = stt_mat.sublist("file pressure");
          sublist.set<std::string>("file name", ic_for_region->sublist("IC: File Pressure").get<std::string>("File"));
          sublist.set<std::string>("label", ic_for_region->sublist("IC: File Pressure").get<std::string>("Label"));
        } else {
          Exceptions::amanzi_throw(Errors::Message("An initial condition for pressure must be specified. It must either be IC: Uniform Pressure, IC: Linear Pressure, or IC: File Pressure."));
        }

        // write the initial conditions for saturation, since this is not a primary variable, this is not required
        // so we don't throw if these initial conditions are missing
        if ( ic_for_region->isSublist("IC: Uniform Saturation") ) {
          Teuchos::ParameterList& sublist = stt_mat.sublist("uniform saturation");
          sublist.set<double>("value",ic_for_region->sublist("IC: Uniform Saturation").get<double>("Value"));
        } else if (ic_for_region->isSublist("IC: Linear Saturation")) {
          Teuchos::ParameterList& sublist = stt_mat.sublist("linear saturation");
        }
	// write the initial conditions for velocity
	// we don't throw if these initial conditions are missing
        if ( ic_for_region->isSublist("IC: Uniform Velocity") ) {
	  Teuchos::Array<double> vel = ic_for_region->sublist("IC: Uniform Velocity").get<Teuchos::Array<double> >("Velocity Vector");
	  stt_mat.set<double>("Constant velocity x",vel[0]);
	  if (spatial_dimension_>1) stt_mat.set<double>("Constant velocity y", vel[1]);
	  if (spatial_dimension_>2) stt_mat.set<double>("Constant velocity z", vel[2]);
	  
	}

	if (plist->sublist("Execution Control").get<std::string>("Transport Model") != std::string("Off")  ||
	    plist->sublist("Execution Control").get<std::string>("Chemistry Model") != std::string("Off") ) {
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
	    
	    conc = ic_for_region->sublist("Solute IC").sublist(phase_name).sublist(phase_comp_name).sublist(comp_names[ii]).sublist("IC: Uniform Concentration").get<double>("Free Ion Guess", 1.0e-9);
	    ss.clear();
	    ss.str("");
	    ss << "Free Ion Guess " << comp_names_map[ comp_names[ii] ];
	    stt_mat.set<double>(ss.str(), conc);
	    ss.clear();
	    ss.str("");
	  }
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


/* ******************************************************************
 * Empty
 ****************************************************************** */
Teuchos::ParameterList create_Verbosity_List(const std::string& vlevel) {
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


/* ******************************************************************
 * Empty
 ****************************************************************** */
Teuchos::ParameterList CreateChemistryList(Teuchos::ParameterList* plist) {
  Teuchos::ParameterList chemistry_list;
  if ( plist->isSublist("Chemistry") ) {
    chemistry_list = plist->sublist("Chemistry");
  }
  return chemistry_list;
}  // end CreateChemistryList()


}
}
