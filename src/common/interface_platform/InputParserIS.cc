#include "InputParserIS.hh"
#include "InputParserIS-defaults.hh"
#include "Teuchos_XMLParameterListHelpers.hpp"

#include <sstream>
#include <string>
#include <algorithm>
#include <boost/lambda/lambda.hpp>
#include <boost/bind.hpp>

#include "errors.hh"
#include "exceptions.hh"
#include "dbc.hh"

#define  BOOST_FILESYTEM_NO_DEPRECATED
#include "boost/filesystem/operations.hpp"
#include "boost/filesystem/path.hpp"
#include "boost/format.hpp"
#include "boost/lexical_cast.hpp"

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

  // first make sure the version is correct
  check_AmanziInputVersion(plist);

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
  Teuchos::Array<int> cycle_range(0);

  if ( plist->sublist("Output").sublist("Cycle Macros").isSublist(macro_name) ) {
    Teuchos::ParameterList& clist = plist->sublist("Output").sublist("Cycle Macros").sublist(macro_name);
    if (clist.isParameter("Start_Period_Stop"))
      cycle_range = clist.get<Teuchos::Array<int> >("Start_Period_Stop");
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
Teuchos::Array<int> get_Cycle_Macro_Values(const std::string& macro_name, Teuchos::ParameterList* plist) {
  Teuchos::Array<int> values(0);

  if ( plist->sublist("Output").sublist("Cycle Macros").isSublist(macro_name) ) {
    Teuchos::ParameterList& clist = plist->sublist("Output").sublist("Cycle Macros").sublist(macro_name);
    if (clist.isParameter("Values")) values = clist.get<Teuchos::Array<int> >("Values");
  } else {
    std::stringstream ss;
    ss << "The cycle macro " << macro_name << " does not exist in the input file";
    Exceptions::amanzi_throw(Errors::Message(ss.str().c_str()));
  }

  return values;
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

    std::string verbosity = plist->sublist("Execution Control").get<std::string>("Verbosity",VERBOSITY_DEFAULT);

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
              obs_list.sublist(i->first).remove("Values");
            }
            obs_list.sublist(i->first).remove("Time Macro");
          }
          if (obs_list.sublist(i->first).isParameter("Cycle Macro")) {
            std::string cycle_macro = obs_list.sublist(i->first).get<std::string>("Cycle Macro");

            Teuchos::Array<int> sps = get_Cycle_Macro(cycle_macro, plist);
            Teuchos::Array<int> values = get_Cycle_Macro_Values(cycle_macro, plist);

            if (sps.size() != 3  && values.size() == 0) {
              Errors::Message message("Cycle macro " + cycle_macro + " has neither a valid Start_Period_Stop nor a valid Values parameter");
              Exceptions::amanzi_throw(message);
            }

            if (sps.size() == 3)
              obs_list.sublist(i->first).sublist("cycle start period stop").sublist(cycle_macro).set("start period stop", sps);
            if (values.size() > 0)
              obs_list.sublist(i->first).set("cycles", values);

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
          if (meshfile != file) {
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
          //
          // We also have to be able to tell if we have prepartitioned
          // files or if we have one file that we want to partition

          if (numproc_ > 1) {
            std::string par_file(file);
            par_file.replace(file.size()-4,4,std::string(".par"));

            // attach the right extensions as required by Nemesis file naming conventions
            // in which files are named as mymesh.par.N.r where N = numproc and r is rank

            int rank = numproc_-1;
            int ndigits = (int)floor(log10(numproc_)) + 1;
            std::string fmt = boost::str(boost::format("%%s.%%d.%%0%dd") % ndigits);
            std::string par_file_w_ext = boost::str(boost::format(fmt) % par_file % numproc_ % rank);
            boost::filesystem::path p(par_file_w_ext);

            if (boost::filesystem::exists(p))
              fn_list.set<std::string>("File",par_file); // Nemesis file exists. Use the .par extension
            else
              fn_list.set<std::string>("File",file); // Use original .exo file extension

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

  bool transport_on(false), chemistry_on(false);

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
        transport_on = true;
        mpc_list.set<std::string>("disable Transport_PK","no");

      } else {
        Exceptions::amanzi_throw(Errors::Message("Transport Model must either be On or Off"));
      }
    } else {
      Exceptions::amanzi_throw(Errors::Message("The parameter Transport Model must be specified."));
    }

    if ( exe_sublist.isParameter("Chemistry Model") ) {
      if ( exe_sublist.get<std::string>("Chemistry Model") == "Off" ) {
        mpc_list.set<std::string>("Chemistry Model","Off");
      } else {
        chemistry_on = true;
        std::string chem_model = exe_sublist.get<std::string>("Chemistry Model");
        mpc_list.set<std::string>("Chemistry Model", chem_model);
      }
    } else {
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

    if ( exe_sublist.isParameter("Flow Model") ) {
      if ( exe_sublist.get<std::string>("Flow Model") == "Off") {
        mpc_list.set<std::string>("disable Flow_PK", "yes");
      } else if ( exe_sublist.get<std::string>("Flow Model") == "Richards") {
        mpc_list.set<std::string>("disable Flow_PK", "no");
        mpc_list.set<std::string>("Flow model","Richards");
      } else if (exe_sublist.get<std::string>("Flow Model") == "Single Phase") {
        mpc_list.set<std::string>("disable Flow_PK", "no");
        mpc_list.set<std::string>("Flow model","Steady State Saturated");
      } else {
        Exceptions::amanzi_throw(Errors::Message("Flow Model must either be Richards, Single Phase, or Off"));
      }
    } else {
      Exceptions::amanzi_throw(Errors::Message("The parameter Flow Model must be specified."));
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

        // transport is on, set some defaults
        trp_list.set<int>("discretization order",1);
        trp_list.sublist("VerboseObject") = create_Verbosity_List(verbosity_level);
        trp_list.set<std::string>("enable internal tests", "no");
        trp_list.set<double>("CFL", 1.0);
        trp_list.set<std::string>("flow mode", "transient");
        trp_list.set<std::string>("advection limiter", "Tensorial");

        if (plist->sublist("Execution Control").isSublist("Numerical Control Parameters")) {
          Teuchos::ParameterList& ec_ncp_list = plist->sublist("Execution Control").sublist("Numerical Control Parameters");
          if (ec_ncp_list.isSublist("Unstructured Algorithm")) {
            Teuchos::ParameterList ec_ncp_u_list = ec_ncp_list.sublist("Unstructured Algorithm");
            if (ec_ncp_u_list.isSublist("Transport Process Kernel")) {
              Teuchos::ParameterList& t_list = ec_ncp_u_list.sublist("Transport Process Kernel");
              if (t_list.isParameter("Transport Integration Algorithm")) {
                std::string tia = t_list.get<std::string>("Transport Integration Algorithm");
                if ( tia == "Explicit First-Order" ) {
                  trp_list.set<int>("discretization order",1);
                } else if ( tia == "Explicit Second-Order" ) {
                  trp_list.set<int>("discretization order",2);
                }
              }
            }
          }
        }
      }
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
  double tol = LIN_SOLVE_TOL;
  int maxiter = LIN_SOLVE_MAXITER;
  std::string method = LIN_SOLVE_METHOD;
  std::string prec = LIN_SOLVE_PREC;
  // get values from Execution control list if they exist
  if (plist->sublist("Execution Control").isSublist("Numerical Control Parameters")) {
    Teuchos::ParameterList& ncp_list = plist->sublist("Execution Control").sublist("Numerical Control Parameters");
    if (ncp_list.isSublist("Unstructured Algorithm")) {
      Teuchos::ParameterList& ncpu_list = ncp_list.sublist("Unstructured Algorithm");
      if (ncpu_list.isSublist("Linear Solver")) {
        Teuchos::ParameterList& num_list = ncpu_list.sublist("Linear Solver");
        if (num_list.isParameter("linear solver tolerance"))
          tol = num_list.get<double>("linear solver tolerance");
        if (num_list.isParameter("linear solver maximum iterations"))
          maxiter = num_list.get<int>("linear solver maximum iterations");
        if (num_list.isParameter("linear solver method"))
          method = num_list.get<std::string>("linear solver method");
        if (num_list.isParameter("linear solver method"))
          prec = num_list.get<std::string>("linear solver preconditioner");
      }
    }
  }
  aztecoo_list.set<double>("error tolerance", tol);
  aztecoo_list.set<std::string>("iterative method", method);
  aztecoo_list.set<int>("maximum number of iterations", maxiter);
  aztecoo_list.set<std::string>("preconditioner", prec);

  return solver_list;
}


/* ******************************************************************
 * Empty
 ****************************************************************** */
Teuchos::ParameterList create_Flow_List(Teuchos::ParameterList* plist) {
  Teuchos::ParameterList flw_list;

  if ( plist->isSublist("Execution Control") ) {
    if ( plist->sublist("Execution Control").isParameter("Flow Model") ) {
      std::string flow_model = plist->sublist("Execution Control").get<std::string>("Flow Model");
      if (flow_model == "Single Phase") {
        Teuchos::ParameterList& darcy_problem = flw_list.sublist("Darcy Problem");
        darcy_problem.sublist("VerboseObject") = create_Verbosity_List(verbosity_level);
        darcy_problem.set<double>("atmospheric pressure", ATMOSPHERIC_PRESSURE);
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
        bool use_picard(USE_PICARD);
        Teuchos::ParameterList& ti_mode_list = plist->sublist("Execution Control").sublist("Time Integration Mode");
        if (ti_mode_list.isSublist("Steady")) {
          use_picard = ti_mode_list.sublist("Steady").get<bool>("Use Picard",USE_PICARD);
        } else if (ti_mode_list.isSublist("Initialize To Steady")) {
          use_picard = ti_mode_list.sublist("Initialize To Steady").get<bool>("Use Picard",USE_PICARD);
        }
        if (use_picard) {
          bool have_picard_params_list(false);
          Teuchos::ParameterList picard_params_list;
          if (plist->sublist("Execution Control").isSublist("Numerical Control Parameters")) {
            Teuchos::ParameterList& ncp_list = plist->sublist("Execution Control").sublist("Numerical Control Parameters");
            if (ncp_list.isSublist("Unstructured Algorithm")) {
              Teuchos::ParameterList& ncpu_list = ncp_list.sublist("Unstructured Algorithm");
              if (ncpu_list.isSublist("Steady-State Pseudo-Time Implicit Solver")) {
                have_picard_params_list = true;
                picard_params_list = ncpu_list.sublist("Steady-State Pseudo-Time Implicit Solver");
              }
            }
          }

          Teuchos::ParameterList& picard_list = richards_problem.sublist("initial guess pseudo time integrator");

          if (have_picard_params_list) {
            picard_list.set<bool>("initialize with darcy",picard_params_list.get<bool>("pseudo time integrator initialize with darcy",PIC_INIT_DARCY));
            picard_list.set<double>("clipping saturation value",picard_params_list.get<double>("pseudo time integrator clipping saturation value",PIC_CLIP_SAT));
            picard_list.set<std::string>("time integration method",picard_params_list.get<std::string>("pseudo time integrator time integration method",PIC_METHOD));
            picard_list.set<std::string>("preconditioner",picard_params_list.get<std::string>("pseudo time integrator preconditioner",PIC_PRECOND));
            picard_list.set<std::string>("linear solver",picard_params_list.get<std::string>("pseudo time integrator linear solver",PIC_SOLVE));
            Teuchos::Array<std::string> error_ctrl(1);
            error_ctrl[0] = std::string(PIC_ERROR_METHOD);
            picard_list.set<Teuchos::Array<std::string> >("error control options",picard_params_list.get<Teuchos::Array<std::string> >("pseudo time integrator error control options",error_ctrl));
            picard_list.sublist("Picard").set<double>("convergence tolerance",picard_params_list.get<double>("pseudo time integrator picard convergence tolerance",PICARD_TOLERANCE));
            picard_list.sublist("Picard").set<int>("maximum number of iterations",picard_params_list.get<int>("pseudo time integrator picard maximum number of iterations",PIC_MAX_ITER));
          } else {
            picard_list.set<bool>("initialize with darcy",PIC_INIT_DARCY);
            picard_list.set<double>("clipping saturation value",PIC_CLIP_SAT);
            picard_list.set<std::string>("time integration method",PIC_METHOD);
            picard_list.set<std::string>("preconditioner",PIC_PRECOND);
            picard_list.set<std::string>("linear solver",PIC_SOLVE);
            Teuchos::Array<std::string> error_ctrl(1);
            error_ctrl[0] = std::string(PIC_ERROR_METHOD);
            picard_list.set<Teuchos::Array<std::string> >("error control options",error_ctrl);
            picard_list.sublist("Picard").set<double>("convergence tolerance",PICARD_TOLERANCE);
            picard_list.sublist("Picard").set<int>("maximum number of iterations",PIC_MAX_ITER);
          }
        }


        // create sublists for the steady state time integrator
        Teuchos::ParameterList& steady_time_integrator = richards_problem.sublist("steady state time integrator");
        steady_time_integrator.set<std::string>("time integration method","BDF1");
        Teuchos::ParameterList& sti_bdf1 = steady_time_integrator.sublist("BDF1");
        Teuchos::ParameterList& sti_bdf1_param = sti_bdf1.sublist("BDF1 parameters");

        steady_time_integrator.set<std::string>("preconditioner", ST_PRECOND);
        steady_time_integrator.set<std::string>("linear solver", ST_SOLVER);

        // set defaults
        sti_bdf1_param.set<int>("max iterations",ST_MAX_ITER);
        sti_bdf1_param.set<int>("min iterations",ST_MIN_ITER);
        sti_bdf1_param.set<int>("limit iterations",ST_LIMIT_ITER);
        sti_bdf1_param.set<double>("nonlinear tolerance",STEADY_NONLINEAR_TOLERANCE);
        sti_bdf1_param.set<double>("time step reduction factor",ST_TS_RED_FACTOR);
        sti_bdf1_param.set<double>("time step increase factor",ST_TS_INC_FACTOR);
        sti_bdf1_param.set<double>("max time step", ST_MAX_TS);
        sti_bdf1_param.set<int>("max preconditioner lag iterations", ST_MAX_PREC_LAG);
        sti_bdf1_param.set<double>("error abs tol", ST_ERROR_ABS_TOL);
        sti_bdf1_param.set<double>("error rel tol", ST_ERROR_REL_TOL);
        sti_bdf1_param.set<int>("max divergent iterations",ST_MAX_DIVERGENT_ITERATIONS);
        sti_bdf1_param.set<double>("nonlinear iteration damping factor",ST_NONLIN_DAMP);
        sti_bdf1_param.set<int>("nonlinear iteration initial guess extrapolation order",ST_NONLIN_INIT_GUESS_EXTR_ORD);
        sti_bdf1_param.set<double>("restart tolerance relaxation factor",ST_NONLIN_INIT_TS_FACTOR);
        sti_bdf1_param.set<double>("restart tolerance relaxation factor damping",ST_NONLIN_INIT_TS_FACTOR_DAMP);
        sti_bdf1_param.set<double>("nonlinear iteration divergence factor",ST_DIVERG_FACT);

        if (plist->sublist("Execution Control").isSublist("Numerical Control Parameters")) {
          Teuchos::ParameterList& ncp_list =  plist->sublist("Execution Control").sublist("Numerical Control Parameters");
          if (ncp_list.isSublist("Unstructured Algorithm")) {
            Teuchos::ParameterList& ncpu_list = ncp_list.sublist("Unstructured Algorithm");
            if (ncpu_list.isSublist("Steady-State Implicit Time Integration")) {
              Teuchos::ParameterList& num_list = ncpu_list.sublist("Steady-State Implicit Time Integration");
              sti_bdf1_param.set<int>("max iterations",
                                      num_list.get<int>("steady max iterations",ST_MAX_ITER));
              sti_bdf1_param.set<int>("min iterations",
                                      num_list.get<int>("steady min iterations",ST_MIN_ITER));
              sti_bdf1_param.set<int>("limit iterations",
                                      num_list.get<int>("steady limit iterations",ST_LIMIT_ITER));
              sti_bdf1_param.set<double>("nonlinear tolerance",
                                         num_list.get<double>("steady nonlinear tolerance",STEADY_NONLINEAR_TOLERANCE));
              sti_bdf1_param.set<double>("time step reduction factor",
                                         num_list.get<double>("steady time step reduction factor",ST_TS_RED_FACTOR));
              sti_bdf1_param.set<double>("time step increase factor",
                                         num_list.get<double>("steady time step increase factor",ST_TS_INC_FACTOR));
              sti_bdf1_param.set<double>("max time step", num_list.get<double>("steady max time step",ST_MAX_TS));
              sti_bdf1_param.set<int>("max preconditioner lag iterations",
                                      num_list.get<int>("steady max preconditioner lag iterations",ST_MAX_PREC_LAG));
              sti_bdf1_param.set<double>("error abs tol", num_list.get<double>("steady error abs tol",ST_ERROR_ABS_TOL));
              sti_bdf1_param.set<double>("error rel tol", num_list.get<double>("steady error rel tol",ST_ERROR_REL_TOL));
              sti_bdf1_param.set<int>("max divergent iterations",
                                      num_list.get<int>("steady max divergent iterations",ST_MAX_DIVERGENT_ITERATIONS));
              sti_bdf1_param.set<double>("nonlinear iteration damping factor",
                                         num_list.get<double>("steady nonlinear iteration damping factor",ST_NONLIN_DAMP));
              sti_bdf1_param.set<int>("nonlinear iteration initial guess extrapolation order",
                                      num_list.get<int>("steady nonlinear iteration initial guess extrapolation order",ST_NONLIN_INIT_GUESS_EXTR_ORD));
              sti_bdf1_param.set<double>("restart tolerance relaxation factor",
                                         num_list.get<double>("steady restart tolerance relaxation factor",ST_NONLIN_INIT_TS_FACTOR));
              sti_bdf1_param.set<double>("restart tolerance relaxation factor damping",
                                         num_list.get<double>("steady restart tolerance relaxation factor damping",ST_NONLIN_INIT_TS_FACTOR_DAMP));
              sti_bdf1_param.set<double>("nonlinear iteration divergence factor",
                                         num_list.get<double>("steady nonlinear iteration divergence factor",ST_DIVERG_FACT));

              steady_time_integrator.set<std::string>("preconditioner",
                                                      num_list.get<std::string>("steady preconditioner",ST_PRECOND));
              steady_time_integrator.set<bool>("initialize with darcy",
                                               num_list.get<bool>("steady initialize with darcy",ST_INIT_DARCY));
            }
          }
        }


        // crerate sublists for the transient time integrator
        Teuchos::ParameterList& transient_time_integrator = richards_problem.sublist("transient time integrator");
        transient_time_integrator.set<std::string>("time integration method", "BDF1");
        Teuchos::ParameterList& tti_bdf1 = transient_time_integrator.sublist("BDF1");
        Teuchos::ParameterList& tti_bdf1_param = tti_bdf1.sublist("BDF1 parameters");

        transient_time_integrator.set<std::string>("preconditioner", TR_PRECOND);
        transient_time_integrator.set<std::string>("linear solver", TR_SOLVER);

        // set some probably not so good defaults for the steady computation
        tti_bdf1_param.set<int>("max iterations",TR_MAX_ITER);
        tti_bdf1_param.set<int>("min iterations",TR_MIN_ITER);
        tti_bdf1_param.set<int>("limit iterations",TR_LIMIT_ITER);
        tti_bdf1_param.set<double>("nonlinear tolerance",TRANSIENT_NONLINEAR_TOLERANCE);
        tti_bdf1_param.set<double>("time step reduction factor",TR_TS_RED_FACTOR);
        tti_bdf1_param.set<double>("time step increase factor",TR_TS_INC_FACTOR);
        tti_bdf1_param.set<double>("max time step", TR_MAX_TS);
        tti_bdf1_param.set<int>("max preconditioner lag iterations", TR_MAX_PREC_LAG);
        tti_bdf1_param.set<double>("error abs tol", TR_ERROR_ABS_TOL);
        tti_bdf1_param.set<double>("error rel tol", TR_ERROR_REL_TOL);
        tti_bdf1_param.set<int>("max divergent iterations",TR_MAX_DIVERGENT_ITERATIONS);
        tti_bdf1_param.set<double>("nonlinear iteration damping factor",TR_NONLIN_DAMP);
        tti_bdf1_param.set<int>("nonlinear iteration initial guess extrapolation order",TR_NONLIN_INIT_GUESS_EXTR_ORD);
        tti_bdf1_param.set<double>("restart tolerance relaxation factor",TR_NONLIN_INIT_TS_FACTOR);
        tti_bdf1_param.set<double>("restart tolerance relaxation factor damping",TR_NONLIN_INIT_TS_FACTOR_DAMP);
        tti_bdf1_param.set<double>("nonlinear iteration divergence factor",TR_DIVERG_FACT);

        if (plist->sublist("Execution Control").isSublist("Numerical Control Parameters")) {
          Teuchos::ParameterList& ncp_list = plist->sublist("Execution Control").sublist("Numerical Control Parameters");
          if (ncp_list.isSublist("Unstructured Algorithm")) {
            Teuchos::ParameterList& ncpu_list = ncp_list.sublist("Unstructured Algorithm");
            if (ncpu_list.isSublist("Transient Implicit Time Integration")) {

              Teuchos::ParameterList& num_list = ncpu_list.sublist("Transient Implicit Time Integration");
              tti_bdf1_param.set<int>("max iterations", num_list.get<int>("transient max iterations",TR_MAX_ITER));
              tti_bdf1_param.set<int>("min iterations", num_list.get<int>("transient min iterations",TR_MIN_ITER));
              tti_bdf1_param.set<int>("limit iterations", num_list.get<int>("transient limit iterations",TR_LIMIT_ITER));
              tti_bdf1_param.set<double>("nonlinear tolerance",
                                         num_list.get<double>("transient nonlinear tolerance",TRANSIENT_NONLINEAR_TOLERANCE));
              tti_bdf1_param.set<double>("time step reduction factor",
                                         num_list.get<double>("transient time step reduction factor",TR_TS_RED_FACTOR));
              tti_bdf1_param.set<double>("time step increase factor",
                                         num_list.get<double>("transient time step increase factor",TR_TS_INC_FACTOR));
              tti_bdf1_param.set<double>("max time step", num_list.get<double>("transient max time step",TR_MAX_TS));
              tti_bdf1_param.set<int>("max preconditioner lag iterations",
                                      num_list.get<int>("transient max preconditioner lag iterations",TR_MAX_PREC_LAG));
              tti_bdf1_param.set<double>("error abs tol", num_list.get<double>("transient error abs tol",TR_ERROR_ABS_TOL));
              tti_bdf1_param.set<double>("error rel tol", num_list.get<double>("transient error rel tol",TR_ERROR_REL_TOL));
              tti_bdf1_param.set<int>("max divergent iterations",
                                      num_list.get<int>("transient max divergent iterations",TR_MAX_DIVERGENT_ITERATIONS));
              tti_bdf1_param.set<double>("nonlinear iteration damping factor",
                                         num_list.get<double>("transient nonlinear iteration damping factor",TR_NONLIN_DAMP));
              tti_bdf1_param.set<int>("nonlinear iteration initial guess extrapolation order",
                                      num_list.get<int>("transient nonlinear iteration initial guess extrapolation order",TR_NONLIN_INIT_GUESS_EXTR_ORD));
              tti_bdf1_param.set<double>("restart tolerance relaxation factor",
                                         num_list.get<double>("transient restart tolerance relaxation factor",TR_NONLIN_INIT_TS_FACTOR));
              tti_bdf1_param.set<double>("restart tolerance relaxation factor damping",
                                         num_list.get<double>("transient restart tolerance relaxation factor damping",TR_NONLIN_INIT_TS_FACTOR_DAMP));
              tti_bdf1_param.set<double>("nonlinear iteration divergence factor",
                                         num_list.get<double>("transient nonlinear iteration divergence factor",TR_DIVERG_FACT));

              transient_time_integrator.set<std::string>("preconditioner",
                                                         num_list.get<std::string>("transient preconditioner",TR_PRECOND));
            }
          }
        }

        transient_time_integrator.sublist("VerboseObject") = create_Verbosity_List(verbosity_level);

        // insert the water retention models sublist
        Teuchos::ParameterList &water_retention_models = richards_problem.sublist("Water retention models");
        water_retention_models = create_WRM_List(plist);

        // insert the flow BC sublist
        Teuchos::ParameterList& flow_bc = richards_problem.sublist("boundary conditions");
        flow_bc = create_SS_FlowBC_List(plist);

        // insert sources, if they exist
        Teuchos::ParameterList& flow_src = richards_problem.sublist("source terms");
        flow_src = create_FlowSrc_List(plist);

      } else {
        // something's wrong
      }
    }
  }

  flw_list.sublist("VerboseObject") = create_Verbosity_List(verbosity_level);

  return flw_list;
}

Teuchos::ParameterList create_FlowSrc_List(Teuchos::ParameterList* plist)
{
  Teuchos::ParameterList src_list;
  
  Teuchos::ParameterList& src_sublist = plist->sublist("Sources");

  for (Teuchos::ParameterList::ConstIterator i = src_sublist.begin(); i != src_sublist.end(); i++) {
    // look at sublists
    if (src_sublist.isSublist(src_sublist.name(i))) {
      Teuchos::ParameterList& src = src_sublist.sublist(src_sublist.name(i));
      // get name
      std::string name = src_sublist.name(i);
      // create src sublist
      Teuchos::ParameterList& src_sub_out = src_list.sublist(name);
      // get the regions
      Teuchos::Array<std::string> regions = src.get<Teuchos::Array<std::string> >("Assigned Regions");
      src_sub_out.set<Teuchos::Array<std::string> >("regions",regions);
      // get source function
      Teuchos::ParameterList src_fn; 
      if (src.isSublist("Source: Volume Weighted")) {
        src_sub_out.set<std::string>("source and sink distribution method","volume");
        src_fn = src.sublist("Source: Volume Weighted");
      } else if (src.isSublist("Source: Permeability Weighted")) {
        src_sub_out.set<std::string>("source and sink distribution method","permeability"); 
        src_fn = src.sublist("Source: Permeability Weighted");
      } else {
        src_sub_out.set<std::string>("source and sink distribution method","none");          
        // where does the function come from???
      }
      // create time function 
      Teuchos::ParameterList& src_sub_out_fn = src_sub_out.sublist("sink");
      Teuchos::Array<double> values = src_fn.get<Teuchos::Array<double> >("Values");
      // write the native time function
      if (values.size() == 1) {
        src_sub_out_fn.sublist("function-constant").set<double>("value",values[0]);
      } else if (values.size() > 1) {
       Teuchos::Array<double> times = src_fn.get<Teuchos::Array<double> >("Times");
       Teuchos::Array<std::string> time_fns =  src_fn.get<Teuchos::Array<std::string> >("Time Functions");
      
       Teuchos::ParameterList& ssofn = src_sub_out_fn.sublist("function-tabular");
        
        ssofn.set<Teuchos::Array<double> >("x values", times);
        ssofn.set<Teuchos::Array<double> >("y values", values);

        Teuchos::Array<std::string> forms_(time_fns.size());
          
        for (int i=0; i<time_fns.size(); i++) {
          if (time_fns[i] == "Linear") {
            forms_[i] = "linear";
          } else if (time_fns[i] == "Constant") {
            forms_[i] = "constant";
          } else {
            Exceptions::amanzi_throw(Errors::Message("In the definition of Sources: time function can only be 'Linear' or 'Constant'"));
          }
        }
        ssofn.set<Teuchos::Array<std::string> >("forms",forms_);
      } else {
        // something is wrong with the input
      }
    }
  } 

  src_list.print(std::cout);

  return src_list;
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
      double ell;
      if (rel_perm == "Mualem") {
        ell = vG_list.get<double>("ell",ELL_MUALEM);
      } else if (rel_perm == "Burdine") {
        ell = vG_list.get<double>("ell",ELL_BURDINE);
      }
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
      double ell;
      if (rel_perm == "Mualem") {
        ell = BC_list.get<double>("ell",ELL_MUALEM);
      } else if (rel_perm == "Burdine") {
        ell = BC_list.get<double>("ell",ELL_BURDINE);
      }
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

  double aggthr(ML_AGG_THR);
  std::string smthtyp(ML_SMOOTHER);
  int ncycles(ML_NCYC);
  int nsmooth(ML_NSMOOTH);

  if (plist->sublist("Execution Control").isSublist("Numerical Control Parameters")) {
    Teuchos::ParameterList& ncp_list = plist->sublist("Execution Control").sublist("Numerical Control Parameters");
    if (ncp_list.isSublist("Unstructured Algorithm")) {
      Teuchos::ParameterList& ncpu_list = ncp_list.sublist("Unstructured Algorithm");
      if (ncpu_list.isSublist("Preconditioners")) {
        Teuchos::ParameterList& ncpup_list = ncpu_list.sublist("Preconditioners");
        if (ncpup_list.isSublist("Trilinos ML")) {
          Teuchos::ParameterList& ml_list = ncpup_list.sublist("Trilinos ML");
          if (ml_list.isParameter("ML aggregation threshold")) {
            aggthr = ml_list.get<double>("ML aggregation threshold");
          }
          if (ml_list.isParameter("ML smoother type")) {
            smthtyp = ml_list.get<std::string>("ML smoother type");
          }
          if (ml_list.isParameter("ML cycle applications")) {
            ncycles = ml_list.get<int>("ML cycle applications");
          }
          if (ml_list.isParameter("ML smoother sweeps")) {
            nsmooth = ml_list.get<int>("ML smoother sweeps");
          }
        }
      }
    }
  }

  Teuchos::ParameterList& ml_list = dpc_list.sublist("ML Parameters");
  ml_list.set<int>("ML output", ML_OUTPUT);
  ml_list.set<int>("max levels", ML_MAXLVLS);
  ml_list.set<std::string>("prec type",ML_PRECTYPE);
  ml_list.set<int>("cycle applications", ncycles);
  ml_list.set<std::string>("aggregation: type", ML_AGGTYPE);
  ml_list.set<double>("aggregation: damping factor", ML_AGGDAMP);
  ml_list.set<double>("aggregation: threshold", aggthr);
  ml_list.set<std::string>("eigen-analysis: type",ML_EIGENANAL_TYPE);
  ml_list.set<int>("eigen-analysis: iterations", ML_EIGENANAL_ITERS);
  ml_list.set<int>("smoother: sweeps", nsmooth);
  ml_list.set<double>("smoother: damping factor", ML_SMOOTH_DAMP);
  ml_list.set<std::string>("smoother: pre or post", ML_SMOOTH_PRE_POST);
  ml_list.set<std::string>("smoother: type", smthtyp);
  ml_list.set<double>("smoother: damping factor", ML_SMOOTH_DAMP);
  ml_list.set<std::string>("coarse: type", ML_CSOLVE_TYPE);
  ml_list.set<int>("coarse: max size", ML_CSOLVE_MAX_SIZE);

  return dpc_list;
}

/* ******************************************************************
 * Block ILU preconditioner sublist
 ****************************************************************** */
Teuchos::ParameterList create_BILU_List(Teuchos::ParameterList* plist)
{
  Teuchos::ParameterList bilu_list;

  double bilu_relax_value(ILU_RLXVAL);
  double bilu_abs_thresh(ILU_ABSTHR);
  double bilu_rel_thresh(ILU_RELTHR);
  int bilu_level_of_fill(ILU_LVLFILL);
  int bilu_overlap(ILU_OLV);

  if (plist->sublist("Execution Control").isSublist("Numerical Control Parameters")) {
    Teuchos::ParameterList& ncp_list = plist->sublist("Execution Control").sublist("Numerical Control Parameters");
    if (ncp_list.isSublist("Unstructured Algorithm")) {
      Teuchos::ParameterList& ncpu_list = ncp_list.sublist("Unstructured Algorithm");
      if (ncpu_list.isSublist("Preconditioners")) {
        Teuchos::ParameterList& ncpup_list = ncpu_list.sublist("Preconditioners");
        if (ncpup_list.isSublist("Block ILU")) {
          Teuchos::ParameterList& ilu_list = ncpup_list.sublist("Block ILU");
          if (ilu_list.isParameter("Block ILU relax value")) {
            bilu_relax_value = ilu_list.get<double>("Block ILU relax value");
          }
          if (ilu_list.isParameter("Block ILU relative threshold")) {
            bilu_rel_thresh = ilu_list.get<double>("Block ILU relative threshold");
          }
          if (ilu_list.isParameter("Block ILU absolute threshold")) {
            bilu_abs_thresh = ilu_list.get<double>("Block ILU absolute threshold");
          }
          if (ilu_list.isParameter("Block ILU level of fill")) {
            bilu_level_of_fill = ilu_list.get<int>("Block ILU level of fill");
          }
          if (ilu_list.isParameter("Block ILU overlap")) {
            bilu_overlap = ilu_list.get<int>("Block ILU overlap");
          }
        }
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

  double tol(AMG_TOL);
  int ncycles(AMG_NCYC);
  int nsmooth(AMG_NSMOOTH);
  double strong_threshold(AMG_STR_THR);

  if (plist->sublist("Execution Control").isSublist("Numerical Control Parameters")) {
    Teuchos::ParameterList& ncp_list = plist->sublist("Execution Control").sublist("Numerical Control Parameters");
    if (ncp_list.isSublist("Unstructured Algorithm")) {
      Teuchos::ParameterList& ncpu_list = ncp_list.sublist("Unstructured Algorithm");
      if (ncpu_list.isSublist("Preconditioners")) {
        Teuchos::ParameterList& ncpup_list = ncpu_list.sublist("Preconditioners");
        if (ncpup_list.isSublist("Hypre AMG")) {
          Teuchos::ParameterList& hypre_list = ncpup_list.sublist("Hypre AMG");
          if (hypre_list.isParameter("Hypre AMG cycle applications")) {
            ncycles = hypre_list.get<int>("Hypre AMG cycle applications");
          }
          if (hypre_list.isParameter("Hypre AMG smoother sweeps")) {
            nsmooth = hypre_list.get<int>("Hypre AMG smoother sweeps");
          }
          if (hypre_list.isParameter("Hypre AMG tolerance")) {
            tol = hypre_list.get<double>("Hypre AMG tolerance");
          }
          if (hypre_list.isParameter("Hypre AMG strong threshold")) {
            strong_threshold = hypre_list.get<double>("Hypre AMG strong threshold");
          }
        }
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

          for (int i=0; i<time_fns.size(); i++) {

            if (time_fns[i] == "Linear") {
              forms_[i] = "linear";
              //               values_plot[2*i + 1] = 0.5 *( flux[i] + flux[i+1]);
            } else if (time_fns[i] == "Constant") {
              forms_[i] = "constant";
              //               values_plot[2*i + 1] = flux[i];
            } else {
              Exceptions::amanzi_throw(Errors::Message("In the definition of BCs: tabular function can only be Linear or Constant"));
            }

          }

          Teuchos::Array<std::string> forms = forms_;
          tbcs.set<Teuchos::Array<std::string> >("forms", forms);

          //           for (int i=0; i < times_plot.size(); i++)
          //                 std::cout << times_plot[i] <<" "<<values_plot[i] <<"\n";
          //           }

          //           exit(0);

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

  stt_list.set<double>("atmospheric pressure",ATMOSPHERIC_PRESSURE);

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
      
      double specific_yield;
      bool use_specific_yield(false);
      if (matprop_list.sublist(matprop_list.name(i)).isSublist("Specific Yield: Uniform")) {
	specific_yield = matprop_list.sublist(matprop_list.name(i)).sublist("Specific Yield: Uniform").get<double>("Value");
	use_specific_yield = true;
      }

      double specific_storage;
      bool use_specific_storage(false);
      if (matprop_list.sublist(matprop_list.name(i)).isSublist("Specific Storage: Uniform")) {
	specific_storage = matprop_list.sublist(matprop_list.name(i)).sublist("Specific Storage: Uniform").get<double>("Value");
	use_specific_storage = true;
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

	if (use_specific_storage) stt_mat.set<double>("Constant specific storage",specific_storage);
	if (use_specific_yield)   stt_mat.set<double>("Constant specific yield",specific_yield);


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


void output_boundary_conditions( Teuchos::ParameterList* plist ) {

  Teuchos::ParameterList flow_list;
  int bc_counter = 0;

  if ( plist->isSublist("Flow") ) {
    flow_list =  plist->sublist("Flow");

    Teuchos::ParameterList richards_list;
    if ( flow_list.isSublist("Richards Problem") ) {
      richards_list = flow_list.sublist("Richards Problem");

      Teuchos::ParameterList bc_list;
      if ( richards_list.isSublist("boundary conditions") ) {
        bc_list = richards_list.sublist("boundary conditions");
      }

      Teuchos::ParameterList mass_flux_list, pressure_list, seepage_list, head_list;
      if (bc_list.isSublist("mass flux") ) {
        mass_flux_list = bc_list.sublist("mass flux");
        for (Teuchos::ParameterList::ConstIterator i = mass_flux_list.begin(); i != mass_flux_list.end(); i++) {

          if (mass_flux_list.isSublist(mass_flux_list.name(i))) {
            Teuchos::ParameterList& bc = mass_flux_list.sublist(mass_flux_list.name(i));

            if ((bc.sublist("outward mass flux")).isSublist("function-tabular")) {
              std::stringstream ss;
              ss << "BCmassflux" << bc_counter++;

              Teuchos::ParameterList& f_tab = (bc.sublist("outward mass flux")).sublist("function-tabular");

              Teuchos::Array<double> times = f_tab.get<Teuchos::Array<double> >("x values");
              Teuchos::Array<double> values = f_tab.get<Teuchos::Array<double> >("y values");
              Teuchos::Array<std::string> time_fns = f_tab.get<Teuchos::Array<std::string> >("forms");

              int np = times.size()*2 - 1;
              Teuchos::Array<double> times_plot(np);
              Teuchos::Array<double> values_plot(np);


              for (int i=0; i < times.size() - 1; i++) {
                times_plot[2*i] = times[i];
                values_plot[2*i] = values[i];
                times_plot[2*i + 1] = 0.5*(times[i] + times[i+1]);
              }
              times_plot[np - 1] = times[times.size() - 1];
              values_plot[np - 1] = values[times.size() - 1];

              for (int i=0; i<time_fns.size(); i++) {
                if (time_fns[i] == "linear") {
                  values_plot[2*i + 1] = 0.5 *( values[i] + values[i+1]);
                } else if (time_fns[i] == "constant") {
                  values_plot[2*i + 1] = values[i];
                  times_plot[2*i + 1] = times[i+1];
                } else {
                  Exceptions::amanzi_throw(Errors::Message("In the definition of BCs: tabular function can only be Linear or Constant"));
                }
              }


              std::string filename = ss.str() + ".dat";
              std::ofstream ofile(filename.c_str());

              ofile << "# "<<"time "<< "flux"<<std::endl;
              for (int i=0; i < np; i++) {
                ofile <<times_plot[i] << " " << values_plot[i] << std::endl;
              }

              ofile.close();
            }
          }
        }
      }
      if (bc_list.isSublist("pressure") ) {
        pressure_list = bc_list.sublist("pressure");
        for (Teuchos::ParameterList::ConstIterator i = pressure_list.begin(); i != pressure_list.end(); i++) {

          if (pressure_list.isSublist(pressure_list.name(i))) {
            Teuchos::ParameterList& bc = pressure_list.sublist(pressure_list.name(i));
            if ((bc.sublist("boundary pressure")).isSublist("function-tabular")) {
              std::stringstream ss;
              ss << "BCpressure" << bc_counter++;


              Teuchos::ParameterList& f_tab = (bc.sublist("boundary pressure")).sublist("function-tabular");

              Teuchos::Array<double> times = f_tab.get<Teuchos::Array<double> >("x values");
              Teuchos::Array<double> values = f_tab.get<Teuchos::Array<double> >("y values");
              Teuchos::Array<std::string> time_fns = f_tab.get<Teuchos::Array<std::string> >("forms");

              int np = times.size()*2 - 1;
              Teuchos::Array<double> times_plot(np);
              Teuchos::Array<double> values_plot(np);


              for (int i=0; i < times.size() - 1; i++) {
                times_plot[2*i] = times[i];
                values_plot[2*i] = values[i];
                times_plot[2*i + 1] = 0.5*(times[i] + times[i+1]);
              }
              times_plot[np - 1] = times[times.size() - 1];
              values_plot[np - 1] = values[times.size() - 1];

              for (int i=0; i<time_fns.size(); i++) {
                if (time_fns[i] == "linear") {
                  values_plot[2*i + 1] = 0.5 *( values[i] + values[i+1]);
                } else if (time_fns[i] == "constant") {
                  values_plot[2*i + 1] = values[i];
                  times_plot[2*i + 1] = times[i+1];
                } else {
                  Exceptions::amanzi_throw(Errors::Message("In the definition of BCs: tabular function can only be Linear or Constant"));
                }
              }


              std::string filename = ss.str() + ".dat";
              std::ofstream ofile(filename.c_str());

              ofile << "# time "<<"pressure"<<std::endl;
              for (int i=0; i < np; i++) {
                ofile << times_plot[i] << " " << values_plot[i] << std::endl;
              }

              ofile.close();
            }
          }
        }
      }
      if (bc_list.isSublist("seepage face") ) {
        seepage_list = bc_list.sublist("seepage face");
        for (Teuchos::ParameterList::ConstIterator i = seepage_list.begin(); i != seepage_list.end(); i++) {

          if (seepage_list.isSublist(seepage_list.name(i))) {
            Teuchos::ParameterList& bc = seepage_list.sublist(seepage_list.name(i));
            if ((bc.sublist("outward mass flux")).isSublist("function-tabular")) {
              std::stringstream ss;
              ss << "BCseepage" << bc_counter++;


              Teuchos::ParameterList& f_tab = (bc.sublist("outward mass flux")).sublist("function-tabular");

              Teuchos::Array<double> times = f_tab.get<Teuchos::Array<double> >("x values");
              Teuchos::Array<double> values = f_tab.get<Teuchos::Array<double> >("y values");
              Teuchos::Array<std::string> time_fns = f_tab.get<Teuchos::Array<std::string> >("forms");

              int np = times.size()*2 - 1;
              Teuchos::Array<double> times_plot(np);
              Teuchos::Array<double> values_plot(np);


              for (int i=0; i < times.size() - 1; i++) {
                times_plot[2*i] = times[i];
                values_plot[2*i] = values[i];
                times_plot[2*i + 1] = 0.5*(times[i] + times[i+1]);
              }
              times_plot[np - 1] = times[times.size() - 1];
              values_plot[np - 1] = values[times.size() - 1];

              for (int i=0; i<time_fns.size(); i++) {
                if (time_fns[i] == "linear") {
                  values_plot[2*i + 1] = 0.5 *( values[i] + values[i+1]);
                } else if (time_fns[i] == "constant") {
                  values_plot[2*i + 1] = values[i];
                  times_plot[2*i + 1] = times[i+1];
                } else {
                  Exceptions::amanzi_throw(Errors::Message("In the definition of BCs: tabular function can only be Linear or Constant"));
                }
              }


              std::string filename = ss.str() + ".dat";
              std::ofstream ofile(filename.c_str());

              ofile << "# time "<<"flux"<<std::endl;
              for (int i=0; i < np; i++) {
                ofile << times_plot[i] << " " << values_plot[i] << std::endl;
              }

              ofile.close();
            }
          }
        }
      }
      if (bc_list.isSublist("static head") ) {
        head_list = bc_list.sublist("static head");
        for (Teuchos::ParameterList::ConstIterator i = head_list.begin(); i != head_list.end(); i++) {

          if (head_list.isSublist(head_list.name(i))) {
            Teuchos::ParameterList& bc = head_list.sublist(head_list.name(i));
            if ((bc.sublist("water table elevation")).isSublist("function-tabular")) {
              std::stringstream ss;
              ss << "BChead" << bc_counter++;


              Teuchos::ParameterList& f_tab = (bc.sublist("water table elevation")).sublist("function-tabular");

              Teuchos::Array<double> times = f_tab.get<Teuchos::Array<double> >("x values");
              Teuchos::Array<double> values = f_tab.get<Teuchos::Array<double> >("y values");
              Teuchos::Array<std::string> time_fns = f_tab.get<Teuchos::Array<std::string> >("forms");

              int np = times.size()*2 - 1;
              Teuchos::Array<double> times_plot(np);
              Teuchos::Array<double> values_plot(np);


              for (int i=0; i < times.size() - 1; i++) {
                times_plot[2*i] = times[i];
                values_plot[2*i] = values[i];
                times_plot[2*i + 1] = 0.5*(times[i] + times[i+1]);
              }
              times_plot[np - 1] = times[times.size() - 1];
              values_plot[np - 1] = values[times.size() - 1];

              for (int i=0; i<time_fns.size(); i++) {
                if (time_fns[i] == "linear") {
                  values_plot[2*i + 1] = 0.5 *( values[i] + values[i+1]);
                } else if (time_fns[i] == "constant") {
                  values_plot[2*i + 1] = values[i];
                  times_plot[2*i + 1] = times[i+1];
                } else {
                  Exceptions::amanzi_throw(Errors::Message("In the definition of BCs: tabular function can only be Linear or Constant"));
                }
              }


              std::string filename = ss.str() + ".dat";
              std::ofstream ofile(filename.c_str());

              ofile << "# time "<<"head"<<std::endl;
              for (int i=0; i < np; i++) {
                ofile << times_plot[i] << " " << values_plot[i] << std::endl;
              }

              ofile.close();
            }
          }
        }
      }
    }
  }
  //         exit(0);
}

void check_AmanziInputVersion(Teuchos::ParameterList* plist) {

  std::string version = plist->get<std::string>("Amanzi Input Format Version","FAIL");
  if (version == "FAIL") {
    Exceptions::amanzi_throw(Errors::Message("The input file does not specify an input format version"));
  }

  int major, minor, micro;

  std::stringstream ss;
  ss << version;
  std::string ver;

  try {
    getline(ss,ver,'.');
    major = boost::lexical_cast<int>(ver);

    getline(ss,ver,'.');
    minor = boost::lexical_cast<int>(ver);

    getline(ss,ver);
    micro = boost::lexical_cast<int>(ver);
  }
  catch (...) {
    Exceptions::amanzi_throw(Errors::Message("The version string in the input file '"+version+"' has the wrong format, please use X.Y.Z, where X, Y, and Z are integers."));
  }

  if ((major != AMANZI_INPUT_VERSION_MAJOR) ||
      (minor != AMANZI_INPUT_VERSION_MINOR) ||
      (micro != AMANZI_INPUT_VERSION_MICRO)) {
    std::stringstream ss_ver_reqd;
    ss_ver_reqd << AMANZI_INPUT_VERSION_MAJOR << "." << AMANZI_INPUT_VERSION_MINOR << "." << AMANZI_INPUT_VERSION_MICRO;
    std::stringstream ss_ver_inp;
    ss_ver_inp << major << "." << minor << "." << micro;

    Exceptions::amanzi_throw(Errors::Message("The input format version "+ss_ver_inp.str()+" does not match the required version "+ss_ver_reqd.str()));
  }
}





}
}
