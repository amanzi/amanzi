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


template <typename T>
T amanzi_get(Teuchos::ParameterList& plist, const std::string &name) {
  if (plist.isParameter(name)) {
    return plist.get<T>(name);
  } else {
    Exceptions::amanzi_throw(Errors::Message("Trying to read required parameter '" + name + "'  "));
  }
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

  // walkabout list is optional
  tmp_list = create_Walkabout_Data_List(plist);
  if (tmp_list.begin() != tmp_list.end()) {
    new_list.sublist("Walkabout Data") = tmp_list;
  }

  // visualization list is optional
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


Teuchos::ParameterList get_Cycle_Macro(const std::string& macro_name, Teuchos::ParameterList* plist) {
  Teuchos::ParameterList cycle_macro;

  if ( plist->sublist("Output").sublist("Cycle Macros").isSublist(macro_name) ) {
    if (plist->sublist("Output").sublist("Cycle Macros").sublist(macro_name).isParameter("Start_Period_Stop")) {
      Teuchos::Array<int> cycle_range;
      cycle_range = plist->sublist("Output").sublist("Cycle Macros").sublist(macro_name)
          .get<Teuchos::Array<int> >("Start_Period_Stop");

      cycle_macro.set<Teuchos::Array<int> >("Start_Period_Stop",cycle_range);

    }
    if (plist->sublist("Output").sublist("Cycle Macros").sublist(macro_name).isParameter("Values")) {
      Teuchos::Array<int> values;
      values = plist->sublist("Output").sublist("Cycle Macros").sublist(macro_name)
          .get<Teuchos::Array<int> >("Values");
      cycle_macro.set<Teuchos::Array<int> >("Values",values);
    }
  } else {
    std::stringstream ss;
    ss << "The cycle macro " << macro_name << " does not exist in the input file";
    Exceptions::amanzi_throw(Errors::Message(ss.str().c_str()));
  }

  return cycle_macro;
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

  // create an all region
  if (!plist->sublist("Regions").isSublist("All")) {
    Teuchos::ParameterList& allreg = plist->sublist("Regions").sublist("All")
        .sublist("Region: Box");

    Teuchos::Array<double> low, high;
    low.push_back(-1e99);
    high.push_back(1e99);

    if (spatial_dimension_ >= 2) {
      low.push_back(-1e99);
      high.push_back(1e99);
    }

    if (spatial_dimension_ == 3) {
      low.push_back(-1e99);
      high.push_back(1e99);
    }

    allreg.set<Teuchos::Array<double> >("Low Coordinate", low);
    allreg.set<Teuchos::Array<double> >("High Coordinate", high);
  }

  // check if Transport is Off
  std::string transport_model = plist->sublist("Execution Control").get<std::string>("Transport Model");
  std::string chemistry_model = plist->sublist("Execution Control").get<std::string>("Chemistry Model");

  phase_name = "Aqueous";

  // don't know the history of these variables, clear them just to be safe.
  comp_names.clear();
  mineral_names_.clear();
  sorption_site_names_.clear();

  Teuchos::ParameterList& phase_list = plist->sublist("Phase Definitions");
  Teuchos::ParameterList::ConstIterator item;
  for (item = phase_list.begin(); item != phase_list.end(); ++item) {
    if (transport_model != "Off"  || chemistry_model != "Off" ) {
      if (phase_list.name(item) == "Aqueous" ) {
        Teuchos::ParameterList aqueous_list = phase_list.sublist("Aqueous");
        if (aqueous_list.isSublist("Phase Components")) {
          Teuchos::ParameterList phase_components = aqueous_list.sublist("Phase Components");
          // for now there should only be one sublist here, we allow it to be named something
          // the user chooses, e.g. Water
          Teuchos::ParameterList::ConstIterator pcit = phase_components.begin();
          ++pcit;
          if (pcit != phase_components.end()) {
            Exceptions::amanzi_throw(Errors::Message("Currently Amanzi only supports one phase component, e.g. Water"));
          }
          pcit = phase_components.begin();
          if (!pcit->second.isList()) {
            Exceptions::amanzi_throw(Errors::Message("The Phase Components list must only have one sublist, but you have specified a parameter instead."));
          }
          phase_comp_name = pcit->first;
          Teuchos::ParameterList& water_components = phase_components.sublist(phase_comp_name);
          if ( water_components.isParameter("Component Solutes")) {
            comp_names = water_components.get<Teuchos::Array<std::string> >("Component Solutes");
          }
        }  // end phase components
      }  // end Aqueous phase
    }
    if (chemistry_model != "Off") {
      if (phase_list.name(item) == "Solid") {
        Teuchos::ParameterList solid_list = phase_list.sublist("Solid");
        // this is the order that the chemistry expects
        if (solid_list.isParameter("Minerals")) {
          mineral_names_ = solid_list.get<Teuchos::Array<std::string> >("Minerals");
        }
        if (solid_list.isParameter("Sorption Sites")) {
          sorption_site_names_ = solid_list.get<Teuchos::Array<std::string> >("Sorption Sites");
        }
      }  // end Solid phase
    }
    if ( (phase_list.name(item) != "Aqueous" ) && (phase_list.name(item) != "Solid") ) {
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
  }

  if ( plist->isSublist("Execution Control") ) {

    std::string verbosity = plist->sublist("Execution Control").get<std::string>("Verbosity",VERBOSITY_DEFAULT);

    if ( verbosity == "None" || verbosity == "none") {
      verbosity_level = "none";
    } else if ( verbosity == "Low" || verbosity == "low") {
      verbosity_level = "low";
    } else if ( verbosity == "Medium" ||verbosity == "medium") {
      verbosity_level = "medium";
    } else if ( verbosity == "High" || verbosity == "high") {
      verbosity_level = "high";
    } else if ( verbosity == "Extreme" || verbosity == "extreme") {
      verbosity_level = "extreme";
    } else {
      Exceptions::amanzi_throw(Errors::Message("Verbosity must be one of None, Low, Medium, High, or Extreme."));
    }
  }

  // dispersion (this is going to be used to translate to the transport list as well as the state list)
  // check if we need to write a dispersivity sublist
  need_dispersion_ = false;
  if (plist->isSublist("Material Properties")) {
    for (Teuchos::ParameterList::ConstIterator it = plist->sublist("Material Properties").begin();
         it != plist->sublist("Material Properties").end(); ++it) {
      if ( (it->second).isList()) {
        Teuchos::ParameterList & mat_sublist = plist->sublist("Material Properties").sublist(it->first);
        for (Teuchos::ParameterList::ConstIterator jt = mat_sublist.begin(); jt != mat_sublist.end(); ++jt) {
          if ( (jt->second).isList() ) {
            const std::string pname = jt->first;
            if ( pname.find("Dispersion Tensor") == 0 ||
                 pname.find("Molecular Diffusion") == 0 ||
                 pname.find("Tortuosity") == 0 ) {
              need_dispersion_ = true;
            }
          }
        }
      }
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
      Teuchos::ParameterList rlist = plist->sublist("Output").sublist("Checkpoint Data");

      restart_list.set<std::string>("file name base", rlist.get<std::string>("File Name Base",std::string("checkpoint")));
      restart_list.set<int>("file name digits", rlist.get<int>("File Name Digits",5));

      // check if the cycle range is defined via a macro
      if ( rlist.isParameter("Cycle Macro") ) {
        std::string cycle_macro = rlist.get<std::string>("Cycle Macro");

        Teuchos::Array<int> range = get_Cycle_Macro(cycle_macro, plist).get<Teuchos::Array<int> >("Start_Period_Stop");

        restart_list.set<Teuchos::Array<int> >("cycles start period stop", range);
      }
    }
  }

  return restart_list;
}

/* ******************************************************************
 * Empty
 ****************************************************************** */
Teuchos::ParameterList create_Walkabout_Data_List(Teuchos::ParameterList* plist) {

  Teuchos::ParameterList walkabout_list;

  if ( plist->isSublist("Output") ) {

    if ( plist->sublist("Output").isSublist("Walkabout Data") ) {
      Teuchos::ParameterList rlist = plist->sublist("Output").sublist("Walkabout Data");

      walkabout_list.set<std::string>("file name base", rlist.get<std::string>("File Name Base",std::string("walkabout")));
      walkabout_list.set<int>("file name digits", rlist.get<int>("File Name Digits",5));

      // check if the cycle range is defined via a macro
      if ( rlist.isParameter("Cycle Macro") ) {
        std::string cycle_macro = rlist.get<std::string>("Cycle Macro");

        Teuchos::Array<int> range = get_Cycle_Macro(cycle_macro, plist).get<Teuchos::Array<int> >("Start_Period_Stop");

        walkabout_list.set<Teuchos::Array<int> >("cycles start period stop", range);
      }
    }
  }

  return walkabout_list;
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

      // file name
      if (vis_list.isParameter("File Name Base")) {
        vis_list.set<std::string>("file name base", vis_list.get<std::string>("File Name Base"));
        vis_list.remove("File Name Base");
      }
      if (vis_list.isParameter("File Name Digits")) {
        vis_list.remove("File Name Digits");
      }

      // Cycle Macro
      Teuchos::Array<int> all_cycles;
      all_cycles.clear();
      if ( vis_list.isParameter("Cycle Macros") ) {
        std::vector<std::string> cycle_macros;
       	cycle_macros = vis_list.get<Teuchos::Array<std::string> >("Cycle Macros").toVector();

	int j(0);
        for (int i=0; i < cycle_macros.size(); i++) {
          //Teuchos::Array<int> cm = get_Cycle_Macro(cycle_macro,plist);
          Teuchos::ParameterList cycle_macro_list = get_Cycle_Macro(cycle_macros[i], plist);
          if (cycle_macro_list.isParameter("Start_Period_Stop")) {
            std::stringstream ss;
            ss << "cycles start period stop " << j;
            vis_list.set(ss.str(),cycle_macro_list.get<Teuchos::Array<double> >("Start_Period_Stop"));
            ++j;
          } else if (cycle_macro_list.isParameter("Values")) {
            Teuchos::Array<int> cycles;
            cycles = cycle_macro_list.get<Teuchos::Array<int> >("Values");

            std::list<int> all_list, cur_list;
            for (Teuchos::Array<int>::iterator at = all_cycles.begin(); at != all_cycles.end(); ++at) {
              all_list.push_back(*at);
            }
            for (Teuchos::Array<int>::iterator t = cycles.begin(); t != cycles.end(); ++t) {
              cur_list.push_back(*t);
            }
            all_list.sort();
            cur_list.sort();

            all_list.merge(cur_list);
            all_list.unique();

            all_cycles.clear();
            for (std::list<int>::iterator al = all_list.begin(); al != all_list.end(); ++al) {
              all_cycles.push_back(*al);
            }
          } else {
            Exceptions::amanzi_throw(Errors::Message("Cycle Macros must hace either the Values of Start_Period_Stop parameter."));
          }
	}
        // delete the cycle macro
        vis_list.remove("Cycle Macros");
      }

      // Time Macros
      Teuchos::Array<double> all_times;
      all_times.clear();
      if ( vis_list.isParameter("Time Macros") ) {
        std::vector<std::string> time_macros;
        time_macros = vis_list.get<Teuchos::Array<std::string> >("Time Macros").toVector();

        int j(0);
        for (int i=0; i < time_macros.size(); i++) {
          // Create a local parameter to store the time macro
          Teuchos::ParameterList time_macro_list = get_Time_Macro(time_macros[i], plist);
          if (time_macro_list.isParameter("Start_Period_Stop")) {
            std::stringstream ss;
            ss << "times start period stop " << j;
            vis_list.set(ss.str(),time_macro_list.get<Teuchos::Array<double> >("Start_Period_Stop"));
            ++j;
          }
          if (time_macro_list.isParameter("Values")) {
            Teuchos::Array<double> times;
            times = time_macro_list.get<Teuchos::Array<double> >("Values");

            std::list<double> all_list, cur_list;
            for (Teuchos::Array<double>::iterator at = all_times.begin(); at != all_times.end(); ++at) {
              all_list.push_back(*at);
            }
            for (Teuchos::Array<double>::iterator t = times.begin(); t != times.end(); ++t) {
              cur_list.push_back(*t);
            }
            all_list.sort();
            cur_list.sort();

            all_list.merge(cur_list);
            all_list.unique();

            all_times.clear();
            for (std::list<double>::iterator al = all_list.begin(); al != all_list.end(); ++al) {
              all_times.push_back(*al);
            }
          }
        }
        vis_list.remove("Time Macros");
      }
      if (all_times.size() != 0) {
        vis_list.set("times", all_times);
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
              obs_list.sublist(i->first).set("times start period stop", time_macro_list.get<Teuchos::Array<double> >("Start_Period_Stop"));
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

            Teuchos::ParameterList cycle_macro_list = get_Cycle_Macro(cycle_macro, plist);

            Teuchos::Array<int> sps, values;

            if (cycle_macro_list.isParameter("Start_Period_Stop")) {
              sps = cycle_macro_list.get<Teuchos::Array<int> >("Start_Period_Stop");
            }
            if (cycle_macro_list.isParameter("Values")) {
              values = cycle_macro_list.get<Teuchos::Array<int> >("Values");
            }
            if (sps.size() != 3  && values.size() == 0) {
              Errors::Message message("Cycle macro " + cycle_macro + " has neither a valid Start_Period_Stop nor a valid Values parameter");
              Exceptions::amanzi_throw(message);
            }

            if (sps.size() == 3)
              obs_list.sublist(i->first).set("cycles start period stop", sps);
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
Teuchos::ParameterList create_TimePeriodControl_List(Teuchos::ParameterList* plist) {
  Teuchos::ParameterList tpc_list;

  Teuchos::Array<double> start_times;
  Teuchos::Array<double> initial_time_step;

  std::map<double, double> time_map;
  double default_initial_time_step(RESTART_TIME_STEP);

  if ( plist->isSublist("Execution Control") ) {
    Teuchos::ParameterList exe_sublist = plist->sublist("Execution Control");

    // get the default initial time step
    if (exe_sublist.isSublist("Time Period Control")) {
      default_initial_time_step = exe_sublist.sublist("Time Period Control").get<double>("Default Initial Time Step",RESTART_TIME_STEP);
    }

    if ( exe_sublist.isParameter("Flow Model") ) {
      std::string flow_model = exe_sublist.get<std::string>("Flow Model");
      if (flow_model != "Off") { // we need to process boudary conditions to possibly add additional time periods

        Teuchos::ParameterList& bc_sublist = plist->sublist("Boundary Conditions");

        for (Teuchos::ParameterList::ConstIterator i = bc_sublist.begin(); i != bc_sublist.end(); i++) {
          // look at sublists
          if (bc_sublist.isSublist(bc_sublist.name(i))) {
            Teuchos::ParameterList& bc = bc_sublist.sublist(bc_sublist.name(i));

            Teuchos::ParameterList bc_list;

            if ( bc.isSublist("BC: Flux") ) {
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
      }
    }
  }

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
Teuchos::ParameterList create_MPC_List(Teuchos::ParameterList* plist) {
  Teuchos::ParameterList mpc_list;

  bool transport_on(false), chemistry_on(false), flow_on(false);

  if ( plist->isSublist("Execution Control") ) {
    Teuchos::ParameterList exe_sublist = plist->sublist("Execution Control");

    mpc_list.sublist("Time Integration Mode") = exe_sublist.sublist("Time Integration Mode");

    //if (exe_sublist.isSublist("Time Period Control")) {
    mpc_list.sublist("Time Period Control") = create_TimePeriodControl_List(plist);
    //}

    // now interpret the modes
    if ( exe_sublist.isParameter("Transport Model") ) {
      if ( exe_sublist.get<std::string>("Transport Model") == "Off" || exe_sublist.get<std::string>("Transport Model") == "off") {
        mpc_list.set<std::string>("disable Transport_PK","yes");
      } else if ( exe_sublist.get<std::string>("Transport Model") == "On" || exe_sublist.get<std::string>("Transport Model") == "on") {
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

    if (transport_on || chemistry_on) {
      mpc_list.set<Teuchos::Array<std::string> >("component names", comp_names);
    }



    if ( exe_sublist.isParameter("Flow Model") ) {
      if ( exe_sublist.get<std::string>("Flow Model") == "Off" || exe_sublist.get<std::string>("Flow Model") == "off") {
        mpc_list.set<std::string>("disable Flow_PK", "yes");
        flow_on = false;
      } else if ( exe_sublist.get<std::string>("Flow Model") == "Richards" ) {
        mpc_list.set<std::string>("disable Flow_PK", "no");
        mpc_list.set<std::string>("Flow model","Richards");
        flow_on = true;
      } else if (exe_sublist.get<std::string>("Flow Model") == "Single Phase") {
        mpc_list.set<std::string>("disable Flow_PK", "no");
        mpc_list.set<std::string>("Flow model","Steady State Saturated");
        flow_on = true;
      } else if (exe_sublist.get<std::string>("Flow Model") == "Transient with Static Flow") {
	mpc_list.set<std::string>("disable Flow_PK", "no");
	mpc_list.set<std::string>("Flow model","Static Steady State Saturated");
	flow_on = true;
      } else {
        Exceptions::amanzi_throw(Errors::Message("Flow Model must either be Richards, Single Phase, Transient with Static Flow, or Off"));
      }
    } else {
      Exceptions::amanzi_throw(Errors::Message("The parameter Flow Model must be specified."));
    }

    if (flow_on) {
      double ti_rescue(TI_RESCUE_REDUCTION_FACTOR);

      if (plist->sublist("Execution Control").isSublist("Numerical Control Parameters")) {
        if (plist->sublist("Execution Control").sublist("Numerical Control Parameters").isSublist("Unstructured Algorithm")) {
          Teuchos::ParameterList& ncpu_list = plist->sublist("Execution Control").sublist("Numerical Control Parameters").sublist("Unstructured Algorithm");
          if (ncpu_list.isSublist("MPC")) {
            ti_rescue = ncpu_list.sublist("MPC").get<double>("time integration rescue reduction factor",TI_RESCUE_REDUCTION_FACTOR);
          }
        }
      }
      mpc_list.set<double>("time integration rescue reduction factor",TI_RESCUE_REDUCTION_FACTOR);
    }

    if (plist->sublist("Execution Control").isSublist("Restart from Checkpoint Data File") &&
        plist->sublist("Initial Conditions").isParameter("Init from Checkpoint File")) {
      // this is an error, you can either restart or re-init, but not both
      Exceptions::amanzi_throw(Errors::Message("You can either restart from a checkpoint or initialize from a checkpoint, but not both."));
    }

    if ( plist->sublist("Execution Control").isSublist("Restart from Checkpoint Data File") ) {
      mpc_list.sublist("Restart from Checkpoint Data File") =
          plist->sublist("Execution Control").sublist("Restart from Checkpoint Data File");
    }

    if ( plist->sublist("Initial Conditions").isParameter("Init from Checkpoint File")) {
      Teuchos::ParameterList & rest_list = mpc_list.sublist("Restart from Checkpoint Data File");
      std::string file = plist->sublist("Initial Conditions").get<std::string>("Init from Checkpoint File");
      rest_list.set<std::string>("Checkpoint Data File Name",file);
      rest_list.set<bool>("initialize from checkpoint data file and do not restart",true);
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

        // get the expert parameters
        double CFL(1.0);
        if (plist->sublist("Execution Control").isSublist("Numerical Control Parameters")) {
          if (plist->sublist("Execution Control").sublist("Numerical Control Parameters").isSublist("Unstructured Algorithm")) {
            if (plist->sublist("Execution Control").sublist("Numerical Control Parameters").sublist("Unstructured Algorithm").isSublist("Transport Process Kernel")) {
              Teuchos::ParameterList t_exp_params = plist->sublist("Execution Control").sublist("Numerical Control Parameters").sublist("Unstructured Algorithm").sublist("Transport Process Kernel");
              if (t_exp_params.isParameter("CFL")) {
                CFL = t_exp_params.get<double>("CFL");
              }
            }
          }
        }

        // transport is on, set some defaults
        trp_list.set<int>("spatial discretization order", 1);
        trp_list.set<int>("temporal discretization order", 1);
        trp_list.sublist("VerboseObject") = create_Verbosity_List(verbosity_level);
        trp_list.set<std::string>("enable internal tests", "no");
        trp_list.set<double>("CFL", CFL);
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
                  trp_list.set<int>("spatial discretization order", 1);
                  trp_list.set<int>("temporal discretization order", 1);
                } else if ( tia == "Explicit Second-Order" ) {
                  trp_list.set<int>("spatial discretization order", 2);
                  trp_list.set<int>("temporal discretization order", 2);
                }
              }
            }
          }
        }

        // now write the dispersion lists if needed
        if (need_dispersion_) {
          Teuchos::ParameterList &disp_list = trp_list.sublist("Dispersivity");

          if (plist->isSublist("Material Properties")) {
            for (Teuchos::ParameterList::ConstIterator it = plist->sublist("Material Properties").begin();
                 it != plist->sublist("Material Properties").end(); ++it) {
              disp_list.set<std::string>("numerical method","two point flux approximation");
              disp_list.set<std::string>("solver","PCG with Hypre AMG");

              if ( (it->second).isList()) {
                std::string mat_name = it->first;
                Teuchos::ParameterList & mat_sublist = plist->sublist("Material Properties").sublist(mat_name);
                Teuchos::ParameterList & disp_sublist = disp_list.sublist(mat_name);
                // set a few default paramters
                disp_sublist.set<std::string>("model","Bear");
                // translate the other paramters
                disp_sublist.set<Teuchos::Array<std::string> >("regions",mat_sublist.get<Teuchos::Array<std::string> >("Assigned Regions"));
                if (!mat_sublist.isSublist("Dispersion Tensor: Uniform Isotropic")) {
                  Exceptions::amanzi_throw(Errors::Message("Dispersion is enabled, you must specify Dispersion Tensor: Uniform Isotropic for all materials. Disable it by purging all Material Property sublists of the Dispersion Tensor:, Molecular Diffusion:, and Tortuosity: sublists."));
                }
                disp_sublist.set<double>("alphaL", mat_sublist.sublist("Dispersion Tensor: Uniform Isotropic").get<double>("alphaL"));
                disp_sublist.set<double>("alphaT", mat_sublist.sublist("Dispersion Tensor: Uniform Isotropic").get<double>("alphaT"));
                if (!mat_sublist.isSublist("Molecular Diffusion: Uniform")) {
                  Exceptions::amanzi_throw(Errors::Message("Dispersion is enabled, you must specify Molecular Diffusion: Uniform for all materials. Disable it by purging all Material Property sublists of the Dispersion Tensor:, Molecular Diffusion:, and Tortuosity: sublists."));
                }
                disp_sublist.set<double>("D", mat_sublist.sublist("Molecular Diffusion: Uniform").get<double>("Value"));
                if (!mat_sublist.isSublist("Tortuosity: Uniform")) {
                  Exceptions::amanzi_throw(Errors::Message("Dispersion is enabled, you must specify Tortuosity: Uniform for all materials. Disable it by purging all Material Property sublists of the Dispersion Tensor:, Molecular Diffusion:, and Tortuosity: sublists."));
                }
                disp_sublist.set<double>("tortuosity", mat_sublist.sublist("Tortuosity: Uniform").get<double>("Value"));
              }
            }
          }
        }

        // now generate the source lists
        Teuchos::ParameterList src_list = create_TransportSrc_List(plist);
        if (src_list.begin() != src_list.end()) { // the source lists are not empty
          trp_list.sublist("source terms") = src_list;
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
          Teuchos::ParameterList& tbc_list = trp_list.sublist("boundary conditions").sublist("concentration");

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
                      if ( comps.sublist(*i).isSublist("BC: Uniform Concentration") ) {
                        Teuchos::ParameterList& bc = tbc_list.sublist(compss.str()).sublist(bc_root_str);
                        bc.set<Teuchos::Array<std::string> >("regions",regs);
                        Teuchos::ParameterList& bcsub = comps.sublist(*i).sublist("BC: Uniform Concentration");

                        Teuchos::Array<double> values = bcsub.get<Teuchos::Array<double> >("Values");
                        Teuchos::Array<double> times = bcsub.get<Teuchos::Array<double> >("Times");
                        Teuchos::Array<std::string> time_fns = bcsub.get<Teuchos::Array<std::string> >("Time Functions");

                        Teuchos::ParameterList &bcfn = bc.sublist("boundary concentration").sublist("function-tabular");
                        bcfn.set<Teuchos::Array<double> >("y values", values);
                        bcfn.set<Teuchos::Array<double> >("x values", times);
                        bcfn.set<Teuchos::Array<std::string> >("forms", translate_forms(time_fns));
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
      }
    }
  }
  return trp_list;
}



Teuchos::Array<std::string> translate_forms(Teuchos::Array<std::string> & forms) {
  Teuchos::Array<std::string>  target_forms;
  for (Teuchos::Array<std::string>::const_iterator i = forms.begin();
       i != forms.end(); ++i) {

    if (*i == "Constant") {
      target_forms.push_back("constant");
    } else if (*i == "Linear") {
      target_forms.push_back("linear");
    } else {
      Exceptions::amanzi_throw(Errors::Message("Cannot translate the tabular function form "+*i));
    }
  }
  return target_forms;
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
  Teuchos::ParameterList& pcg_list = solver_list.sublist("PCG with Hypre AMG");

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
        if (num_list.isParameter("linear solver iterative method"))
          method = num_list.get<std::string>("linear solver iterative method");
        if (num_list.isParameter("linear solver method"))
          prec = num_list.get<std::string>("linear solver preconditioner");
      }
    }
  }
  aztecoo_list.set<double>("error tolerance", tol);
  aztecoo_list.set<std::string>("iterative method", method);
  aztecoo_list.set<int>("maximum number of iterations", maxiter);
  aztecoo_list.set<std::string>("preconditioner", prec);

  // add default PCG solver
  pcg_list.set<double>("error tolerance", tol);
  pcg_list.set<std::string>("iterative method", "pcg");
  pcg_list.set<int>("maximum number of iterations", maxiter);
  pcg_list.set<std::string>("preconditioner", prec);

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

      Teuchos::ParameterList *flow_list;

      // get the expert parameters
      std::string disc_method("optimized mfd scaled");
      std::string rel_perm("upwind with Darcy flux");
      double atm_pres(ATMOSPHERIC_PRESSURE);
      std::string nonlinear_solver("NKA");

      if (plist->sublist("Execution Control").isSublist("Numerical Control Parameters")) {
        if (plist->sublist("Execution Control").sublist("Numerical Control Parameters").isSublist("Unstructured Algorithm")) {
          if (plist->sublist("Execution Control").sublist("Numerical Control Parameters").sublist("Unstructured Algorithm").isSublist("Flow Process Kernel")) {
            Teuchos::ParameterList fl_exp_params = plist->sublist("Execution Control").sublist("Numerical Control Parameters").sublist("Unstructured Algorithm").sublist("Flow Process Kernel");
            if (fl_exp_params.isParameter("Discretization Method")) {
              disc_method = fl_exp_params.get<std::string>("Discretization Method");
            }
            if (fl_exp_params.isParameter("Relative Permeability")) {
              rel_perm = fl_exp_params.get<std::string>("Relative Permeability");
            }
            if (fl_exp_params.isParameter("atmospheric pressure")) {
              atm_pres = fl_exp_params.get<double>("atmospheric pressure");
            }
          }
          if (plist->sublist("Execution Control").sublist("Numerical Control Parameters").sublist("Unstructured Algorithm").isSublist("Nonlinear Solver")) {
            Teuchos::ParameterList fl_exp_params = plist->sublist("Execution Control").sublist("Numerical Control Parameters").sublist("Unstructured Algorithm").sublist("Nonlinear Solver");
            if (fl_exp_params.isParameter("Nonlinear Solver Type")) {
              nonlinear_solver = fl_exp_params.get<std::string>("Nonlinear Solver Type");
            }
          }
        }
      }
      // discretization method must be two point flux approximation for if newton is used
      if (nonlinear_solver == std::string("Newton") || nonlinear_solver == std::string("inexact Newton")) {
        disc_method = std::string("two point flux approximation");
      }

      if (flow_model == "Single Phase" || flow_model == "Richards" || flow_model == "Transient with Static Flow") {

        if (flow_model == "Single Phase" || flow_model == "Transient with Static Flow") {
          Teuchos::ParameterList& darcy_problem = flw_list.sublist("Darcy Problem");
          darcy_problem.sublist("VerboseObject") = create_Verbosity_List(verbosity_level);
          darcy_problem.set<double>("atmospheric pressure", atm_pres);
          darcy_problem.set<std::string>("discretization method", disc_method);

          flow_list = &darcy_problem; // we use this below to insert sublists that are shared by Richards and Darcy
        } else if (flow_model == "Richards") {
          Teuchos::ParameterList& richards_problem = flw_list.sublist("Richards Problem");
          richards_problem.set<std::string>("relative permeability", rel_perm);
          // this one should come from the input file...
          richards_problem.sublist("VerboseObject") = create_Verbosity_List(verbosity_level);
          richards_problem.set<double>("atmospheric pressure", atm_pres);
          richards_problem.set<std::string>("discretization method", disc_method);

          // see if we need to generate a Picard list

          flow_list = &richards_problem; // we use this below to insert sublists that are shared by Richards and Darcy

          // insert the water retention models sublist (these are only relevant for Richards)
          Teuchos::ParameterList &water_retention_models = richards_problem.sublist("Water retention models");
          water_retention_models = create_WRM_List(plist);
        }

        // insert the flow BC sublist
        Teuchos::ParameterList flow_bc; // = flow_list->sublist("boundary conditions");
        flow_bc = create_SS_FlowBC_List(plist);
        if ( flow_bc.begin() != flow_bc.end() ) {
          flow_list->sublist("boundary conditions") = flow_bc;
        }

        // insert sources, if they exist
        Teuchos::ParameterList flow_src; // = flow_list->sublist("source terms");
        flow_src = create_FlowSrc_List(plist);
        if (flow_src.begin() != flow_src.end()) {
          flow_list->sublist("source terms") = flow_src;
        }

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

          Teuchos::ParameterList& picard_list = flow_list->sublist("initial guess pseudo time integrator");
          picard_list.sublist("Picard").set<double>("error abs tol", ST_ERROR_ABS_TOL);
          picard_list.sublist("Picard").set<double>("error rel tol", ST_ERROR_REL_TOL);
          picard_list.sublist("Picard").set<double>("time step increase factor", ST_TS_INC_FACTOR);

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

        // only include a steady state time integrator list if not transient
        if (! ti_mode_list.isSublist("Transient")) {
          // create sublists for the steady state time integrator
          Teuchos::ParameterList& steady_time_integrator = flow_list->sublist("steady state time integrator");

          // error control options
          Teuchos::Array<std::string> err_opts;
          err_opts.push_back(std::string("pressure"));
          steady_time_integrator.set<Teuchos::Array<std::string> >("error control options",err_opts);

          // linear solver
          steady_time_integrator.set<std::string>("linear solver", ST_SOLVER);
          steady_time_integrator.set<std::string>("preconditioner", ST_PRECOND);

          // pressure-lambda constraints
          Teuchos::ParameterList &sti_plamb = steady_time_integrator.sublist("pressure-lambda constraints");
          sti_plamb.set<std::string>("method","projection");
          sti_plamb.set<std::string>("linear solver", ST_PLAMB_SOLVER);

          // time integration method
          steady_time_integrator.set<std::string>("time integration method","BDF1");
          Teuchos::ParameterList& sti_bdf1 = steady_time_integrator.sublist("BDF1");

          // use standard timestep controller type
          sti_bdf1.set<std::string>("timestep controller type", ST_TS_CONTROLLER);
          Teuchos::ParameterList &sti_bdf1_std = sti_bdf1.sublist("timestep controller standard parameters");
          sti_bdf1_std.set<int>("max iterations", ST_MAX_ITER);
          sti_bdf1_std.set<int>("min iterations", ST_MIN_ITER);
          sti_bdf1_std.set<double>("time step increase factor", ST_TS_INC_FACTOR);
          sti_bdf1_std.set<double>("time step reduction factor", ST_TS_RED_FACTOR);
          sti_bdf1_std.set<double>("max time step", ST_MAX_TS);
          sti_bdf1_std.set<double>("min time step", ST_MIN_TS);

          // solver type
          sti_bdf1.set<std::string>("solver type", "nka");
          Teuchos::ParameterList &sti_bdf1_nka = sti_bdf1.sublist("nka parameters");
          sti_bdf1_nka.set<double>("nonlinear tolerance", STEADY_NONLINEAR_TOLERANCE);
          sti_bdf1_nka.set<double>("diverged tolerance", ST_NKA_DIVGD_TOL);
          sti_bdf1_nka.set<double>("max du growth factor", ST_DIVERG_FACT);
          sti_bdf1_nka.set<int>("max divergent iterations", ST_MAX_DIVERGENT_ITERATIONS);
          sti_bdf1_nka.set<int>("max nka vectors", ST_NKA_NUMVEC);
	  sti_bdf1_nka.set<int>("limit iterations", ST_LIMIT_ITER);

          // remaining BDF1 parameters
          sti_bdf1.set<int>("max preconditioner lag iterations", ST_MAX_PREC_LAG);
          Teuchos::ParameterList &olist = steady_time_integrator.sublist("obsolete parameters");
          // olist.set<double>("nonlinear iteration damping factor", ST_NONLIN_DAMP);
          // olist.set<int>("nonlinear iteration initial guess extrapolation order",ST_NONLIN_INIT_GUESS_EXTR_ORD);
          // olist.set<double>("restart tolerance relaxation factor", ST_NONLIN_INIT_TS_FACTOR);
          // olist.set<double>("restart tolerance relaxation factor damping", ST_NONLIN_INIT_TS_FACTOR_DAMP);
          // olist.set<double>("error abs tol", ST_ERROR_ABS_TOL);
          // olist.set<double>("error rel tol", ST_ERROR_REL_TOL);
          // olist.set<std::string>("time stepping strategy", ST_TS_STRATEGY);

          if (plist->sublist("Execution Control").isSublist("Numerical Control Parameters")) {
            Teuchos::ParameterList& ncp_list =  plist->sublist("Execution Control").sublist("Numerical Control Parameters");
            if (ncp_list.isSublist("Unstructured Algorithm")) {
              Teuchos::ParameterList& ncpu_list = ncp_list.sublist("Unstructured Algorithm");
              if (ncpu_list.isSublist("Steady-State Implicit Time Integration")) {
                Teuchos::ParameterList& num_list = ncpu_list.sublist("Steady-State Implicit Time Integration");
                sti_bdf1_std.set<int>("max iterations", num_list.get<int>("steady max iterations", ST_MAX_ITER));
                sti_bdf1_std.set<int>("min iterations", num_list.get<int>("steady min iterations", ST_MIN_ITER));
                sti_bdf1_nka.set<int>("limit iterations", num_list.get<int>("steady limit iterations", ST_LIMIT_ITER));
                sti_bdf1_nka.set<double>("nonlinear tolerance",
                                     num_list.get<double>("steady nonlinear tolerance", STEADY_NONLINEAR_TOLERANCE));
                sti_bdf1_std.set<double>("time step reduction factor",
                                     num_list.get<double>("steady time step reduction factor", ST_TS_RED_FACTOR));
                sti_bdf1_std.set<double>("time step increase factor",
                                     num_list.get<double>("steady time step increase factor", ST_TS_INC_FACTOR));
                sti_bdf1_std.set<double>("max time step", num_list.get<double>("steady max time step", ST_MAX_TS));
                sti_bdf1.set<int>("max preconditioner lag iterations",
                                  num_list.get<int>("steady max preconditioner lag iterations", ST_MAX_PREC_LAG));
                // sti_bdf1.set<double>("error abs tol", num_list.get<double>("steady error abs tol", ST_ERROR_ABS_TOL));
                // sti_bdf1.set<double>("error rel tol", num_list.get<double>("steady error rel tol", ST_ERROR_REL_TOL));
                sti_bdf1_nka.set<int>("max divergent iterations",
                                  num_list.get<int>("steady max divergent iterations", ST_MAX_DIVERGENT_ITERATIONS));
                sti_bdf1.set<double>("nonlinear iteration damping factor",
                                     num_list.get<double>("steady nonlinear iteration damping factor", ST_NONLIN_DAMP));
                sti_bdf1.set<int>("nonlinear iteration initial guess extrapolation order",
                                  num_list.get<int>("steady nonlinear iteration initial guess extrapolation order", ST_NONLIN_INIT_GUESS_EXTR_ORD));
                sti_bdf1.set<double>("restart tolerance relaxation factor",
                                     num_list.get<double>("steady restart tolerance relaxation factor", ST_NONLIN_INIT_TS_FACTOR));
                sti_bdf1.set<double>("restart tolerance relaxation factor damping",
                                     num_list.get<double>("steady restart tolerance relaxation factor damping", ST_NONLIN_INIT_TS_FACTOR_DAMP));
                sti_bdf1_nka.set<double>("max du growth factor",
                                     num_list.get<double>("steady nonlinear iteration divergence factor", ST_DIVERG_FACT));

                steady_time_integrator.set<std::string>("preconditioner",
                                                        num_list.get<std::string>("steady preconditioner", ST_PRECOND));

                if (flow_model == "Single Phase") {
                  sti_bdf1_std.set<double>("time step increase factor",num_list.get<double>("steady time step increase factor",ST_SP_DT_INCR_FACTOR));
                }
		// initialization
		if (num_list.get<bool>("steady initialize with darcy", ST_INIT_DARCY_BOOL)) {
		  Teuchos::ParameterList &sti_init = steady_time_integrator.sublist("initialization");
		  sti_init.set<std::string>("method", "saturated solver");
		  sti_init.set<std::string>("linear solver", ST_INIT_SOLVER);
		}
              }
            }
          }
          if (nonlinear_solver == std::string("Newton")) {
            sti_bdf1.set<int>("max preconditioner lag iterations", 0);
          }
        }

        // only include the transient list if not in steady mode
        if ( ! ti_mode_list.isSublist("Steady")) {
          // create sublists for the transient state time integrator
          Teuchos::ParameterList& transient_time_integrator = flow_list->sublist("transient time integrator");

          // error control options
          Teuchos::Array<std::string> err_opts;
          err_opts.push_back(std::string("pressure"));
          transient_time_integrator.set<Teuchos::Array<std::string> >("error control options",err_opts);

          // linear solver
          transient_time_integrator.set<std::string>("linear solver", TR_SOLVER);
          transient_time_integrator.set<std::string>("preconditioner", TR_PRECOND);

          // pressure-lambda constraints
          Teuchos::ParameterList &tti_plamb = transient_time_integrator.sublist("pressure-lambda constraints");
          tti_plamb.set<std::string>("method","projection");
          tti_plamb.set<std::string>("linear solver", TR_PLAMB_SOLVER);

          // time integration method
          transient_time_integrator.set<std::string>("time integration method","BDF1");
          Teuchos::ParameterList& tti_bdf1 = transient_time_integrator.sublist("BDF1");

          // use standard timestep controller type
          tti_bdf1.set<std::string>("timestep controller type", TR_TS_CONTROLLER);
          Teuchos::ParameterList &tti_bdf1_std = tti_bdf1.sublist("timestep controller standard parameters");
          tti_bdf1_std.set<int>("max iterations", TR_MAX_ITER);
          tti_bdf1_std.set<int>("min iterations", TR_MIN_ITER);
          tti_bdf1_std.set<double>("time step increase factor", TR_TS_INC_FACTOR);
          tti_bdf1_std.set<double>("time step reduction factor", TR_TS_RED_FACTOR);
          tti_bdf1_std.set<double>("max time step", TR_MAX_TS);
          tti_bdf1_std.set<double>("min time step", TR_MIN_TS);

          // solver type
          tti_bdf1.set<std::string>("solver type", "nka");
          Teuchos::ParameterList &tti_bdf1_nka = tti_bdf1.sublist("nka parameters");
          tti_bdf1_nka.set<double>("nonlinear tolerance", TRANSIENT_NONLINEAR_TOLERANCE);
          tti_bdf1_nka.set<double>("diverged tolerance", TR_NKA_DIVGD_TOL);
          tti_bdf1_nka.set<double>("max du growth factor", TR_DIVERG_FACT);
          tti_bdf1_nka.set<int>("max divergent iterations", TR_MAX_DIVERGENT_ITERATIONS);
          tti_bdf1_nka.set<int>("max nka vectors", TR_NKA_NUMVEC);
	  tti_bdf1_nka.set<int>("limit iterations", TR_LIMIT_ITER);

          // remaining parameters
          tti_bdf1.set<int>("max preconditioner lag iterations", TR_MAX_PREC_LAG);
          Teuchos::ParameterList &olist = transient_time_integrator.sublist("obsolete parameters");
          // olist.set<int>("maximum number of iterations", TR_LIMIT_ITER); // this is not limit iters
          // olist.set<double>("nonlinear iteration damping factor", TR_NONLIN_DAMP);
          // olist.set<int>("nonlinear iteration initial guess extrapolation order", TR_NONLIN_INIT_GUESS_EXTR_ORD);
          // olist.set<double>("restart tolerance relaxation factor", TR_NONLIN_INIT_TS_FACTOR);
          // olist.set<double>("restart tolerance relaxation factor damping", TR_NONLIN_INIT_TS_FACTOR_DAMP);
          // olist.set<double>("error abs tol", TR_ERROR_ABS_TOL);
          // olist.set<double>("error rel tol", TR_ERROR_REL_TOL);
          // olist.set<std::string>("time stepping strategy", TR_TS_STRATEGY);

          if (plist->sublist("Execution Control").isSublist("Numerical Control Parameters")) {
            Teuchos::ParameterList& ncp_list = plist->sublist("Execution Control").sublist("Numerical Control Parameters");
            if (ncp_list.isSublist("Unstructured Algorithm")) {
              Teuchos::ParameterList& ncpu_list = ncp_list.sublist("Unstructured Algorithm");
              if (ncpu_list.isSublist("Transient Implicit Time Integration")) {

                Teuchos::ParameterList& num_list = ncpu_list.sublist("Transient Implicit Time Integration");

                tti_bdf1_std.set<int>("max iterations", num_list.get<int>("transient max iterations", TR_MAX_ITER));
                tti_bdf1_std.set<int>("min iterations", num_list.get<int>("transient min iterations", TR_MIN_ITER));
                tti_bdf1_nka.set<int>("limit iterations", num_list.get<int>("transient limit iterations", TR_LIMIT_ITER));
                tti_bdf1_nka.set<double>("nonlinear tolerance",
                                     num_list.get<double>("transient nonlinear tolerance", TRANSIENT_NONLINEAR_TOLERANCE));
                tti_bdf1_std.set<double>("time step reduction factor",
                                     num_list.get<double>("transient time step reduction factor", TR_TS_RED_FACTOR));
                tti_bdf1_std.set<double>("time step increase factor",
                                     num_list.get<double>("transient time step increase factor", TR_TS_INC_FACTOR));
                tti_bdf1_std.set<double>("max time step", num_list.get<double>("transient max time step", TR_MAX_TS));
                tti_bdf1.set<int>("max preconditioner lag iterations",
                                  num_list.get<int>("transient max preconditioner lag iterations", TR_MAX_PREC_LAG));
                // tti_bdf1.set<double>("error abs tol", num_list.get<double>("transient error abs tol", TR_ERROR_ABS_TOL));
                // tti_bdf1.set<double>("error rel tol", num_list.get<double>("transient error rel tol", TR_ERROR_REL_TOL));
                tti_bdf1_nka.set<int>("max divergent iterations",
                                  num_list.get<int>("transient max divergent iterations", TR_MAX_DIVERGENT_ITERATIONS));
                tti_bdf1.set<double>("nonlinear iteration damping factor",
                                     num_list.get<double>("transient nonlinear iteration damping factor", TR_NONLIN_DAMP));
                tti_bdf1.set<int>("nonlinear iteration initial guess extrapolation order",
                                  num_list.get<int>("transient nonlinear iteration initial guess extrapolation order", TR_NONLIN_INIT_GUESS_EXTR_ORD));
                tti_bdf1.set<double>("restart tolerance relaxation factor",
                                     num_list.get<double>("transient restart tolerance relaxation factor", TR_NONLIN_INIT_TS_FACTOR));
                tti_bdf1.set<double>("restart tolerance relaxation factor damping",
                                     num_list.get<double>("transient restart tolerance relaxation factor damping", TR_NONLIN_INIT_TS_FACTOR_DAMP));
                tti_bdf1_nka.set<double>("max du growth factor",
                                     num_list.get<double>("transient nonlinear iteration divergence factor", TR_DIVERG_FACT));

                transient_time_integrator.set<std::string>("preconditioner",
                                                           num_list.get<std::string>("transient preconditioner", TR_PRECOND));

                if (flow_model == "Single Phase") {
                  tti_bdf1_std.set<double>("time step increase factor",
                                           num_list.get<double>("transient time step increase factor", TR_SP_DT_INCR_FACTOR));
                }
		// create an initialization sublist
		if (num_list.get<bool>("transient initialize with darcy", TR_INIT_DARCY_BOOL)) {
		  Teuchos::ParameterList &tti_init = transient_time_integrator.sublist("initialization");
		  tti_init.set<std::string>("method","saturated solver");
		  tti_init.set<std::string>("linear solver", TR_INIT_SOLVER);
		} else {
		  Teuchos::ParameterList &tti_init = transient_time_integrator.sublist("initialization");
		  tti_init.set<std::string>("method","projection");
		  tti_init.set<std::string>("linear solver", TR_INIT_SOLVER);
		}
	      }
            }
          }

          if (nonlinear_solver == std::string("Newton")) {
            tti_bdf1.set<int>("max preconditioner lag iterations", 0);
          }


          transient_time_integrator.sublist("VerboseObject") = create_Verbosity_List(verbosity_level);
        }
      }
    }
  }

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
        src_sub_out.set<std::string>("spatial distribution method","volume");
        src_fn = src.sublist("Source: Volume Weighted");
      } else if (src.isSublist("Source: Permeability Weighted")) {
        src_sub_out.set<std::string>("spatial distribution method","permeability");
        src_fn = src.sublist("Source: Permeability Weighted");
      } else if (src.isSublist("Source: Uniform")) {
        src_sub_out.set<std::string>("spatial distribution method","none");
        src_fn = src.sublist("Source: Uniform");
      } else {
        Exceptions::amanzi_throw(Errors::Message("In the definition of Sources: you must either specify 'Source: Volume Weighted', 'Source: Permeability Weighted', or 'Source: Uniform'"));
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

  return src_list;
}





Teuchos::ParameterList create_TransportSrc_List(Teuchos::ParameterList* plist)
{
  Teuchos::ParameterList src_list;

  Teuchos::ParameterList& src_sublist = plist->sublist("Sources");

  for (Teuchos::ParameterList::ConstIterator i = src_sublist.begin(); i != src_sublist.end(); i++) {
    // look at sublists
    if (src_sublist.isSublist(src_sublist.name(i))) {
      Teuchos::ParameterList& src = src_sublist.sublist(src_sublist.name(i));
      // get name
      std::string name = src_sublist.name(i);

      // get the regions
      Teuchos::Array<std::string> regions = src.get<Teuchos::Array<std::string> >("Assigned Regions");

      std::string dist_method("none");
      if (src.isSublist("Source: Volume Weighted")) {
        dist_method = "volume";
      } else if (src.isSublist("Source: Permeability Weighted")) {
        dist_method = "permeability";
      }


      // go to the phase list
      if (src.isSublist("Solute SOURCE")) {
        if (src.sublist("Solute SOURCE").isSublist("Aqueous")) {
          if (src.sublist("Solute SOURCE").sublist("Aqueous").isSublist(phase_comp_name)) {

            Teuchos::ParameterList src_bc_list = src.sublist("Solute SOURCE").sublist("Aqueous").sublist(phase_comp_name);

            // loop over all the source definitions
            for (Teuchos::ParameterList::ConstIterator ibc = src_bc_list.begin(); ibc != src_bc_list.end(); ibc++) {
              Teuchos::ParameterList& solute_src = src_bc_list.sublist(src_bc_list.name(ibc));
              std::string solute_name = src_bc_list.name(ibc);

              // create src sublist
              Teuchos::ParameterList& src_sub_out = src_list.sublist("concentration").sublist(solute_name).sublist("source for " + name);
              src_sub_out.set<Teuchos::Array<std::string> >("regions",regions);

              // get source function
              Teuchos::ParameterList src_fn;
              if (solute_src.isSublist("Source: Uniform Concentration")) {
                src_sub_out.set<std::string>("spatial distribution method","none");
                src_fn = solute_src.sublist("Source: Uniform Concentration");
              } else if (solute_src.isSublist("Source: Flow Weighted Concentration")) {
                src_sub_out.set<std::string>("spatial distribution method",dist_method);
                src_fn = solute_src.sublist("Source: Flow Weighted Concentration");
              } else {
                Exceptions::amanzi_throw(Errors::Message("In the definition of Sources: you must either specify 'Source: Uniform Concentration' or 'Source: Flow Weighted Concentration'."));
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
        }
      }
    }
  }

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

        wrm_sublist.set<std::string>("water retention model", "van Genuchten");
        wrm_sublist.set<std::string>("region", *i);
        wrm_sublist.set<double>("van Genuchten m", m);
        wrm_sublist.set<double>("van Genuchten l", ell);
        wrm_sublist.set<double>("van Genuchten alpha", alpha);
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

      double lambda = BC_list.get<double>("lambda");
      double alpha  = BC_list.get<double>("alpha");
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

        wrm_sublist.set<std::string>("water retention model", "Brooks Corey");
        wrm_sublist.set<std::string>("region", *i);
        wrm_sublist.set<double>("Brooks Corey lambda", lambda);
        wrm_sublist.set<double>("Brooks Corey alpha", alpha);
        wrm_sublist.set<double>("Brooks Corey l", ell);
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

  dpc_list.set<std::string>("discretization method", "optimized mfd scaled");
  dpc_list.set<std::string>("preconditioner type", "ml");

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

  Teuchos::ParameterList& ml_list = dpc_list.sublist("ml parameters");
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

  bilu_list.set<std::string>("discretization method", "optimized mfd scaled");
  bilu_list.set<std::string>("preconditioner type", "block ilu");

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

  Teuchos::ParameterList& p_list = bilu_list.sublist("block ilu parameters");
  p_list.set<double>("fact: relax value",bilu_relax_value);
  p_list.set<double>("fact: absolute threshold",bilu_abs_thresh);
  p_list.set<double>("fact: relative threshold",bilu_rel_thresh);
  p_list.set<int>("fact: level-of-fill",bilu_level_of_fill);
  p_list.set<int>("overlap",bilu_overlap);
  p_list.set<std::string>("schwarz: combine mode","Add");

  return bilu_list;
}



/* ******************************************************************
 * HypreBoomerAMG preconditioner sublist
 ****************************************************************** */
Teuchos::ParameterList create_HypreAMG_List(Teuchos::ParameterList* plist)
{
  Teuchos::ParameterList dpc_list;

  dpc_list.set<std::string>("discretization method", "optimized mfd scaled");
  dpc_list.set<std::string>("preconditioner type", "boomer amg");

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

  Teuchos::ParameterList& amg_list = dpc_list.sublist("boomer amg parameters");
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

        tbc.set<bool>("rainfall", bc_flux.get<bool>("rainfall",false));

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
        std::string                 coordsys = bc_dir.get<std::string>("Coordinate System",BCHYDRST_COORD);

        std::stringstream ss;
        ss << "BC " << bc_counter++;

        Teuchos::ParameterList& tbc = ssf_list.sublist("static head").sublist(ss.str());
        tbc.set<Teuchos::Array<std::string> >("regions", regions );

        if (coordsys == "Absolute") {
          tbc.set<bool>("relative to top",false);
        } else if (coordsys == "Relative") {
          tbc.set<bool>("relative to top",true);
        } else {
          // we have a default for this value... "Absolute", if for some reason this does not
          // get read, then we must bail
          Exceptions::amanzi_throw(Errors::Message("In 'BC: Hydrostatic': must specify a value for parameter 'Coordinate System', valid values are 'Absolute' and 'Relative'"));
        }


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



Teuchos::ParameterList create_State_List(Teuchos::ParameterList* plist) {
  Teuchos::ParameterList stt_list;

  // first we write initial conditions for scalars and vectors, not region-specific

  Teuchos::ParameterList & stt_ic = stt_list.sublist("initial conditions");
  //
  // --- gravity
  Teuchos::Array<double> gravity(spatial_dimension_);
  for (int i=0; i!=spatial_dimension_-1; ++i) gravity[i] = 0.0;
  gravity[spatial_dimension_-1] =  - GRAVITY_MAGNITUDE;
  stt_ic.sublist("gravity").set<Teuchos::Array<double> >("value", gravity);
  //
  // --- viscosity
  //
  Teuchos::ParameterList& phase_list = plist->sublist("Phase Definitions");
  double viscosity = phase_list.sublist(phase_name).sublist("Phase Properties").sublist("Viscosity: Uniform").get<double>("Viscosity");
  stt_ic.sublist("fluid_viscosity").set<double>("value", viscosity);
  //
  // --- density
  //
  double density = phase_list.sublist(phase_name).sublist("Phase Properties").sublist("Density: Uniform").get<double>("Density");
  stt_ic.sublist("fluid_density").set<double>("value", density);
  // this is stupid, but for some reason we also have an array for water density, so here it goes...
  stt_ic.sublist("water_density").sublist("function").sublist("All")
      .set<std::string>("region","All")
      .set<std::string>("component","cell")
      .sublist("function").sublist("function-constant")
      .set<double>("value", density);
  //
  // --- region specific initial conditions from material properties
  //
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

    // read the parameters that need to be written for each region in the Assigned regions list

    // create regions string
    std::string reg_str;
    for (Teuchos::Array<std::string>::const_iterator ireg=regions.begin(); ireg!=regions.end(); ++ireg) {
      reg_str = reg_str + *ireg;
    }

    // porosity...
    double porosity;
    if (matprop_list.sublist(matprop_list.name(i)).isSublist("Porosity: Uniform")) {
      porosity = matprop_list.sublist(matprop_list.name(i)).sublist("Porosity: Uniform").get<double>("Value");
    } else {
      Exceptions::amanzi_throw(Errors::Message("Porosity must be specified as Intrinsic Porosity: Uniform, for every region."));
    }
    Teuchos::ParameterList &porosity_ic = stt_ic.sublist("porosity");
    porosity_ic.sublist("function").sublist(reg_str)
        .set<Teuchos::Array<std::string> >("regions",regions)
        .set<std::string>("component","cell")
        .sublist("function").sublist("function-constant")
        .set<double>("value", porosity);

    // permeability...
    double perm_x, perm_y, perm_z;
    std::string perm_file, perm_attribute, perm_format;
    bool perm_init_from_file = false;
    Teuchos::ParameterList& mplist = matprop_list.sublist(matprop_list.name(i));
    if ((mplist.isSublist("Intrinsic Permeability: Uniform") || 
	 mplist.isSublist("Intrinsic Permeability: Anisotropic Uniform")) && 
	(mplist.isSublist("Hydraulic Conductivity: Uniform") || 
	 mplist.isSublist("Hydraulic Conductivity: Anisotropic Uniform"))) {
      Exceptions::amanzi_throw(Errors::Message("Permeability can only be specified either Intrinsic Permeability or Hydraulic Conductivity, but not both."));
    }

    if (mplist.isSublist("Intrinsic Permeability: Uniform")) {
      perm_x = perm_y = perm_z = mplist.sublist("Intrinsic Permeability: Uniform").get<double>("Value");
    } else if (mplist.isSublist("Intrinsic Permeability: Anisotropic Uniform")) {
      perm_x = mplist.sublist("Intrinsic Permeability: Anisotropic Uniform").get<double>("x");
      perm_y = mplist.sublist("Intrinsic Permeability: Anisotropic Uniform").get<double>("y");
      if (mplist.sublist("Intrinsic Permeability: Anisotropic Uniform").isParameter("z") ) {
        if (spatial_dimension_ == 3) {
          perm_z = mplist.sublist("Intrinsic Permeability: Anisotropic Uniform").get<double>("z");
        } else {
          Exceptions::amanzi_throw(Errors::Message("Intrinsic Permeability: Anisotropic Uniform defines a value for z, while the spatial dimension of the problem not 3."));
        }
      }
    } else if (mplist.isSublist("Hydraulic Conductivity: Uniform")) {
      perm_x = perm_y = perm_z = mplist.sublist("Hydraulic Conductivity: Uniform").get<double>("Value");
      // now scale with rho, g and mu to get correct permeability values
      perm_x *= viscosity/(density*GRAVITY_MAGNITUDE);
      perm_y *= viscosity/(density*GRAVITY_MAGNITUDE);
      perm_z *= viscosity/(density*GRAVITY_MAGNITUDE);
    } else if (mplist.isSublist("Hydraulic Conductivity: Anisotropic Uniform")) {
      perm_x = mplist.sublist("Hydraulic Conductivity: Anisotropic Uniform").get<double>("x");
      perm_y = mplist.sublist("Hydraulic Conductivity: Anisotropic Uniform").get<double>("y");
      if (mplist.sublist("Hydraulic Conductivity: Anisotropic Uniform").isParameter("z")) {
        if (spatial_dimension_ == 3) {
          perm_z = mplist.sublist("Hydraulic Conductivity: Anisotropic Uniform").get<double>("z");
        } else {
          Exceptions::amanzi_throw(Errors::Message("Hydraulic Conductivity: Anisotropic Uniform defines a value for z, while the spatial dimension of the problem not 3."));
        }
      }
      // now scale with rho, g and mu to get correct permeablity values
      perm_x *= viscosity/(density*GRAVITY_MAGNITUDE);
      perm_y *= viscosity/(density*GRAVITY_MAGNITUDE);
      perm_z *= viscosity/(density*GRAVITY_MAGNITUDE);
    } else if (mplist.isSublist("Intrinsic Permeability: File")) {
      Teuchos::ParameterList &aux_list = mplist.sublist("Intrinsic Permeability: File");
      bool input_error = false;
      if (aux_list.isParameter("File")) {
        perm_file = aux_list.get<std::string>("File");
      } else {
        input_error = true;
      }
      if (aux_list.isParameter("Attribute")) {
        perm_attribute = aux_list.get<std::string>("Attribute");
      } else {
        input_error = true;
      }
      if (aux_list.isParameter("Format")) {
        perm_format = mplist.sublist("Intrinsic Permeability: File").get<std::string>("Format");
      } else {
        input_error = true;
      }
      if (input_error) {
        Exceptions::amanzi_throw(Errors::Message("The list 'Input Permeability: File' could not be parsed, a required parameter is missing. Check the input specification."));
      }
      perm_init_from_file = true;
    } else {
      Exceptions::amanzi_throw(Errors::Message("Permeability can only be specified as 'Intrinsic Permeability: Uniform', 'Intrinsic Permeability: Anisotropic Uniform', 'Hydraulic Conductivity: Uniform', 'Hydraulic Conductivity: Anisotropic Uniform', or 'Intrinsic Permeability: File'"));
    }

    Teuchos::ParameterList &permeability_ic = stt_ic.sublist("permeability");
    if (perm_init_from_file) {
      if (perm_format == std::string("exodus")) {
        // first make sure the file actually exists

        boost::filesystem::path p;
        if (numproc_>1) {
          // attach the right extensions as required by Nemesis file naming conventions
          // in which files are named as mymesh.par.N.r where N = numproc and r is rank
          // and check if files exist
          std::string perm_file_par = perm_file.substr(0,perm_file.size()-4) + std::string(".par");
          int rank = numproc_-1;
          int ndigits = static_cast<int>(floor(log10(numproc_))) + 1;
          std::string fmt = boost::str(boost::format("%%s.%%d.%%0%dd") % ndigits);
          std::string par_file_w_ext = boost::str(boost::format(fmt) % perm_file_par % numproc_ % rank);
          p = par_file_w_ext;
          if (!boost::filesystem::exists(p)) {
            Exceptions::amanzi_throw(Errors::Message("Permeability initialization from file: not all the partitioned files '" + perm_file_par + ".<numranks>.<rank>' exist, and possibly none of them exist."));
          }
          perm_file = perm_file_par;
        } else { // numproc_ == 1
          std::string suffix = perm_file.substr(perm_file.size()-4);
          if (suffix != std::string(".exo")) {          
            Exceptions::amanzi_throw(Errors::Message("Permeability initialization from file: in serial the exodus mesh file must have the suffix .exo"));
          }
          p = perm_file;
          if (!boost::filesystem::exists(p)) {
            Exceptions::amanzi_throw(Errors::Message("Permeability initialization from file: the file '" + perm_file + "' does not exist."));
          }
        }
        permeability_ic.sublist("exodus file initialization")
            .set<std::string>("file",perm_file)
            .set<std::string>("attribute",perm_attribute);
      } else {
        Exceptions::amanzi_throw(Errors::Message("Permeabily initialization from file, incompatible format specified: '" + perm_format + "', only 'exodus' is supported."));
      }
     
    } else {
      Teuchos::ParameterList& aux_list =
        permeability_ic.sublist("function").sublist(reg_str)
        .set<Teuchos::Array<std::string> >("regions",regions)
        .set<std::string>("component","cell")
        .sublist("function");
      aux_list.set<int>("Number of DoFs",spatial_dimension_)
        .set<std::string>("Function type","composite function");
      aux_list.sublist("DoF 1 Function").sublist("function-constant")
        .set<double>("value", perm_x);
      if (spatial_dimension_ >= 2) {
        aux_list.sublist("DoF 2 Function").sublist("function-constant")
          .set<double>("value", perm_y);
      }
      if (spatial_dimension_ == 3) {
        aux_list.sublist("DoF 3 Function").sublist("function-constant")
          .set<double>("value", perm_z);
      }
    }

    // specific_yield...
    if (matprop_list.sublist(matprop_list.name(i)).isSublist("Specific Yield: Uniform")) {
      Teuchos::ParameterList& spec_yield_ic = stt_ic.sublist("specific_yield");
      double specific_yield = matprop_list.sublist(matprop_list.name(i)).sublist("Specific Yield: Uniform").get<double>("Value");
      spec_yield_ic.sublist("function").sublist(reg_str)
          .set<Teuchos::Array<std::string> >("regions",regions)
          .set<std::string>("component","cell")
          .sublist("function").sublist("function-constant")
          .set<double>("value", specific_yield);
    }

    // specific_storage...
    if (matprop_list.sublist(matprop_list.name(i)).isSublist("Specific Storage: Uniform")) {
      Teuchos::ParameterList& spec_stor_ic = stt_ic.sublist("specific_storage");
      double specific_storage = matprop_list.sublist(matprop_list.name(i)).sublist("Specific Storage: Uniform").get<double>("Value");
      spec_stor_ic.sublist("function").sublist(reg_str)
          .set<Teuchos::Array<std::string> >("regions",regions)
          .set<std::string>("component","cell")
          .sublist("function").sublist("function-constant")
          .set<double>("value", specific_storage);
    }

    // particle_density...
    if (matprop_list.sublist(matprop_list.name(i)).isSublist("Particle Density: Uniform")) {
      Teuchos::ParameterList& part_dens_ic = stt_ic.sublist("particle_density");
      double particle_density = matprop_list.sublist(matprop_list.name(i)).sublist("Particle Density: Uniform").get<double>("Value");
      part_dens_ic.sublist("function").sublist(reg_str)
          .set<Teuchos::Array<std::string> >("regions",regions)
          .set<std::string>("component","cell")
          .sublist("function").sublist("function-constant")
          .set<double>("value", particle_density);

    }
  }

  // and now the initialization of fields that are initialized via an Initial Condition in the Akuna spec

  // loop over the initial conditions
  Teuchos::ParameterList& ic_list = plist->sublist("Initial Conditions");

  if (! ic_list.isParameter("Init from Checkpoint File")) {
    // only process initial conditions if we are not initializing from
    // a checkpoint file
    for (Teuchos::ParameterList::ConstIterator iic = ic_list.begin(); iic != ic_list.end(); ++iic) {
      // get the regions
      Teuchos::Array<std::string> regions = ic_list.sublist(ic_list.name(iic)).get<Teuchos::Array<std::string> >("Assigned Regions");

      Teuchos::ParameterList* ic_for_region = &ic_list.sublist(ic_list.name(iic));

      // create regions string
      std::string reg_str;
      for (Teuchos::Array<std::string>::const_iterator ireg=regions.begin(); ireg!=regions.end(); ++ireg) {
        reg_str = reg_str + *ireg;
      }

      // pressure...
      if (ic_for_region->isSublist("IC: Uniform Pressure") || ic_for_region->isSublist("IC: Linear Pressure")) {

        Teuchos::ParameterList &pressure_ic = stt_ic.sublist("pressure");
        if ( ic_for_region->isSublist("IC: Uniform Pressure")) {
          double p = ic_for_region->sublist("IC: Uniform Pressure").get<double>("Value");

          pressure_ic.sublist("function").sublist(reg_str)
              .set<Teuchos::Array<std::string> >("regions",regions)
              .set<std::string>("component","cell")
              .sublist("function").sublist("function-constant")
              .set<double>("value", p);

        } else if (ic_for_region->isSublist("IC: Linear Pressure")) {
          Teuchos::Array<double> grad = ic_for_region->sublist("IC: Linear Pressure").get<Teuchos::Array<double> >("Gradient Value");
          Teuchos::Array<double> refcoord = ic_for_region->sublist("IC: Linear Pressure").get<Teuchos::Array<double> >("Reference Coordinate");
          double refval =  ic_for_region->sublist("IC: Linear Pressure").get<double>("Reference Value");

          Teuchos::Array<double> grad_with_time(grad.size()+1);
          grad_with_time[0] = 0.0;
          for (int j=0; j!=grad.size(); ++j) {
            grad_with_time[j+1] = grad[j];
          }

          Teuchos::Array<double> refcoord_with_time(refcoord.size()+1);
          refcoord_with_time[0] = 0.0;
          for (int j=0; j!=refcoord.size(); ++j) {
            refcoord_with_time[j+1] = refcoord[j];
          }

          pressure_ic.sublist("function").sublist(reg_str)
              .set<Teuchos::Array<std::string> >("regions",regions)
              .set<std::string>("component","cell")
              .sublist("function").sublist("function-linear")
              .set<double>("y0", refval)
              .set<Teuchos::Array<double> >("x0",refcoord_with_time)
              .set<Teuchos::Array<double> >("gradient",grad_with_time);

        } else if (ic_for_region->isSublist("IC: File Pressure")) {
          Exceptions::amanzi_throw(Errors::Message("IC: File Pressure cannot currently be used to initialize pressure in a region."));
        }
      }


      // saturation...
      if (ic_for_region->isSublist("IC: Uniform Saturation") || ic_for_region->isSublist("IC: Linear Saturation")) {

        Teuchos::ParameterList &saturation_ic = stt_ic.sublist("water_saturation");
        if ( ic_for_region->isSublist("IC: Uniform Saturation")) {
          double s = ic_for_region->sublist("IC: Uniform Saturation").get<double>("Value");
          saturation_ic.sublist("function").sublist(reg_str)
              .set<Teuchos::Array<std::string> >("regions",regions)
              .set<std::string>("component","cell")
              .sublist("function").sublist("function-constant")
              .set<double>("value", s);
        } else if ( ic_for_region->isSublist("IC: Linear Saturation")) {
          Teuchos::Array<double> grad = ic_for_region->sublist("IC: Linear Saturation").get<Teuchos::Array<double> >("Gradient Value");
          Teuchos::Array<double> refcoord = ic_for_region->sublist("IC: Linear Saturation").get<Teuchos::Array<double> >("Reference Coordinate");
          double refval =  ic_for_region->sublist("IC: Linear Saturation").get<double>("Reference Value");

          Teuchos::Array<double> grad_with_time(grad.size()+1);
          grad_with_time[0] = 0.0;
          for (int j=0; j!=grad.size(); ++j) {
            grad_with_time[j+1] = grad[j];
          }

          Teuchos::Array<double> refcoord_with_time(refcoord.size()+1);
          refcoord_with_time[0] = 0.0;
          for (int j=0; j!=refcoord.size(); ++j) {
            refcoord_with_time[j+1] = refcoord[j];
          }

          saturation_ic.sublist("function").sublist(reg_str)
              .set<Teuchos::Array<std::string> >("regions",regions)
              .set<std::string>("component","cell")
              .sublist("function").sublist("function-linear")
              .set<double>("y0", refval)
              .set<Teuchos::Array<double> >("x0",refcoord_with_time)
              .set<Teuchos::Array<double> >("gradient",grad_with_time);
        } else if (ic_for_region->isSublist("IC: File Saturation")) {
          Exceptions::amanzi_throw(Errors::Message("IC: File Saturation cannot currently be used to initialize saturation in a region."));
        }
      }

      // darcy_flux...
      if (ic_for_region->isSublist("IC: Uniform Velocity")) {
        Teuchos::ParameterList &darcy_flux_ic =  stt_ic.sublist("darcy_flux");
        Teuchos::Array<double> vel_vec = ic_for_region->sublist("IC: Uniform Velocity").get<Teuchos::Array<double> >("Velocity Vector");

        if (vel_vec.size() != spatial_dimension_) {
          Exceptions::amanzi_throw(Errors::Message("The velocity vector defined in in the IC: Uniform Velocity list does not match the spatial dimension of the problem."));
        }

        Teuchos::ParameterList& d_list =
            darcy_flux_ic.set<bool>("dot with normal", true)
            .sublist("function").sublist(reg_str)
            .set<Teuchos::Array<std::string> >("regions",regions)
            .set<std::string>("component","face")
            .sublist("function")
            .set<int>("Number of DoFs", vel_vec.size())
            .set<std::string>("Function type", "composite function");

        for (int ii=0; ii != vel_vec.size(); ++ii) {
          std::stringstream dof_str;
          dof_str << "DoF " << ii+1 << " Function";

          d_list.sublist(dof_str.str())
              .sublist("function-constant")
              .set<double>("value", vel_vec[ii]);
        }
      }


      // total_component_concentration...
      if (ic_for_region->isSublist("Solute IC")) {

        Teuchos::ParameterList &concentration_ic = stt_ic.sublist("total_component_concentration");
        if (plist->sublist("Execution Control").get<std::string>("Transport Model") != std::string("Off")  ||
            plist->sublist("Execution Control").get<std::string>("Chemistry Model") != std::string("Off") ) {
          // write the initial conditions for the solutes, note that we hardcode for there only being one phase, with one phase component

          Teuchos::ParameterList & dof_list = concentration_ic.sublist("function").sublist(reg_str)
              .set<Teuchos::Array<std::string> >("regions",regions)
              .set<std::string>("component","cell")
              .sublist("function")
              .set<int>("Number of DoFs", comp_names.size())
              .set<std::string>("Function type", "composite function");


          for (int ii=0; ii<comp_names.size(); ii++) {
            if (ic_for_region->sublist("Solute IC").sublist(phase_name).sublist(phase_comp_name).isSublist(comp_names[ii])) {

              double conc = ic_for_region->sublist("Solute IC").sublist(phase_name).sublist(phase_comp_name).sublist(comp_names[ii]).sublist("IC: Uniform Concentration").get<double>("Value");

              std::stringstream dof_str;
              dof_str << "DoF " << ii+1 << " Function";

              dof_list.sublist(dof_str.str())
                  .sublist("function-constant")
                  .set<double>("value",conc);
            }
          }
        }
      }
    }
  }

  return stt_list;
}



/* ******************************************************************
 * populates parameters in the State list.
 ****************************************************************** */
Teuchos::ParameterList create_State_List_old(Teuchos::ParameterList* plist) {
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
    if (comp_names.size() > 0) {
      stt_list.set<Teuchos::Array<std::string> >("Component Solutes", comp_names);
    }
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
  } else if (vlevel == "extreme") {
    vlist.set<std::string>("Verbosity Level","extreme");
  } else if (vlevel == "none") {
    vlist.set<std::string>("Verbosity Level","none");
  }

  return vlist;
}


/* ******************************************************************
 * Empty
 ****************************************************************** */
Teuchos::ParameterList CreateChemistryList(Teuchos::ParameterList* plist) {
  Teuchos::ParameterList chem_list;
  if ( plist->isSublist("Chemistry") ) {
    chem_list = plist->sublist("Chemistry");

    Teuchos::ParameterList & chem_ic = chem_list.sublist("initial conditions");

    //
    // read the minerals
    //

    // Teuchos::Array<std::string> minerals;
    // if (plist->sublist("Phase Definitions").isSublist("Solid")) {
    //minerals = plist->sublist("Phase Definitions").sublist("Solid")
    // .get<Teuchos::Array<std::string> >("Minerals");
    if (mineral_names_.size() > 0) {
      chem_list.set<Teuchos::Array<std::string> >("Minerals", mineral_names_);
    }
    if (sorption_site_names_.size() > 0) {
      chem_list.set<Teuchos::Array<std::string> >("Sorption Sites", sorption_site_names_);
    }

    chem_list.set<int>("Number of component concentrations", comp_names.size() );

    //
    // --- region specific initial conditions
    //
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

      if (mineral_names_.size() > 0) {
        double mvf(0.0), msa(0.0);

        // mineral volume fractions
        Teuchos::ParameterList &mineral_volfrac_ic = chem_ic.sublist("mineral_volume_fractions");
        // mineral specific surface area
        Teuchos::ParameterList &mineral_surfarea_ic = chem_ic.sublist("mineral_specific_surface_area");


        for (Teuchos::Array<std::string>::const_iterator ir=regions.begin(); ir!=regions.end(); ir++) {

          // mineral volume fractions and specific surface area

          Teuchos::ParameterList & aux1_list =
              mineral_volfrac_ic.sublist("function").sublist(*ir)
              .set<std::string>("region",*ir)
              .set<std::string>("component","cell")
              .sublist("function");

          Teuchos::ParameterList & aux2_list =
              mineral_surfarea_ic.sublist("function").sublist(*ir)
              .set<std::string>("region",*ir)
              .set<std::string>("component","cell")
              .sublist("function");

          aux1_list.set<int>("Number of DoFs", mineral_names_.size())
              .set("Function type", "composite function");

          aux2_list.set<int>("Number of DoFs", mineral_names_.size())
              .set("Function type", "composite function");

          for (int j = 0; j<mineral_names_.size(); ++j) {
            std::stringstream ss;
            ss << "DoF " << j+1 << " Function";

            mvf = matprop_list.sublist(matprop_list.name(i)).sublist("Mineralogy").sublist(mineral_names_[j])
                .get<double>("Volume Fraction");

            msa = matprop_list.sublist(matprop_list.name(i)).sublist("Mineralogy").sublist(mineral_names_[j])
                .get<double>("Specific Surface Area");

            aux1_list.sublist(ss.str()).sublist("function-constant")
                .set<double>("value", mvf);

            aux2_list.sublist(ss.str()).sublist("function-constant")
                .set<double>("value", msa);
          }
        }
      }




      if ( matprop_list.sublist(matprop_list.name(i)).isParameter("Cation Exchange Capacity")) {

        double cec = matprop_list.sublist(matprop_list.name(i)).get<double>("Cation Exchange Capacity");

        Teuchos::ParameterList &ion_exchange_sites_ic = chem_ic.sublist("ion_exchange_sites");
        for (Teuchos::Array<std::string>::const_iterator ir=regions.begin(); ir!=regions.end(); ir++) {
          Teuchos::ParameterList & aux1_list =
              ion_exchange_sites_ic.sublist("function").sublist(*ir)
              .set<std::string>("region",*ir)
              .set<std::string>("component","cell")
              .sublist("function");

          // this needs to be more general... for now we initialize with one DoF
          aux1_list.set<int>("Number of DoFs", 1)
              .set("Function type", "composite function");

          aux1_list.sublist("DoF 1 Function").sublist("function-constant")
              .set<double>("value", cec);
        }


        Teuchos::ParameterList &ion_exchange_ref_cation_conc_ic = chem_ic.sublist("ion_exchange_ref_cation_conc");
        for (Teuchos::Array<std::string>::const_iterator ir=regions.begin(); ir!=regions.end(); ir++) {
          Teuchos::ParameterList & aux1_list =
              ion_exchange_ref_cation_conc_ic.sublist("function").sublist(*ir)
              .set<std::string>("region",*ir)
              .set<std::string>("component","cell")
              .sublist("function");

          // this needs to be more general... for now we initialize with one DoF
          aux1_list.set<int>("Number of DoFs", 1)
              .set("Function type", "composite function");

          aux1_list.sublist("DoF 1 Function").sublist("function-constant")
              .set<double>("value", 1.0);  // this should be read from the input file??... TODO
        }
      }

      if ( matprop_list.sublist(matprop_list.name(i)).isSublist("Sorption Isotherms")) {

        Teuchos::ParameterList & isotherm_kd_ic = chem_ic.sublist("isotherm_kd");
        Teuchos::ParameterList & isotherm_langmuir_b_ic = chem_ic.sublist("isotherm_langmuir_b");
        Teuchos::ParameterList & isotherm_freundlich_n_ic = chem_ic.sublist("isotherm_freundlich_n");

        Teuchos::ParameterList & sorption_isotherms_list = matprop_list.sublist(matprop_list.name(i)).sublist("Sorption Isotherms");

        for (Teuchos::Array<std::string>::const_iterator ir=regions.begin(); ir!=regions.end(); ir++) {

          // Kd

          {
            Teuchos::ParameterList & aux1_list =
                isotherm_kd_ic.sublist("function").sublist(*ir)
                .set<std::string>("region",*ir)
                .set<std::string>("component","cell")
                .sublist("function");

            aux1_list.set<int>("Number of DoFs", comp_names.size())
                .set("Function type", "composite function");

            for ( int ic = 0; ic != comp_names.size(); ++ic) {

              std::stringstream ss;
              ss << "DoF " << ic + 1 << " Function";

              double kd(0.0);
              if (sorption_isotherms_list.isSublist(comp_names[ic])) {
                kd = sorption_isotherms_list.sublist(comp_names[ic]).get<double>("Kd",0.0);
              }

              aux1_list.sublist(ss.str()).sublist("function-constant")
                  .set<double>("value", kd);

            }
          }

          // Langmuir

          {
            Teuchos::ParameterList & aux2_list =
                isotherm_langmuir_b_ic.sublist("function").sublist(*ir)
                .set<std::string>("region",*ir)
                .set<std::string>("component","cell")
                .sublist("function");

            aux2_list.set<int>("Number of DoFs", comp_names.size())
                .set("Function type", "composite function");

            for ( int ic = 0; ic != comp_names.size(); ++ic) {

              std::stringstream ss;
              ss << "DoF " << ic + 1 << " Function";

              double langmuir_b(1.0);
              if (sorption_isotherms_list.isSublist(comp_names[ic])) {
                langmuir_b = sorption_isotherms_list.sublist(comp_names[ic]).get<double>("Langmuir b",1.0);
              }

              aux2_list.sublist(ss.str()).sublist("function-constant")
                  .set<double>("value", langmuir_b);

            }
          }

          // Freundlich

          {
            Teuchos::ParameterList & aux3_list =
                isotherm_freundlich_n_ic.sublist("function").sublist(*ir)
                .set<std::string>("region",*ir)
                .set<std::string>("component","cell")
                .sublist("function");

            aux3_list.set<int>("Number of DoFs", comp_names.size())
                .set("Function type", "composite function");

            for ( int ic = 0; ic != comp_names.size(); ++ic) {

              std::stringstream ss;
              ss << "DoF " << ic + 1 << " Function";

              double freundlich_n(1.0);
              if (sorption_isotherms_list.isSublist(comp_names[ic])) {
                freundlich_n = sorption_isotherms_list.sublist(comp_names[ic]).get<double>("Freundlich n",1.0);
              }

              aux3_list.sublist(ss.str()).sublist("function-constant")
                  .set<double>("value", freundlich_n);

            }
          }

        }

      }

      if (sorption_site_names_.size() > 0) {



        Teuchos::ParameterList & sorption_sites_ic = chem_ic.sublist("sorption_sites");

        Teuchos::ParameterList & sorption_sites_list = matprop_list.sublist(matprop_list.name(i)).sublist("Surface Complexation Sites");

        for (Teuchos::Array<std::string>::const_iterator ir=regions.begin(); ir!=regions.end(); ir++) {

          Teuchos::ParameterList & aux1_list =
              sorption_sites_ic.sublist("function").sublist(*ir)
              .set<std::string>("region",*ir)
              .set<std::string>("component","cell")
              .sublist("function");

          aux1_list.set<int>("Number of DoFs", sorption_site_names_.size())
              .set("Function type", "composite function");

          for ( int ic = 0; ic != sorption_site_names_.size(); ++ic) {

            std::stringstream ss;
            ss << "DoF " << ic + 1 << " Function";

            double value(0.0);
            if (sorption_sites_list.isSublist(sorption_site_names_[ic])) {
              value = sorption_sites_list.sublist(sorption_site_names_[ic]).get<double>("Site Density",0.0);
            }

            aux1_list.sublist(ss.str()).sublist("function-constant")
                .set<double>("value", value);

          }
        }
      }

    }


    Teuchos::ParameterList & ic_list = plist->sublist("Initial Conditions");

    for (Teuchos::ParameterList::ConstIterator ic = ic_list.begin(); ic != ic_list.end(); ++ic) {
      if (ic_list.isSublist(ic->first)) {

        Teuchos::ParameterList & ics = ic_list.sublist(ic->first);

        Teuchos::Array<std::string> ass_regions = ics.get<Teuchos::Array<std::string> >("Assigned Regions");

        Teuchos::ParameterList &free_ion_species_ic = chem_ic.sublist("free_ion_species");

        for (Teuchos::Array<std::string>::const_iterator ir=ass_regions.begin(); ir!=ass_regions.end(); ir++) {

          Teuchos::ParameterList & aux1_list =
              free_ion_species_ic.sublist("function").sublist(*ir)
              .set<std::string>("region",*ir)
              .set<std::string>("component","cell")
              .sublist("function");

          aux1_list.set<int>("Number of DoFs", comp_names.size())
              .set("Function type", "composite function");

          for (int j = 0; j<comp_names.size(); ++j) {
            std::stringstream ss;
            ss << "DoF " << j+1 << " Function";

            double value(1.0e-9);
            value = ics.sublist("Solute IC").sublist(phase_name).sublist(phase_comp_name).sublist(comp_names[j]).sublist("IC: Uniform Concentration").get<double>("Free Ion Guess",1.0e-9);

            aux1_list.sublist(ss.str()).sublist("function-constant")
                .set<double>("value", value);
          }
        }
      }

    }
  }
  return chem_list;
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
    Exceptions::amanzi_throw(Errors::Message("The input file does not specify an \"Amanzi Input Format Version\""));
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

  if ((major != AMANZI_OLD_INPUT_VERSION_MAJOR) ||
      (minor != AMANZI_OLD_INPUT_VERSION_MINOR) ||
      (micro != AMANZI_OLD_INPUT_VERSION_MICRO)) {
    std::stringstream ss_ver_reqd;
    ss_ver_reqd << AMANZI_OLD_INPUT_VERSION_MAJOR << "." << AMANZI_OLD_INPUT_VERSION_MINOR << "." << AMANZI_OLD_INPUT_VERSION_MICRO;
    std::stringstream ss_ver_inp;
    ss_ver_inp << major << "." << minor << "." << micro;

    Exceptions::amanzi_throw(Errors::Message("The input format version "+ss_ver_inp.str()+" does not match the required version "+ss_ver_reqd.str()));
  }
}





}
}
