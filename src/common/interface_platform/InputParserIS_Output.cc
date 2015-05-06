#include <sstream>
#include <string>
#include <boost/bind.hpp>

#include "errors.hh"
#include "exceptions.hh"
#include "dbc.hh"

#include "InputParserIS.hh"
#include "InputParserIS_Defs.hh"

namespace Amanzi {
namespace AmanziInput {

/* ******************************************************************
* Populate visualization list.
******************************************************************* */
Teuchos::ParameterList InputParserIS::CreateVisualizationDataList_(Teuchos::ParameterList* plist)
{
  Teuchos::ParameterList vis_list;
  Teuchos::Array<double> visualizationPoints;

  if (plist->isSublist("Output")) {
    if (plist->sublist("Output").isSublist("Visualization Data")) {
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
      if (vis_list.isParameter("Cycle Macros")) {
        std::vector<std::string> cycle_macros;
       	cycle_macros = vis_list.get<Teuchos::Array<std::string> >("Cycle Macros").toVector();

	int j(0);
        for (int i = 0; i < cycle_macros.size(); i++) {
          //Teuchos::Array<int> cm = CreateCycleMacro_(cycle_macro,plist);
          Teuchos::ParameterList cycle_macro_list = CreateCycleMacro_(cycle_macros[i], plist);
          if (cycle_macro_list.isParameter("Start_Period_Stop")) {
            std::stringstream ss;
            ss << "cycles start period stop " << j;
            vis_list.set(ss.str(),cycle_macro_list.get<Teuchos::Array<int> >("Start_Period_Stop"));
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
      if (vis_list.isParameter("Time Macros")) {
        std::vector<std::string> time_macros;
        time_macros = vis_list.get<Teuchos::Array<std::string> >("Time Macros").toVector();

        int j(0);
        for (int i = 0; i < time_macros.size(); i++) {
          // Create a local parameter to store the time macro
          Teuchos::ParameterList time_macro_list = CreateTimeMacro_(time_macros[i], plist);
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
* Checkpoint list.
****************************************************************** */
Teuchos::ParameterList InputParserIS::CreateCheckpointDataList_(Teuchos::ParameterList* plist)
{
  Teuchos::ParameterList restart_list;

  if (plist->isSublist("Output")) {
    if (plist->sublist("Output").isSublist("Checkpoint Data")) {
      Teuchos::ParameterList rlist = plist->sublist("Output").sublist("Checkpoint Data");

      restart_list.set<std::string>("file name base", rlist.get<std::string>("File Name Base",std::string("checkpoint")));
      restart_list.set<int>("file name digits", rlist.get<int>("File Name Digits",5));

      // check if the cycle range is defined via a macro
      if (rlist.isParameter("Cycle Macro")) {
        std::string cycle_macro = rlist.get<std::string>("Cycle Macro");

        Teuchos::Array<int> range = CreateCycleMacro_(cycle_macro, plist).get<Teuchos::Array<int> >("Start_Period_Stop");

        restart_list.set<Teuchos::Array<int> >("cycles start period stop", range);
      }
    }
  }

  return restart_list;
}


/* ******************************************************************
* Walkabout list.
****************************************************************** */
Teuchos::ParameterList InputParserIS::CreateWalkaboutDataList_(Teuchos::ParameterList* plist)
{
  Teuchos::ParameterList walkabout_list;

  if (plist->isSublist("Output")) {
    if (plist->sublist("Output").isSublist("Walkabout Data")) {
      Teuchos::ParameterList rlist = plist->sublist("Output").sublist("Walkabout Data");

      walkabout_list.set<std::string>("file name base", rlist.get<std::string>("File Name Base",std::string("walkabout")));
      walkabout_list.set<int>("file name digits", rlist.get<int>("File Name Digits",5));

      // check if the cycle range is defined via a macro
      if (rlist.isParameter("Cycle Macro")) {
        std::string cycle_macro = rlist.get<std::string>("Cycle Macro");

        Teuchos::Array<int> range = CreateCycleMacro_(cycle_macro, plist).get<Teuchos::Array<int> >("Start_Period_Stop");

        walkabout_list.set<Teuchos::Array<int> >("cycles start period stop", range);
      }
    }
  }

  return walkabout_list;
}


/* ******************************************************************
* Observation list.
****************************************************************** */
Teuchos::ParameterList InputParserIS::CreateObservationDataList_(Teuchos::ParameterList* plist)
{
  using namespace boost;
  using boost::bind;

  // Create a parameter list for holding data
  Teuchos::ParameterList obs_list;
  Teuchos::Array<double> observationPoints;

  if (plist->isSublist("Output")) {
    if (plist->sublist("Output").isSublist("Observation Data")) {
      // Iinitialize a structure with the XML data
      Teuchos::ParameterList olist = plist->sublist("Output").sublist("Observation Data");
      // If the node has value refering to the name of the output file, grab it
      if (olist.isParameter("Observation Output Filename")) {
        obs_list.set<std::string>("Observation Output Filename", olist.get<std::string>("Observation Output Filename"));
      } else {
        Exceptions::amanzi_throw(Errors::Message("The required parameter Observation Output Filename was not specified."));
      }
      obs_list.set<int>("precision", olist.get<int>("precision", 16));
      // Iterate through the array
      for (Teuchos::ParameterList::ConstIterator i = olist.begin(); i != olist.end(); i++) {
        // If the current iteration node is a "tree"
        if (olist.isSublist(i->first)) {
          // copy the observation data sublist into the local list
          obs_list.sublist(i->first) = olist.sublist(i->first);
          
          Teuchos::Array<double> all_times;
          all_times.clear();
          if (obs_list.sublist(i->first).isParameter("Time Macros")) {
            std::vector<std::string> time_macros;
            time_macros = obs_list.sublist(i->first).get<Teuchos::Array<std::string> >("Time Macros").toVector();
            // Create a local parameter list and store the time macro (3 doubles)
            int j(0);
            for (int k = 0; k < time_macros.size(); k++) {
              // Create a local parameter to store the time macro
              Teuchos::ParameterList time_macro_list = CreateTimeMacro_(time_macros[k], plist);
              if (time_macro_list.isParameter("Start_Period_Stop")) {
                std::stringstream ss;
                ss << "times start period stop " << j;
                obs_list.sublist(i->first).set(ss.str(),time_macro_list.get<Teuchos::Array<double> >("Start_Period_Stop"));
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
            obs_list.sublist(i->first).remove("Time Macros");
          }
          if (obs_list.sublist(i->first).isParameter("Cycle Macro")) {
            std::string cycle_macro = obs_list.sublist(i->first).get<std::string>("Cycle Macro");

            Teuchos::ParameterList cycle_macro_list = CreateCycleMacro_(cycle_macro, plist);

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
          if (obs_list.sublist(i->first).isParameter("Region")) {
            std::string name = obs_list.sublist(i->first).get<std::string>("Region");
            obs_list.sublist(i->first).set<std::string>("region", name);
            obs_list.sublist(i->first).remove("Region");
            vv_obs_regions.push_back(name);
          }
          if (obs_list.sublist(i->first).isParameter("Variable")) {
            std::string name = obs_list.sublist(i->first).get<std::string>("Variable");
            obs_list.sublist(i->first).set<std::string>("variable", name);
            obs_list.sublist(i->first).remove("Variable");
          }
          if (obs_list.sublist(i->first).isParameter("Functional")) {
            std::string name = obs_list.sublist(i->first).get<std::string>("Functional");
            obs_list.sublist(i->first).set<std::string>("functional", name);
            obs_list.sublist(i->first).remove("Functional");
          }
        }
      }
    }
  }
  return obs_list;
}


/* ******************************************************************
* Create a time macro.
****************************************************************** */
Teuchos::ParameterList InputParserIS::CreateTimeMacro_(
    const std::string& macro_name, Teuchos::ParameterList* plist) 
{
  Teuchos::ParameterList time_macro;

  if (plist->sublist("Output").sublist("Time Macros").isSublist(macro_name)) {
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
  }
  else {
    std::stringstream ss;
    ss << "The time macro " << macro_name << " does not exist in the input file";
    Exceptions::amanzi_throw(Errors::Message(ss.str().c_str()));
  }

  return time_macro;
}


/* ******************************************************************
* Create cycle macro
****************************************************************** */
Teuchos::ParameterList InputParserIS::CreateCycleMacro_(
    const std::string& macro_name, Teuchos::ParameterList* plist)
{
  Teuchos::ParameterList cycle_macro;

  if (plist->sublist("Output").sublist("Cycle Macros").isSublist(macro_name)) {
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
  }
  else {
    std::stringstream ss;
    ss << "The cycle macro " << macro_name << " does not exist in the input file";
    Exceptions::amanzi_throw(Errors::Message(ss.str().c_str()));
  }

  return cycle_macro;
}



/* ******************************************************************
* Empty
****************************************************************** */
Teuchos::Array<std::string> InputParserIS::CreateVariableMacro_(
    Teuchos::Array<std::string>& macro_name, Teuchos::ParameterList* plist)
{
  std::vector<std::string> vars;

  for (int i = 0; i < macro_name.size(); i++) {
    if (plist->sublist("Output").sublist("Variable Macros").isSublist(macro_name[i])) {
      Teuchos::ParameterList& macro_list = plist->sublist("Output").sublist("Variable Macros").sublist(macro_name[i]);

      if (macro_list.isParameter("Phase")) {
        std::string macro_phase = macro_list.get<std::string>("Phase");
        if (macro_phase == "All") {
          vars.push_back(phases_[0].solute_name);
        }
        else {  // not All, must equal phase_comp_name
          if (macro_list.isParameter("Component")) {
            std::string macro_comp = macro_list.get<std::string>("Component");
            if (macro_comp == "All") {
              vars.push_back(phases_[0].solute_name);
            }
            else { // not All, must equal
              if (macro_comp != phases_[0].solute_name) {
                std::stringstream ss;
                ss << "The phase component name " << macro_comp << " is refered to in a variable macro but is not defined";
                Exceptions::amanzi_throw(Errors::Message(ss.str().c_str()));
              }
              vars.push_back(macro_comp);
            }
          }
        }
      }

      if (macro_list.isParameter("Solute")) {
        std::string macro_solute = macro_list.get<std::string>("Solute");
        if (macro_solute == "All") {
          for (int i = 0; i<comp_names_.size(); i++) vars.push_back(comp_names_[i]);
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

  for (int i = 0; i < vars.size(); i++) {
    ret_vars[i] = vars[i];
  }

  return ret_vars;
}

}  // namespace AmanziInput
}  // namespace Amanzi
