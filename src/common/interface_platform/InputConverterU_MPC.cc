/*
  This is the input component of the Amanzi code. 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <algorithm>
#include <sstream>
#include <string>

//TPLs
#include <xercesc/dom/DOM.hpp>
#include <xercesc/util/XMLString.hpp>

// Amanzi's
#include "errors.hh"
#include "exceptions.hh"
#include "dbc.hh"

#include "InputConverterU.hh"
#include "InputConverterU_Defs.hh"

namespace Amanzi {
namespace AmanziInput {

XERCES_CPP_NAMESPACE_USE

/* ******************************************************************
* Create MPC list, version 2, dubbed cycle driver.
****************************************************************** */
Teuchos::ParameterList InputConverterU::TranslateCycleDriver_()
{
  Teuchos::ParameterList out_list;

  if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH) {
    *vo_->os() << "Translating cycle driver" << std::endl;
  }

  // parse available PKs
  std::map<std::string, bool> pk_state;
  std::map<std::string, std::string> pk_model;

  int transient_model(0);
  char* tagname;

  DOMNodeList* node_list = doc_->getElementsByTagName(XMLString::transcode("process_kernels"));
  DOMNode* node = node_list->item(0);
  DOMNodeList* children = node->getChildNodes();

  for (int i = 0; i < children->getLength(); ++i) {
    DOMNode* inode = children->item(i);
    if (inode->getNodeType() != DOMNode::ELEMENT_NODE) continue;
    tagname = XMLString::transcode(inode->getNodeName());

    char* state = GetAttributeValueC_(static_cast<DOMElement*>(inode), "state");
    pk_state[tagname] = (strcmp(state, "on") == 0);

    if (strcmp(tagname, "flow") == 0) {
      char* model = GetAttributeValueC_(static_cast<DOMElement*>(inode), "model");
      pk_model["flow"] = TrimString_(model);
      XMLString::release(&model);

      transient_model += 4 * pk_state[tagname];
    } else if (strcmp(tagname, "chemistry") == 0) {
      char* model = GetAttributeValueC_(static_cast<DOMElement*>(inode), "engine");
      pk_model["chemistry"] = TrimString_(model);
      XMLString::release(&model);

      transient_model += pk_state[tagname];
    } else if (strcmp(tagname, "transport") == 0) {
      transient_model += 2 * pk_state[tagname];
    }
    XMLString::release(&state);
    XMLString::release(&tagname);
  }

  // parse defaults of execution_controls 
  bool flag;
  node_list = doc_->getElementsByTagName(XMLString::transcode("execution_controls"));
  node = getUniqueElementByTagNames_(node_list->item(0), "execution_control_defaults", flag);

  double t0;
  char *method, *method_d, *dt0_d, *filename;

  method_d = GetAttributeValueC_(static_cast<DOMElement*>(node), "method");
  dt0_d = GetAttributeValueC_(static_cast<DOMElement*>(node), "init_dt");

  // parse execution_control
  std::map<double, std::string> tp_method, tp_mode;
  std::map<double, double> tp_dt0, tp_t1;
  std::map<double, int> tp_max_cycles;

  children = node_list->item(0)->getChildNodes();
  for (int i = 0; i < children->getLength(); ++i) {
    DOMNode* inode = children->item(i);
    if (inode->getNodeType() != DOMNode::ELEMENT_NODE) continue;

    tagname = XMLString::transcode(inode->getNodeName());
    if (strcmp(tagname, "execution_control") != 0) continue;

    // generate fake name (will be created later)
    std::stringstream ss;
    ss << "TP " << i;
    std::string name = ss.str();

    t0 = TimeCharToValue_(GetAttributeValueC_(static_cast<DOMElement*>(inode), "start"));

    tp_mode[t0] = GetAttributeValueC_(static_cast<DOMElement*>(inode), "mode");
    tp_t1[t0] = TimeCharToValue_(GetAttributeValueC_(static_cast<DOMElement*>(inode), "end"));
    tp_method[t0] = GetAttributeValueC_(static_cast<DOMElement*>(inode), "method", false, method_d);
    tp_dt0[t0] = TimeCharToValue_(GetAttributeValueC_(static_cast<DOMElement*>(inode), "init_dt", false, dt0_d));
    tp_max_cycles[t0] = GetAttributeValueD_(static_cast<DOMElement*>(inode), "max_cycles", false, 10000000);

    filename = GetAttributeValueC_(static_cast<DOMElement*>(inode), "restart", false, NULL);
  }

  // sort time periods
  // std::sort(tp_t1.begin(), tp_t1.end());
  // std::sort(tp_mode.begin(), tp_mode.end());

  // create steady-state TP
  int tp_id(0);
  Teuchos::ParameterList pk_tree_list;

  std::map<double, std::string>::iterator it = tp_mode.begin();
  if (it->second == "steady") {
    std::ostringstream ss;
    ss << "TP " << tp_id;

    Teuchos::ParameterList& tmp_list = out_list.sublist("time periods").sublist(ss.str());
    tmp_list.sublist("PK Tree").sublist("Flow Steady").set<std::string>("PK type", pk_model["flow"]);
    tmp_list.set<double>("start period time", t0);
    tmp_list.set<double>("end period time", tp_t1.begin()->second);
    tmp_list.set<int>("maximum cycle number", tp_max_cycles.begin()->second);
    tmp_list.set<double>("initial time step", tp_dt0.begin()->second);

    tp_id++;
    it++;
  }

  // old version 
  while (it != tp_mode.end()) {
    switch (transient_model) {
    case 1:
      pk_tree_list.sublist("Chemistry").set<std::string>("PK type", "chemistry");
      break;
    case 2:
      pk_tree_list.sublist("Transport").set<std::string>("PK type", "transport");
      break;
    case 3:
      {
        Teuchos::ParameterList& tmp_list = pk_tree_list.sublist("Reactive Transport");
        tmp_list.set<std::string>("PK type", "reactive transport");
        tmp_list.sublist("Transport").set<std::string>("PK type", "transport");
        tmp_list.sublist("Chemistry").set<std::string>("PK type", "chemistry");  
        break;
      }
    case 4:
      pk_tree_list.sublist("Flow").set<std::string>("PK type", pk_model["flow"]);    
      break;
    case 5:
      {
        Teuchos::ParameterList& tmp_list = pk_tree_list.sublist("Flow and Chemistry");
        tmp_list.set<std::string>("PK type", "flow reactive transport");
        tmp_list.sublist("Chemistry").set<std::string>("PK type", "chemistry");
        tmp_list.sublist("Flow").set<std::string>("PK type", pk_model["flow"]); 
        break;
      }
    case 6:
      {
        Teuchos::ParameterList& tmp_list = pk_tree_list.sublist("Flow and Transport");
        tmp_list.set<std::string>("PK type", "flow reactive transport");
        tmp_list.sublist("Transport").set<std::string>("PK type", "transport");
        tmp_list.sublist("Flow").set<std::string>("PK type", pk_model["flow"]);
        break;
      }
    case 7:
      {
        Teuchos::ParameterList& tmp_list = pk_tree_list.sublist("Flow and Reactive Transport");
        tmp_list.set<std::string>("PK type", "flow reactive transport");
        tmp_list.sublist("Reactive Transport").set<std::string>("PK type", "reactive transport");
        tmp_list.sublist("Reactive Transport").sublist("Transport").set<std::string>("PK type", "transport");
        tmp_list.sublist("Reactive Transport").sublist("Chemistry").set<std::string>("PK type", "chemistry");
        tmp_list.sublist("Flow").set<std::string>("PK type", pk_model["flow"]);
        break;
      }
    default:
      Exceptions::amanzi_throw(Errors::Message("This model does not supported by new MPC driver."));
    }

    std::ostringstream ss;
    ss << "TP " << tp_id;

    Teuchos::ParameterList& tmp_list = out_list.sublist("time periods").sublist(ss.str());
    tmp_list.sublist("PK Tree") = pk_tree_list;
    tmp_list.set<double>("start period time", it->first);
    tmp_list.set<double>("end period time", tp_t1[it->first]);
    tmp_list.set<int>("maximum cycle number", tp_max_cycles[it->first]);
    tmp_list.set<double>("initial time step", tp_dt0[it->first]);

    tp_id++;
    it++;
  }

  if (transient_model & 2 || transient_model & 1) {
    out_list.set<Teuchos::Array<std::string> >("component names", comp_names_all_);
  }

  out_list.sublist("Time Period Control") = TranslateTimePeriodControls_();
  if (filename != NULL) {
    out_list.sublist("Restart").set<std::string>("File Name", filename);
  }
  out_list.sublist("VerboseObject") = verb_list_;

  return out_list;
}


/* ******************************************************************
* Empty
****************************************************************** */
Teuchos::ParameterList InputConverterU::TranslateTimePeriodControls_()
{
  Teuchos::ParameterList out_list;

  if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH) {
    *vo_->os() << "Translating time period controls" << std::endl;
  }

  Teuchos::Array<double> start_times;
  Teuchos::Array<double> initial_time_step;
  Teuchos::Array<double> maximum_time_step;

  std::map<double, double> time_map;
  std::map<double, double> max_dt_map;
  double default_initial_time_step(RESTART_TIME_STEP);
  double default_max_time_step(MAXIMUM_TIME_STEP);
  
  double max_cycle_number;

  /*
  Teuchos::ParameterList exe_sublist = *Teuchos::sublist(plist, "Execution Control", true);

  // get the default initial time step
  if (exe_sublist.isSublist("Time Period Control")) {
    default_initial_time_step = exe_sublist.sublist("Time Period Control").get<double>("Default Initial Time Step",RESTART_TIME_STEP);
    default_max_time_step = exe_sublist.sublist("Time Period Control").get<double>("Default Maximum Time Step", MAXIMUM_TIME_STEP);
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

            for (Teuchos::Array<double>::const_iterator times_it = times.begin(); times_it != times.end(); ++times_it) {
              // skip the first one, there is no jump
              if (times_it != times.begin()) {
                time_map[*times_it] = default_initial_time_step;
                max_dt_map[*times_it] = -1;
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

            for (Teuchos::Array<double>::const_iterator times_it = times.begin(); times_it != times.end(); ++times_it) {
              // skip the first one, there is no jump
              if (times_it != times.begin()) {
                time_map[*times_it] = default_initial_time_step;
                max_dt_map[*times_it] = -1;
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
    if (exe_sublist.sublist("Time Period Control").isParameter("Maximum Time Step")){
      maximum_time_step = exe_sublist.sublist("Time Period Control").get<Teuchos::Array<double> >("Maximum Time Step");
    }

    if (maximum_time_step.size() != initial_time_step.size()) {
      maximum_time_step.resize(initial_time_step.size());
      for (int i=0; i < maximum_time_step.size(); ++i){
        maximum_time_step[i] = default_max_time_step;
      }
    }

    Teuchos::Array<double>::const_iterator initial_time_step_it = initial_time_step.begin();
    Teuchos::Array<double>::const_iterator max_time_step_it = maximum_time_step.begin();
    for (Teuchos::Array<double>::const_iterator start_times_it = start_times.begin();
         start_times_it != start_times.end(); ++start_times_it) {
      time_map[*start_times_it] = *initial_time_step_it;
      ++initial_time_step_it;

      max_dt_map[*start_times_it] = *max_time_step_it;
      ++max_time_step_it;
    }
  }

  // delete the start, switch, and end times, since the user must specify initial time steps for those seperately
  if (exe_sublist.isSublist("Time Integration Mode")) {
    if (exe_sublist.sublist("Time Integration Mode").isSublist("Initialize To Steady")) {
      double start_time = exe_sublist.sublist("Time Integration Mode").sublist("Initialize To Steady").get<double>("Start");
      double switch_time = exe_sublist.sublist("Time Integration Mode").sublist("Initialize To Steady").get<double>("Switch");
      double end_time = exe_sublist.sublist("Time Integration Mode").sublist("Initialize To Steady").get<double>("End");

      time_map.erase(start_time);  max_dt_map.erase(start_time);
      time_map.erase(switch_time); max_dt_map.erase(switch_time);
      time_map.erase(end_time);    max_dt_map.erase(end_time);
    }
  }

  ASSERT(start_times.size() == initial_time_step.size());
  ASSERT(start_times.size() == maximum_time_step.size());

  start_times.clear();
  initial_time_step.clear();
  maximum_time_step.clear();

  for (std::map<double,double>::const_iterator map_it = time_map.begin(), max_it = max_dt_map.begin();
       map_it != time_map.end(); ++map_it, ++max_it) {
    start_times.push_back(map_it->first);
    initial_time_step.push_back(map_it->second);
    if (max_it->second < 0){
      if (max_it == max_dt_map.begin()) 
	maximum_time_step.push_back(default_max_time_step);
      else {
        int sz = maximum_time_step.size();
        if (sz > 0) maximum_time_step.push_back(maximum_time_step[sz-1]);
      }
    }
    else {
      maximum_time_step.push_back(max_it->second);
    }
  }

  out_list.set<Teuchos::Array<double> >("Start Times", start_times);
  out_list.set<Teuchos::Array<double> >("Initial Time Step", initial_time_step);
  out_list.set<Teuchos::Array<double> >("Maximum Time Step", maximum_time_step);
  */

  return out_list;
}

}  // namespace AmanziInput
}  // namespace Amanzi
