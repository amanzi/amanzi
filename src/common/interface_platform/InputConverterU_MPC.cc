/*
  This is the input component of the Amanzi code. 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Daniil Svyatskiy (original version)
           Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <algorithm>
#include <sstream>
#include <string>

//TPLs
#include <xercesc/dom/DOM.hpp>

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

  Errors::Message msg;
  if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH)
    *vo_->os() << "Translating cycle driver" << std::endl;

  MemoryManager mm;
  DOMNodeList *node_list, *children;
  DOMNode* node;
  DOMElement* element;

  // do we need to call new version of CD?
  bool flag;
  node_list = doc_->getElementsByTagName(mm.transcode("process_kernels"));
  element = static_cast<DOMElement*>(node_list->item(0));
  children = element->getElementsByTagName(mm.transcode("pk"));
  if (children->getLength() > 0) return TranslateCycleDriverNew_();

  // parse defaults of execution_controls 
  node_list = doc_->getElementsByTagName(mm.transcode("execution_controls"));
  node = GetUniqueElementByTagsString_(node_list->item(0), "execution_control_defaults", flag);
  element = static_cast<DOMElement*>(node);

  int max_cycles, max_cycles_steady;
  double t0, t1, dt0, t0_steady, t1_steady, dt0_steady, dt_max, dt_max_steady;
  char *method, *tagname;
  bool flag_steady(false); 
  std::string mode_d, method_d, dt0_d, dt_cut_d, dt_inc_d, filename, dt_max_d;

  mode_d = GetAttributeValueS_(element, "mode", TYPE_NONE, false, "");
  method_d = GetAttributeValueS_(element, "method", TYPE_NONE, false, "");
  dt0_d = GetAttributeValueS_(element, "init_dt", TYPE_TIME, false, "0.0");
  dt_max_d = GetAttributeValueS_(element, "max_dt", TYPE_TIME, false, "1.0e+99");
  dt_cut_d = GetAttributeValueS_(element, "reduction_factor", TYPE_TIME, false, "0.8");
  dt_inc_d = GetAttributeValueS_(element, "increase_factor", TYPE_TIME, false, "1.2");

  // parse execution_control
  std::map<double, std::string> tp_mode;
  std::map<double, double> tp_dt0, tp_t1, tp_dt_max;
  std::map<double, int> tp_max_cycles;

  element = static_cast<DOMElement*>(node_list->item(0));
  children = element->getElementsByTagName(mm.transcode("execution_control"));
  for (int i = 0; i < children->getLength(); ++i) {
    DOMNode* inode = children->item(i);
    if (inode->getNodeType() != DOMNode::ELEMENT_NODE) continue;
    element = static_cast<DOMElement*>(inode);

    t0 = GetAttributeValueD_(element, "start", TYPE_TIME);
    t1 = GetAttributeValueD_(element, "end", TYPE_TIME, false, t0 - 10.0);
    dt0 = TimeStringToValue_(GetAttributeValueS_(element, "init_dt", TYPE_TIME, false, dt0_d));
    dt_max = TimeStringToValue_(GetAttributeValueS_(element, "max_dt", TYPE_TIME, false, dt_max_d));
    max_cycles = GetAttributeValueL_(element, "max_cycles", TYPE_NUMERICAL, false, -1);
    std::string mode = GetAttributeValueS_(element, "mode", TYPE_NONE, false, mode_d);

    dt_cut_[mode] = TimeStringToValue_(GetAttributeValueS_(
        element, "reduction_factor", TYPE_TIME, false, dt_cut_d));
    dt_inc_[mode] = TimeStringToValue_(GetAttributeValueS_(
        element, "increase_factor", TYPE_TIME, false, dt_inc_d));

    if (mode == "steady") {
      t0_steady = t0;
      t1_steady = t1;
      dt0_steady = dt0;
      dt_max_steady = dt_max;
      max_cycles_steady = max_cycles;
      flag_steady = true;
    } else {
      if (tp_mode.find(t0) != tp_mode.end()) {
        msg << "Transient \"execution_controls\" cannot have the same start time.\n";
        Exceptions::amanzi_throw(msg);
      }  

      tp_mode[t0] = mode;
      tp_t1[t0] = t1;
      tp_dt0[t0] = dt0;
      tp_dt_max[t0] = dt_max;
      tp_max_cycles[t0] = max_cycles;

      filename = GetAttributeValueS_(element, "restart", TYPE_NONE, false, "");
    }

    if (init_filename_.size() == 0)
        init_filename_ = GetAttributeValueS_(element, "initialize", TYPE_NONE, false, "");
  }

  // populate optional end-times 
  flag = true;
  std::map<double, double>::iterator it1, it2;
  it1 = tp_t1.begin();
  if (flag_steady && t1_steady < t0_steady) { 
    if (it1 != tp_t1.end())
      t1_steady = it1->first;
    else
      flag = false;
  }

  if (it1 != tp_t1.end()) {
    for (++(it2 = it1); it2 != tp_t1.end(); it2++) {
      it1->second = it2->first;
      it1 = it2;
    }
    if (it1->second < it1->first) flag = false;
  }

  if (!flag) {
    msg << "Execution_constrols have incorrect attribute \"end\".\n";
    Exceptions::amanzi_throw(msg);
  }

  // old version 
  // -- parse available PKs
  int transient_model(0);
  std::map<std::string, bool> pk_state;

  node_list = doc_->getElementsByTagName(mm.transcode("process_kernels"));
  node = node_list->item(0);
  children = node->getChildNodes();

  for (int i = 0; i < children->getLength(); ++i) {
    DOMNode* inode = children->item(i);
    if (inode->getNodeType() != DOMNode::ELEMENT_NODE) continue;
    tagname = mm.transcode(inode->getNodeName());

    std::string state = GetAttributeValueS_(static_cast<DOMElement*>(inode), "state");
    pk_state[tagname] = (strcmp(state.c_str(), "on") == 0);

    element = static_cast<DOMElement*>(inode);
    if (strcmp(tagname, "flow") == 0) {
      flow_model_ = GetAttributeValueS_(element, "model", "richards, saturated, constant");
      pk_model_["flow"] = (flow_model_ == "richards") ? "richards" : "darcy";
      pk_master_["flow"] = true;
      if (flow_model_ != "constant") transient_model += 4 * pk_state[tagname];

    } else if (strcmp(tagname, "chemistry") == 0) {
      std::string model = GetAttributeValueS_(element, "engine", TYPE_NONE, false, "none");
      pk_model_["chemistry"] = model;
      transient_model += pk_state[tagname];

    } else if (strcmp(tagname, "transport") == 0) {
      transient_model += 2 * pk_state[tagname];
    }
  }

  // -- create steady-state TP
  int tp_id(0);
  Teuchos::ParameterList pk_tree_list;

  if (flag_steady && pk_state["flow"]) {
    if (flow_model_ == "constant") {
      if (t1_steady != t0_steady) {
        msg << "Constant flow must have end time = start time.\n";
        Exceptions::amanzi_throw(msg);
      }
      node = GetUniqueElementByTagsString_(
          "numerical_controls, unstructured_controls, unstr_steady-state_controls, unstr_initialization", flag);
      if (!flag) {
        msg << "Constant flow must have an initialization list, unless state=off.\n";
        Exceptions::amanzi_throw(msg);
      }
    }

    Teuchos::ParameterList& tmp_list = out_list.sublist("time periods").sublist("TP 0");
    tmp_list.sublist("PK Tree").sublist("Flow Steady").set<std::string>("PK type", pk_model_["flow"]);
    tmp_list.set<double>("start period time", t0_steady);
    tmp_list.set<double>("end period time", t1_steady);
    tmp_list.set<double>("initial time step", dt0_steady);
    tmp_list.set<double>("maximum time step", dt_max_steady);
    tmp_list.set<int>("maximum cycle number", max_cycles_steady);

    tp_id++;
  }

  // -- create PK tree for transient TP
  std::string submodel;
  std::map<double, std::string>::iterator it = tp_mode.begin();
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
      pk_tree_list.sublist("Flow").set<std::string>("PK type", pk_model_["flow"]);    
      break;
    case 5:
      {
        Teuchos::ParameterList& tmp_list = pk_tree_list.sublist("Flow and Chemistry");
        tmp_list.set<std::string>("PK type", "flow reactive transport");
        tmp_list.sublist("Chemistry").set<std::string>("PK type", "chemistry");
        tmp_list.sublist("Flow").set<std::string>("PK type", pk_model_["flow"]); 
        break;
      }
    case 6:
      {
        Teuchos::ParameterList& tmp_list = pk_tree_list.sublist("Flow and Transport");
        tmp_list.set<std::string>("PK type", "flow reactive transport");
        tmp_list.sublist("Transport").set<std::string>("PK type", "transport");
        tmp_list.sublist("Flow").set<std::string>("PK type", pk_model_["flow"]);
        break;
      }
    case 7:
      {
        submodel = (pk_model_["chemistry"] == "amanzi") ? "chemistry amanzi" : "chemistry alquimia";
        Teuchos::ParameterList& tmp_list = pk_tree_list.sublist("Flow and Reactive Transport");
        tmp_list.set<std::string>("PK type", "flow reactive transport");
        tmp_list.sublist("Reactive Transport").set<std::string>("PK type", "reactive transport");
        tmp_list.sublist("Reactive Transport").sublist("Transport").set<std::string>("PK type", "transport");
        tmp_list.sublist("Reactive Transport").sublist("Chemistry").set<std::string>("PK type", submodel);
        tmp_list.sublist("Flow").set<std::string>("PK type", pk_model_["flow"]);
        break;
      }
    default:
      Exceptions::amanzi_throw(Errors::Message("This model is not supported by the MPC."));
    }

    std::ostringstream ss;
    ss << "TP " << tp_id;

    Teuchos::ParameterList& tmp_list = out_list.sublist("time periods").sublist(ss.str());
    tmp_list.sublist("PK Tree") = pk_tree_list;
    tmp_list.set<double>("start period time", it->first);
    tmp_list.set<double>("end period time", tp_t1[it->first]);
    tmp_list.set<int>("maximum cycle number", tp_max_cycles[it->first]);
    tmp_list.set<double>("initial time step", tp_dt0[it->first]);
    tmp_list.set<double>("maximum time step", tp_dt_max[it->first]);

    tp_id++;
    it++;
  }

  if (transient_model & 2 || transient_model & 1) {
    out_list.set<Teuchos::Array<std::string> >("component names", comp_names_all_);
    out_list.set<int>("number of liquid components", phases_["water"].size());
  }

  out_list.sublist("time period control") = TranslateTimePeriodControls_();
  if (filename.size() > 0) {
    restart_ = true;
    out_list.sublist("restart").set<std::string>("file name", filename);
  }
  out_list.sublist("VerboseObject") = verb_list_.sublist("VerboseObject");

  return out_list;
}


/* ******************************************************************
* Create new cycle driver list.
****************************************************************** */
Teuchos::ParameterList InputConverterU::TranslateCycleDriverNew_()
{
  Teuchos::ParameterList out_list;

  Errors::Message msg;
  if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH)
    *vo_->os() << "switching to the new format of process_kernel" << std::endl;

  MemoryManager mm;
  DOMNodeList *node_list, *children;
  DOMNode* node;
  DOMElement* element;
  char* text;

  // parse execution_controls_defaults
  bool flag;
  node_list = doc_->getElementsByTagName(mm.transcode("execution_controls"));
  node = GetUniqueElementByTagsString_(node_list->item(0), "execution_control_defaults", flag);
  element = static_cast<DOMElement*>(node);

  double t0, t1, dt0, dt_max;
  char *method, *tagname;
  bool flag_steady(false); 
  std::string method_d, dt0_d, dt_max_d, mode_d, dt_cut_d, dt_inc_d, filename;

  method_d = GetAttributeValueS_(element, "method", TYPE_NONE, false, "");
  dt0_d = GetAttributeValueS_(element, "init_dt", TYPE_TIME, false, "0.0");
  dt_max_d = GetAttributeValueS_(element, "max_dt", TYPE_TIME, false, "1.0e+99");
  mode_d = GetAttributeValueS_(element, "mode", TYPE_NONE, false, "");
  dt_cut_d = GetAttributeValueS_(element, "reduction_factor", TYPE_TIME, false, "0.8");
  dt_inc_d = GetAttributeValueS_(element, "increase_factor", TYPE_TIME, false, "1.2");

  // Logic behind attribute "mode" in the new PK struncture is not clear yet, 
  // so that we set up some defaults.
  dt_cut_["steady"] = 0.8;
  dt_inc_["steady"] = 1.2;

  dt_cut_["transient"] = 0.8;
  dt_inc_["transient"] = 1.2;

  // parse execution_control
  std::map<std::string, double> tp_t0, tp_t1, tp_dt0;
  std::map<std::string, int> tp_max_cycles;
  std::map<std::string, double> tp_max_dt;

  children = node_list->item(0)->getChildNodes();
  for (int i = 0; i < children->getLength(); ++i) {
    DOMNode* inode = children->item(i);
    if (inode->getNodeType() != DOMNode::ELEMENT_NODE) continue;
    element = static_cast<DOMElement*>(inode);

    tagname = mm.transcode(inode->getNodeName());
    if (strcmp(tagname, "execution_control") == 0) {
      t0 = GetAttributeValueD_(element, "start", TYPE_TIME);
      t1 = GetAttributeValueD_(element, "end", TYPE_TIME);
      dt0 = TimeStringToValue_(GetAttributeValueS_(element, "init_dt", TYPE_TIME, false, dt0_d));
      dt_max = TimeStringToValue_(GetAttributeValueS_(element, "max_dt", TYPE_TIME, false, dt_max_d));
      std::string mode = GetAttributeValueS_(element, "mode", TYPE_NONE, false, mode_d);

      tp_t0[mode] = t0;
      tp_t1[mode] = t1;
      tp_dt0[mode] = dt0;
      tp_max_dt[mode] = dt_max;
      tp_max_cycles[mode] = GetAttributeValueL_(element, "max_cycles", TYPE_NUMERICAL, false, -1);
      dt_cut_[mode] = TimeStringToValue_(GetAttributeValueS_(
          element, "reduction_factor", TYPE_TIME, false, dt_cut_d));
      dt_inc_[mode] = TimeStringToValue_(GetAttributeValueS_(
          element, "increase_factor", TYPE_TIME, false, dt_inc_d));

      filename = GetAttributeValueS_(element, "restart", TYPE_NONE, false, "");
    }
  }

  // new version of process_kernels
  // -- parse available PKs
  int tp_id(0);
  std::string model, state, pkname, strong_name, weak_name;
  std::vector<std::string> pks_strong, pks_weak;
  std::map<std::string, bool> pk_state;

  node_list = doc_->getElementsByTagName(mm.transcode("process_kernels"));
  element = static_cast<DOMElement*>(node_list->item(0));
  DOMNodeList* pks = element->getElementsByTagName(mm.transcode("pk"));

  for (int i = 0; i < pks->getLength(); ++i) {
    DOMNode* inode = pks->item(i);
    if (inode->getNodeType() != DOMNode::ELEMENT_NODE) continue;
    children = inode->getChildNodes();
    std::string mode = GetAttributeValueS_(static_cast<DOMElement*>(inode), "mode");

    // collect active pks and coupling of pks
    int transient_model(0);
    pk_state.clear();
    pk_model_.clear();
    pk_master_.clear();
    for (int j = 0; j < children->getLength(); ++j) {
      DOMNode* jnode = children->item(j);
      if (jnode->getNodeType() != DOMNode::ELEMENT_NODE) continue;
      tagname = mm.transcode(jnode->getNodeName());

      element = static_cast<DOMElement*>(jnode);
      if (strcmp(tagname, "flow") == 0) {
        flow_model_ = GetAttributeValueS_(element, "model", "richards, saturated, constant");
        pk_model_["flow"] = (flow_model_ == "richards") ? "richards" : "darcy";
        pk_master_["flow"] = true;
        state = GetAttributeValueS_(element, "state");
        pk_state["flow"] = (strcmp(state.c_str(), "on") == 0);
        transient_model += 4;
      
      } else if (strcmp(tagname, "chemistry") == 0) {
        model = GetAttributeValueS_(element, "engine");
        pk_model_["chemistry"] = model;
        GetAttributeValueS_(element, "state", "on");
        transient_model += 1;

      } else if (strcmp(tagname, "transport") == 0) {
        GetAttributeValueS_(element, "state", "on");
        transient_model += 2;

      } else if (strcmp(tagname, "energy") == 0) {
        model = GetAttributeValueS_(element, "model");
        pk_model_["energy"] = model;
        pk_master_["energy"] = true;
        GetAttributeValueS_(element, "state", "on");
        transient_model += 8;
      } 
    }

    // we allow so far only one strongly coupled MPC
    pks_strong.clear();
    node = GetUniqueElementByTagsString_(inode, "strongly_coupled", flag);
    if (flag) {
      pkname = GetAttributeValueS_(static_cast<DOMElement*>(node), "name");
      pks_strong = CharToStrings_(mm.transcode(node->getTextContent()));
      strong_name = mm.transcode(node->getNodeName());
    }

    // we allow so far only one weakly coupled MPC
    pks_weak.clear();
    node = GetUniqueElementByTagsString_(inode, "weakly_coupled", flag);
    if (flag) {
      pkname = GetAttributeValueS_(static_cast<DOMElement*>(node), "name");
      pks_weak = CharToStrings_(mm.transcode(node->getTextContent()));
      weak_name = mm.transcode(node->getNodeName());
    }

    // create TP
    Teuchos::ParameterList pk_tree_list;
    std::ostringstream ss;
    ss << "TP " << tp_id;

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
      pk_tree_list.sublist("Flow").set<std::string>("PK type", pk_model_["flow"]);    
      break;
    case 5:
      {
        Teuchos::ParameterList& tmp_list = pk_tree_list.sublist("Flow and Chemistry");
        tmp_list.set<std::string>("PK type", "flow reactive transport");
        tmp_list.sublist("Chemistry").set<std::string>("PK type", "chemistry");
        tmp_list.sublist("Flow").set<std::string>("PK type", pk_model_["flow"]); 
        break;
      }
    case 6:
      {
        Teuchos::ParameterList& tmp_list = pk_tree_list.sublist("Flow and Transport");
        tmp_list.set<std::string>("PK type", "flow reactive transport");
        tmp_list.sublist("Transport").set<std::string>("PK type", "transport");
        tmp_list.sublist("Flow").set<std::string>("PK type", pk_model_["flow"]);
        break;
      }
    case 7:
      {
        Teuchos::ParameterList& tmp_list = pk_tree_list.sublist("Flow and Reactive Transport");
        tmp_list.set<std::string>("PK type", "flow reactive transport");
        tmp_list.sublist("Reactive Transport").set<std::string>("PK type", "reactive transport");
        tmp_list.sublist("Reactive Transport").sublist("Transport").set<std::string>("PK type", "transport");
        tmp_list.sublist("Reactive Transport").sublist("Chemistry").set<std::string>("PK type", "chemistry");
        tmp_list.sublist("Flow").set<std::string>("PK type", pk_model_["flow"]);
        break;
      }
    case 12: 
      {
        pk_master_["thermal richards"] = true;

        Teuchos::ParameterList& tmp_list = pk_tree_list.sublist("Flow and Energy");
        tmp_list.set<std::string>("PK type", "thermal richards");
        tmp_list.sublist("Flow").set<std::string>("PK type", pk_model_["flow"]);
        tmp_list.sublist("Energy").set<std::string>("PK type", pk_model_["energy"]);
        break;
      }
    default:
      Exceptions::amanzi_throw(Errors::Message("This model is not supported by the MPC."));
    }

    Teuchos::ParameterList& tmp_list = out_list.sublist("time periods").sublist(ss.str());
    tmp_list.sublist("PK Tree") = pk_tree_list;
    tmp_list.set<double>("start period time", tp_t0[mode]);
    tmp_list.set<double>("end period time", tp_t1[mode]);
    tmp_list.set<int>("maximum cycle number", tp_max_cycles[mode]);
    tmp_list.set<double>("initial time step", tp_dt0[mode]);
    tmp_list.set<double>("maximum time step", tp_max_dt[mode]);

    tp_id++;
  }

  out_list.set<Teuchos::Array<std::string> >("component names", comp_names_all_);

  out_list.sublist("time period control") = TranslateTimePeriodControls_();
  if (filename.size() > 0) {
    out_list.sublist("restart").set<std::string>("file name", filename);
  }
  out_list.sublist("VerboseObject") = verb_list_.sublist("VerboseObject");

  return out_list;
}


/* ******************************************************************
* Generic time period control list that can be atatched to any PK.
* PK specific extension are included at the end.
****************************************************************** */
Teuchos::ParameterList InputConverterU::TranslateTimePeriodControls_()
{
  Teuchos::ParameterList out_list;

  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH)
      *vo_->os() << "Translating time period controls" << std::endl;

  MemoryManager mm;
  DOMNodeList *node_list, *children;
  DOMNode* node;
  DOMElement* element;

  // get the default time steps
  bool flag;
  node = GetUniqueElementByTagsString_("execution_controls, execution_control_defaults", flag);

  double t, dt_init_d, dt_max_d;
  std::map<double, double> init_dt, max_dt;

  element = static_cast<DOMElement*>(node);
  dt_init_d = GetAttributeValueD_(element, "init_dt", TYPE_TIME, false, RESTART_TIMESTEP);
  dt_max_d = GetAttributeValueD_(element, "max_dt", TYPE_TIME, false, MAXIMUM_TIMESTEP);

  children = doc_->getElementsByTagName(mm.transcode("execution_control"));
  for (int i = 0; i < children->getLength(); ++i) {
    DOMNode* inode = children->item(i);
    if (inode->getNodeType() != DOMNode::ELEMENT_NODE) continue;
    element = static_cast<DOMElement*>(inode);

    t = GetAttributeValueD_(element, "start", TYPE_TIME);
    double dt_init = GetAttributeValueD_(element, "init_dt", TYPE_TIME, false, dt_init_d);
    double dt_max = GetAttributeValueD_(element, "max_dt", TYPE_TIME, false, dt_max_d);
    init_dt[t] = dt_init;
    max_dt[t] = dt_max;
  }

  // add start times of all boundary conditions to the list
  std::map<double, double> dt_init_map, dt_max_map;

  std::vector<std::string> bc_names;
  bc_names.push_back("hydrostatic");
  bc_names.push_back("uniform_pressure");
  bc_names.push_back("inward_mass_flux");
  bc_names.push_back("outward_mass_flux");
  bc_names.push_back("inward_volumetric_flux");
  bc_names.push_back("outward_volumetric_flux");
  bc_names.push_back("seepage_face");
  bc_names.push_back("aqueous_conc");
  bc_names.push_back("uniform_conc");
  bc_names.push_back("constraint");
  bc_names.push_back("diffusion_dominated_release");
  bc_names.push_back("uniform_temperature");

  node_list = doc_->getElementsByTagName(mm.transcode("boundary_conditions"));
  if (node_list->getLength() > 0) {
    node = node_list->item(0);

    for (int n = 0; n < bc_names.size(); ++n) {
      children = static_cast<DOMElement*>(node)->getElementsByTagName(mm.transcode(bc_names[n].c_str()));
      for (int i = 0; i < children->getLength(); ++i) {
        DOMNode* inode = children->item(i);
        if (inode->getNodeType() != DOMNode::ELEMENT_NODE) continue;

        t = GetAttributeValueD_(static_cast<DOMElement*>(inode), "start", TYPE_TIME);
        // find position before t
        std::map<double, double>::iterator it = init_dt.upper_bound(t);
        if (it == init_dt.end()) it--;
        dt_init_map[t] = it->second;

        it = max_dt.upper_bound(t);
        if (it == max_dt.end()) it--;
        dt_max_map[t] = it->second;
      }
    }
  }

  // add start times of all sources to the list
  std::vector<std::string> src_names;
  src_names.push_back("volume_weighted");
  src_names.push_back("perm_weighted");
  src_names.push_back("uniform");
  src_names.push_back("flow_weighted_conc");

  node_list = doc_->getElementsByTagName(mm.transcode("sources"));
  if (node_list->getLength() > 0) {
    node = node_list->item(0);

    for (int n = 0; n < src_names.size(); ++n) {
      children = static_cast<DOMElement*>(node)->getElementsByTagName(mm.transcode(src_names[n].c_str()));
      for (int i = 0; i < children->getLength(); ++i) {
        DOMNode* inode = children->item(i);
        if (inode->getNodeType() != DOMNode::ELEMENT_NODE) continue;

        t = GetAttributeValueD_(static_cast<DOMElement*>(inode), "start", TYPE_TIME);
        // find position before t
        std::map<double, double>::iterator it = init_dt.upper_bound(t);
        if (it == init_dt.end()) it--;
        dt_init_map[t] = it->second;

        it = max_dt.upper_bound(t);
        if (it == max_dt.end()) it--;
        dt_max_map[t] = it->second;
      }
    }
  }

  // save times in the XML, skipping TP start times
  std::vector<double> times, dt_init, dt_max;

  for (std::map<double, double>::const_iterator it = dt_init_map.begin(), max_it = dt_max_map.begin();
       it != dt_init_map.end(); ++it, ++max_it) {
    if (init_dt.find(it->first) == init_dt.end()) {
      times.push_back(it->first);
      dt_init.push_back(it->second);
      dt_max.push_back(max_it->second);
    }
  }

  out_list.set<Teuchos::Array<double> >("start times", times);
  out_list.set<Teuchos::Array<double> >("initial time step", dt_init);
  out_list.set<Teuchos::Array<double> >("maximum time step", dt_max);

  if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH)
      *vo_->os() << "created " << dt_max.size() << " special times" << std::endl;

  return out_list;
}


/* ******************************************************************
* Translate PKs list
****************************************************************** */
Teuchos::ParameterList InputConverterU::TranslatePKs_(const Teuchos::ParameterList& cd_list)
{
  Teuchos::ParameterList out_list;

  if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH)
      *vo_->os() << "Translating process kernels" << std::endl;
  Teuchos::OSTab tab = vo_->getOSTab();

  // create PKs list
  Teuchos::ParameterList tp_list = cd_list.sublist("time periods");

  for (Teuchos::ParameterList::ConstIterator it = tp_list.begin(); it !=tp_list.end(); ++it) {
    if ((it->second).isList()) {
      Teuchos::ParameterList& pk_tree = tp_list.sublist(it->first).sublist("PK Tree");
      RegisterPKsList_(pk_tree, out_list);

      if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH) {
        std::string name = pk_tree.begin()->first;
        *vo_->os() << "PK name=\"" << name << "\", factory: \"" 
                   << pk_tree.sublist(name).get<std::string>("PK type") << "\"" << std::endl;
      }
    }
  }

  // parse list of supported PKs
  for (Teuchos::ParameterList::ConstIterator it = out_list.begin(); it != out_list.end(); ++it) {
    if ((it->second).isList()) {
      if (it->first == "Flow Steady") {
        out_list.sublist(it->first) = TranslateFlow_("steady");
      }
      else if (it->first == "Flow") {
        out_list.sublist(it->first) = TranslateFlow_("transient");
      }
      else if (it->first == "Energy") {
        out_list.sublist(it->first) = TranslateEnergy_();
      }
      else if (it->first == "Transport") {
        out_list.sublist(it->first) = TranslateTransport_();
      }
      else if (it->first == "Chemistry") {
        out_list.sublist(it->first) = TranslateChemistry_();
      }
      else if (it->first == "Reactive Transport") {
        Teuchos::Array<std::string> pk_names;
        pk_names.push_back("Chemistry");
        pk_names.push_back("Transport");
        out_list.sublist(it->first).set<Teuchos::Array<std::string> >("PKs order", pk_names);
      }
      else if (it->first == "Flow and Reactive Transport") {
        Teuchos::Array<std::string> pk_names;
        pk_names.push_back("Flow");
        pk_names.push_back("Reactive Transport");
        out_list.sublist(it->first).set<Teuchos::Array<std::string> >("PKs order", pk_names);
        out_list.sublist(it->first).set<int>("master PK index", 0);
      }
      else if (it->first == "Flow and Transport") {
        Teuchos::Array<std::string> pk_names;
        pk_names.push_back("Flow");
        pk_names.push_back("Transport");
        out_list.sublist(it->first).set<Teuchos::Array<std::string> >("PKs order", pk_names);
        out_list.sublist(it->first).set<int>("master PK index", 0);
      }
      else if (it->first == "Flow and Energy") {
        Teuchos::Array<std::string> pk_names;
        pk_names.push_back("Flow");
        pk_names.push_back("Energy");
        out_list.sublist(it->first).set<Teuchos::Array<std::string> >("PKs order", pk_names);
        out_list.sublist(it->first).set<int>("master PK index", 0);

        if (pk_master_.find("thermal richards") != pk_master_.end()) {
          // we use steady defaults so far
          out_list.sublist(it->first).sublist("time integrator") = TranslateTimeIntegrator_(
              "pressure, temperature", "nka", false,
              "unstructured_controls, unstr_thermal_richards_controls",
              TI_TS_REDUCTION_FACTOR, TI_TS_INCREASE_FACTOR);  
          out_list.sublist(it->first).sublist("VerboseObject") = verb_list_.sublist("VerboseObject");
        }
      }
    }
  }

  return out_list;
}


/* ******************************************************************
* Empty
****************************************************************** */
void InputConverterU::RegisterPKsList_(
    Teuchos::ParameterList& pk_tree, Teuchos::ParameterList& pks_list)
{
  for (Teuchos::ParameterList::ConstIterator it = pk_tree.begin(); it !=pk_tree.end();++it) {
    if ((it->second).isList()) {
      pks_list.sublist(it->first);
      RegisterPKsList_(pk_tree.sublist(it->first), pks_list);
    }   
  }
}

}  // namespace AmanziInput
}  // namespace Amanzi
