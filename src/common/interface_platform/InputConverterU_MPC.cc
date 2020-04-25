/*
  Input Converter

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
#include <climits>

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

  int max_cycles, max_cycles_steady;
  double t0, t1, dt0, t0_steady, t1_steady, dt0_steady, dt_max, dt_max_steady;
  char *tagname;
  bool flag_steady(false); 
  std::string mode_d, method_d, dt0_d, dt_cut_d, dt_inc_d, dt_max_d;

  mode_d = GetAttributeValueS_(node, "mode", TYPE_NONE, false, "");
  method_d = GetAttributeValueS_(node, "method", TYPE_NONE, false, "");
  dt0_d = GetAttributeValueS_(node, "init_dt", TYPE_TIME, false, "0.0");
  dt_max_d = GetAttributeValueS_(node, "max_dt", TYPE_TIME, false, "1.0e+99");
  dt_cut_d = GetAttributeValueS_(node, "reduction_factor", TYPE_TIME, false, "0.8");
  dt_inc_d = GetAttributeValueS_(node, "increase_factor", TYPE_TIME, false, "1.2");

  // parse execution_control
  std::string unit;
  std::map<double, std::string> tp_mode;
  std::map<double, double> tp_dt0, tp_t1, tp_dt_max;
  std::map<double, int> tp_max_cycles;

  element = static_cast<DOMElement*>(node_list->item(0));
  children = element->getElementsByTagName(mm.transcode("execution_control"));
  for (int i = 0; i < children->getLength(); ++i) {
    DOMNode* inode = children->item(i);
    if (inode->getNodeType() != DOMNode::ELEMENT_NODE) continue;

    t0 = GetAttributeValueD_(inode, "start", TYPE_TIME, DVAL_MIN, DVAL_MAX, "s");
    t1 = GetAttributeValueD_(inode, "end", TYPE_TIME, DVAL_MIN, DVAL_MAX, "s", false, t0 - 10.0);
    dt0 = ConvertUnits_(GetAttributeValueS_(inode, "init_dt", TYPE_TIME, false, dt0_d), unit);
    dt_max = ConvertUnits_(GetAttributeValueS_(inode, "max_dt", TYPE_TIME, false, dt_max_d), unit);
    max_cycles = GetAttributeValueL_(inode, "max_cycles", TYPE_NUMERICAL, -1, INT_MAX, false, -1);
    std::string mode = GetAttributeValueS_(inode, "mode", TYPE_NONE, false, mode_d);

    if (mode != "steady" && mode != "transient") {
      msg << "\"execution_controls\" has incorrect mode=" << mode << ".\n";
      Exceptions::amanzi_throw(msg);
    }  

    dt_cut_[mode] = ConvertUnits_(GetAttributeValueS_(
        inode, "reduction_factor", TYPE_TIME, false, dt_cut_d), unit);
    dt_inc_[mode] = ConvertUnits_(GetAttributeValueS_(
        inode, "increase_factor", TYPE_TIME, false, dt_inc_d), unit);

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
    }
  }

  std::string filename;
  node = GetUniqueElementByTagsString_("execution_controls, restart", flag);
  if (flag) {
    filename = GetTextContentS_(node, "", false);
    if (filename.size() == 0) ThrowErrorIllformed_("execution_controls", "restart", "filename");
  }

  node = GetUniqueElementByTagsString_("execution_controls, initialize", flag);
  if (flag) {
    init_filename_ = GetTextContentS_(node, "", false);
    if (init_filename_.size() == 0) ThrowErrorIllformed_("execution_controls", "initialize", "filename");
  }

  // time for initial conditions
  if (flag_steady) {
    ic_time_flow_ = t0_steady;
    ic_time_ = t0_steady;
    if (tp_mode.size() > 0) ic_time_ = tp_mode.begin()->first;
  } else {
    ic_time_flow_ = tp_mode.begin()->first;
    ic_time_ = ic_time_flow_;
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

    std::string state = GetAttributeValueS_(inode, "state");
    pk_state[tagname] = (strcmp(state.c_str(), "on") == 0);

    if (strcmp(tagname, "flow") == 0) {
      flow_model_ = GetAttributeValueS_(inode, "model", "richards, saturated, constant");
      pk_model_["flow"] = (flow_model_ == "richards") ? "richards" : "darcy";
      pk_master_["flow"] = true;
      if (flow_model_ != "constant") transient_model += 4 * pk_state[tagname];

    } else if (strcmp(tagname, "chemistry") == 0) {
      std::string model = GetAttributeValueS_(inode, "engine", TYPE_NONE, false, "none");
      pk_model_["chemistry"] = model;
      transient_model += pk_state[tagname];

    } else if (strcmp(tagname, "transport") == 0) {
      transient_model += 2 * pk_state[tagname];
    }
  }

  // -- create steady-state TP
  GetUniqueElementByTagsString_("fracture_network", coupled_flow_);
  coupled_transport_ = coupled_flow_;

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
        msg << "Constant flow must have unstr_steady-state_controls->unstr_initialization list, unless state=off.\n";
        Exceptions::amanzi_throw(msg);
      }
    }

    Teuchos::ParameterList& tmp_list = out_list.sublist("time periods").sublist("TP 0");
    if (!coupled_flow_) {
      tmp_list.sublist("PK tree").sublist("flow steady").set<std::string>("PK type", pk_model_["flow"]);
    } else {
      Teuchos::ParameterList& aux_list = tmp_list.sublist("PK tree").sublist("coupled flow steady");
      aux_list.set<std::string>("PK type", "darcy matrix fracture");
      aux_list.sublist("flow matrix steady").set<std::string>("PK type", "darcy");
      aux_list.sublist("flow fracture steady").set<std::string>("PK type", "darcy");
    }
    tmp_list.set<double>("start period time", t0_steady);
    tmp_list.set<double>("end period time", t1_steady);
    tmp_list.set<double>("initial time step", dt0_steady);
    tmp_list.set<double>("maximum time step", dt_max_steady);
    tmp_list.set<int>("maximum cycle number", max_cycles_steady);

    tp_id++;
  }

  // -- create PK tree for transient TP
  node = GetUniqueElementByTagsString_(
      "numerical_controls, unstructured_controls, unstr_transport_controls, algorithm", flag);
  if (flag) {
    std::string algorithm = TrimString_(mm.transcode(node->getTextContent()));
    transport_implicit_ = (algorithm == "implicit");
  }

  std::string submodel;
  std::map<double, std::string>::iterator it = tp_mode.begin();
  while (it != tp_mode.end()) {
    switch (transient_model) {
    case 1:
      PopulatePKTree_(pk_tree_list, "chemistry");
      break;
    case 2:
      if (!coupled_transport_)
        PopulatePKTree_(pk_tree_list, "transport");
      else
        PopulatePKTree_(pk_tree_list, "coupled transport");
      break;
    case 3:
      if (!coupled_transport_)
        PopulatePKTree_(pk_tree_list, "reactive transport");
      else
        PopulatePKTree_(pk_tree_list, "coupled reactive transport");
      break;
    case 4:
      if (!coupled_flow_)
        PopulatePKTree_(pk_tree_list, "flow");
      else
        PopulatePKTree_(pk_tree_list, "coupled flow");
      break;
    case 5:
      PopulatePKTree_(pk_tree_list, "flow and chemistry");
      break;
    case 6:
      if (!coupled_flow_)
        PopulatePKTree_(pk_tree_list, "flow and transport");
      else
        PopulatePKTree_(pk_tree_list, "coupled flow and transport");
      break;
    case 7:
      PopulatePKTree_(pk_tree_list, "flow and reactive transport");
      break;
    default:
      Exceptions::amanzi_throw(Errors::Message("This model is not supported by the MPC."));
    }

    std::ostringstream ss;
    ss << "TP " << tp_id;

    Teuchos::ParameterList& tmp_list = out_list.sublist("time periods").sublist(ss.str());
    tmp_list.sublist("PK tree") = pk_tree_list;
    tmp_list.set<double>("start period time", it->first);
    tmp_list.set<double>("end period time", tp_t1[it->first]);
    tmp_list.set<int>("maximum cycle number", tp_max_cycles[it->first]);
    tmp_list.set<double>("initial time step", tp_dt0[it->first]);
    tmp_list.set<double>("maximum time step", tp_dt_max[it->first]);

    tp_id++;
    it++;
  }

  if (transient_model & 2 || transient_model & 1) {
    // names
    out_list.set<Teuchos::Array<std::string> >("component names", comp_names_all_);
    out_list.set<int>("number of liquid components", phases_["water"].size());

    // available molar masses
    std::string name;
    std::vector<double> tmp(comp_names_all_.size(), 1.0);

    for (int i = 0; i < comp_names_all_.size(); ++i) {
      name = comp_names_all_[i];
      if (solute_molar_mass_.find(name) != solute_molar_mass_.end()) tmp[i] = solute_molar_mass_[name];
    }
    out_list.set<Teuchos::Array<double> >("component molar masses", tmp);
  }

  out_list.sublist("time period control") = TranslateTimePeriodControls_();
  if (filename.size() > 0) {
    restart_ = true;
    out_list.sublist("restart").set<std::string>("file name", filename);
  }
  out_list.sublist("verbose object") = verb_list_.sublist("verbose object");

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

  // parse execution_controls_defaults
  bool flag;
  node_list = doc_->getElementsByTagName(mm.transcode("execution_controls"));
  node = GetUniqueElementByTagsString_(node_list->item(0), "execution_control_defaults", flag);

  double t0, t1, dt0, dt_max;
  char *tagname;
  std::string method_d, dt0_d, dt_max_d, mode_d, dt_cut_d, dt_inc_d;

  method_d = GetAttributeValueS_(node, "method", TYPE_NONE, false, "");
  dt0_d = GetAttributeValueS_(node, "init_dt", TYPE_TIME, false, "0.0");
  dt_max_d = GetAttributeValueS_(node, "max_dt", TYPE_TIME, false, "1.0e+99");
  mode_d = GetAttributeValueS_(node, "mode", TYPE_NONE, false, "");
  dt_cut_d = GetAttributeValueS_(node, "reduction_factor", TYPE_TIME, false, "0.8");
  dt_inc_d = GetAttributeValueS_(node, "increase_factor", TYPE_TIME, false, "1.2");

  // Logic behind attribute "mode" in the new PK structure is not clear yet, 
  // so that we set up some defaults.
  dt_cut_["steady"] = 0.8;
  dt_inc_["steady"] = 1.2;

  dt_cut_["transient"] = 0.8;
  dt_inc_["transient"] = 1.2;

  // parse execution_control
  std::string unit;
  std::map<std::string, double> tp_t0, tp_t1, tp_dt0;
  std::map<std::string, int> tp_max_cycles;
  std::map<std::string, double> tp_max_dt;

  children = node_list->item(0)->getChildNodes();
  for (int i = 0; i < children->getLength(); ++i) {
    DOMNode* inode = children->item(i);
    if (inode->getNodeType() != DOMNode::ELEMENT_NODE) continue;

    tagname = mm.transcode(inode->getNodeName());
    if (strcmp(tagname, "execution_control") == 0) {
      t0 = GetAttributeValueD_(inode, "start", TYPE_TIME, DVAL_MIN, DVAL_MAX, "s");
      t1 = GetAttributeValueD_(inode, "end", TYPE_TIME, DVAL_MIN, DVAL_MAX, "s");
      dt0 = ConvertUnits_(GetAttributeValueS_(inode, "init_dt", TYPE_TIME, false, dt0_d), unit);
      dt_max = ConvertUnits_(GetAttributeValueS_(inode, "max_dt", TYPE_TIME, false, dt_max_d), unit);
      std::string mode = GetAttributeValueS_(inode, "mode", TYPE_NONE, false, mode_d);

      tp_t0[mode] = t0;
      tp_t1[mode] = t1;
      tp_dt0[mode] = dt0;
      tp_max_dt[mode] = dt_max;
      tp_max_cycles[mode] = GetAttributeValueL_(inode, "max_cycles", TYPE_NUMERICAL, -1, INT_MAX, false, -1);
      dt_cut_[mode] = ConvertUnits_(GetAttributeValueS_(
          inode, "reduction_factor", TYPE_TIME, false, dt_cut_d), unit);
      dt_inc_[mode] = ConvertUnits_(GetAttributeValueS_(
          inode, "increase_factor", TYPE_TIME, false, dt_inc_d), unit);
    }
  }

  std::string filename;
  node = GetUniqueElementByTagsString_("execution_controls, restart", flag);
  if (flag) {
    filename = GetTextContentS_(node, "", false);
    if (filename.size() == 0) ThrowErrorIllformed_("execution_controls", "restart", "filename");
  }

  node = GetUniqueElementByTagsString_("execution_controls, initialize", flag);
  if (flag) {
    init_filename_ = GetTextContentS_(node, "", false);
    if (init_filename_.size() == 0) ThrowErrorIllformed_("execution_controls", "initialize", "filename");
  }

  GetUniqueElementByTagsString_("fracture_network", coupled_flow_);
  coupled_transport_ = coupled_flow_;

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
    std::string mode = GetAttributeValueS_(inode, "mode");

    // collect active pks and coupling of pks
    int transient_model(0);
    pk_state.clear();
    pk_model_.clear();
    pk_master_.clear();
    for (int j = 0; j < children->getLength(); ++j) {
      DOMNode* jnode = children->item(j);
      if (jnode->getNodeType() != DOMNode::ELEMENT_NODE) continue;
      tagname = mm.transcode(jnode->getNodeName());

      if (strcmp(tagname, "flow") == 0) {
        flow_model_ = GetAttributeValueS_(jnode, "model", "richards, saturated, constant");
        pk_model_["flow"] = (flow_model_ == "richards") ? "richards" : "darcy";
        pk_master_["flow"] = true;
        state = GetAttributeValueS_(jnode, "state");
        pk_state["flow"] = (strcmp(state.c_str(), "on") == 0);
        transient_model += 4;
      
      } else if (strcmp(tagname, "chemistry") == 0) {
        model = GetAttributeValueS_(jnode, "engine");
        pk_model_["chemistry"] = model;
        GetAttributeValueS_(jnode, "state", "on");
        transient_model += 1;

      } else if (strcmp(tagname, "transport") == 0) {
        GetAttributeValueS_(jnode, "state", "on");
        transient_model += 2;

      } else if (strcmp(tagname, "energy") == 0) {
        model = GetAttributeValueS_(jnode, "model");
        pk_model_["energy"] = model;
        pk_master_["energy"] = true;
        GetAttributeValueS_(jnode, "state", "on");
        transient_model += 8;

      } else if (strcmp(tagname, "multiphase") == 0) {
        GetAttributeValueS_(jnode, "state", "on");
        model = GetAttributeValueS_(jnode, "model");
        pk_model_["multiphase"] = model;
        pk_master_["multiphase"] = true;
        transient_model += 16;
      } 
    }

    // we allow so far only one strongly coupled MPC
    pks_strong.clear();
    node = GetUniqueElementByTagsString_(inode, "strongly_coupled", flag);
    if (flag) {
      pkname = GetAttributeValueS_(node, "name");
      pks_strong = CharToStrings_(mm.transcode(node->getTextContent()));
      strong_name = mm.transcode(node->getNodeName());
    }

    // we allow so far only one weakly coupled MPC
    pks_weak.clear();
    node = GetUniqueElementByTagsString_(inode, "weakly_coupled", flag);
    if (flag) {
      pkname = GetAttributeValueS_(node, "name");
      pks_weak = CharToStrings_(mm.transcode(node->getTextContent()));
      weak_name = mm.transcode(node->getNodeName());
    }

    // create TP
    Teuchos::ParameterList pk_tree_list;
    std::ostringstream ss;
    ss << "TP " << tp_id;

    switch (transient_model) {
    case 1:
      PopulatePKTree_(pk_tree_list, "chemistry");
      break;
    case 2:
      if (!coupled_transport_)
        PopulatePKTree_(pk_tree_list, "transport");
      else
        PopulatePKTree_(pk_tree_list, "coupled transport");
      break;
    case 3:
      if (!coupled_transport_)
        PopulatePKTree_(pk_tree_list, "reactive transport");
      else
        PopulatePKTree_(pk_tree_list, "coupled reactive transport");
      break;
    case 4:
      if (!coupled_flow_)
        PopulatePKTree_(pk_tree_list, "flow");
      else
        PopulatePKTree_(pk_tree_list, "coupled flow");
      break;
    case 5:
      PopulatePKTree_(pk_tree_list, "flow and chemistry");
      break;
    case 6:
      if (!coupled_flow_)
        PopulatePKTree_(pk_tree_list, "flow and transport");
      else
        PopulatePKTree_(pk_tree_list, "coupled flow and transport");
      break;
    case 7:
      PopulatePKTree_(pk_tree_list, "flow and reactive transport");
      break;
    case 8:
      PopulatePKTree_(pk_tree_list, "energy");
      break;
    case 12: 
      pk_master_["thermal richards"] = true;
      PopulatePKTree_(pk_tree_list, "flow and energy");
      break;
    case 16: 
      pk_master_["multiphase"] = true;
      PopulatePKTree_(pk_tree_list, "multiphase");
      break;
    default:
      Exceptions::amanzi_throw(Errors::Message("This model is not supported by the MPC."));
    }

    Teuchos::ParameterList& tmp_list = out_list.sublist("time periods").sublist(ss.str());
    tmp_list.sublist("PK tree") = pk_tree_list;
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
  out_list.sublist("verbose object") = verb_list_.sublist("verbose object");

  return out_list;
}


/* ******************************************************************
* Resursive generation of PK tree.
****************************************************************** */
void InputConverterU::PopulatePKTree_(
    Teuchos::ParameterList& pk_tree, const std::string pk_name)
{
  std::string submodel, implicit("");
  if (transport_implicit_) implicit = " implicit";

  if (pk_name == "flow") {
    pk_tree.sublist("flow").set<std::string>("PK type", pk_model_["flow"]);    
  }
  else if (pk_name == "chemistry") {
    pk_tree.sublist("chemistry").set<std::string>("PK type", "chemistry" + implicit);
  }
  else if (pk_name == "transport") {
    pk_tree.sublist("transport").set<std::string>("PK type", "transport" + implicit);
  }
  else if (pk_name == "energy") {
    pk_tree.sublist("energy").set<std::string>("PK type", pk_model_["energy"]);
  }
  else if (pk_name == "multiphase") {
    pk_tree.sublist("multiphase").set<std::string>("PK type", pk_model_["multiphase"]);    
  }
  else if (pk_name == "coupled flow") {
    Teuchos::ParameterList& tmp_list = pk_tree.sublist("coupled flow");
    tmp_list.set<std::string>("PK type", "darcy matrix fracture");
    tmp_list.sublist("flow matrix").set<std::string>("PK type", "darcy");
    tmp_list.sublist("flow fracture").set<std::string>("PK type", "darcy");
  }
  else if (pk_name == "coupled transport") {
    Teuchos::ParameterList& tmp_list = pk_tree.sublist("coupled transport");
    tmp_list.set<std::string>("PK type", "transport matrix fracture" + implicit);
    tmp_list.sublist("transport matrix").set<std::string>("PK type", "transport" + implicit);
    tmp_list.sublist("transport fracture").set<std::string>("PK type", "transport" + implicit);
  }
  else if (pk_name == "coupled chemistry") {
    submodel = (pk_model_["chemistry"] == "amanzi") ? "chemistry amanzi" : "chemistry alquimia";
    Teuchos::ParameterList& tmp_list = pk_tree.sublist("coupled chemistry");
    tmp_list.set<std::string>("PK type", "chemistry matrix fracture");
    tmp_list.sublist("chemistry matrix").set<std::string>("PK type", submodel);
    tmp_list.sublist("chemistry fracture").set<std::string>("PK type", submodel);
  }
  else if (pk_name == "reactive transport") {
    submodel = (pk_model_["chemistry"] == "amanzi") ? "chemistry amanzi" : "chemistry alquimia";
    Teuchos::ParameterList& tmp_list = pk_tree.sublist("reactive transport");
    tmp_list.set<std::string>("PK type", "reactive transport");
    tmp_list.sublist("transport").set<std::string>("PK type", "transport");
    tmp_list.sublist("chemistry").set<std::string>("PK type", submodel);  
  }
  else if (pk_name == "coupled reactive transport") {
    Teuchos::ParameterList& tmp_list = pk_tree.sublist("coupled reactive transport");
    tmp_list.set<std::string>("PK type", "reactive transport matrix fracture");
    PopulatePKTree_(tmp_list, "coupled chemistry");
    PopulatePKTree_(tmp_list, "coupled transport");
  }
  else if (pk_name == "flow and chemistry") {
    submodel = (pk_model_["chemistry"] == "amanzi") ? "chemistry amanzi" : "chemistry alquimia";
    Teuchos::ParameterList& tmp_list = pk_tree.sublist("flow and chemistry");
    tmp_list.set<std::string>("PK type", "flow reactive transport");
    tmp_list.sublist("flow").set<std::string>("PK type", pk_model_["flow"]); 
    tmp_list.sublist("chemistry").set<std::string>("PK type", submodel);
  }
  else if (pk_name == "flow and transport") {
    Teuchos::ParameterList& tmp_list = pk_tree.sublist("flow and transport");
    tmp_list.set<std::string>("PK type", "flow reactive transport");
    tmp_list.sublist("flow").set<std::string>("PK type", pk_model_["flow"]);
    tmp_list.sublist("transport").set<std::string>("PK type", "transport");
  }
  else if (pk_name == "flow and energy") {
    Teuchos::ParameterList& tmp_list = pk_tree.sublist("flow and energy");
    tmp_list.set<std::string>("PK type", "thermal richards");
    tmp_list.sublist("flow").set<std::string>("PK type", pk_model_["flow"]);
    tmp_list.sublist("energy").set<std::string>("PK type", pk_model_["energy"]);
  }
  else if (pk_name == "coupled flow and transport") {
    Teuchos::ParameterList& tmp_list = pk_tree.sublist("coupled flow and transport");
    tmp_list.set<std::string>("PK type", "flow reactive transport");
    PopulatePKTree_(tmp_list, "coupled flow");
    PopulatePKTree_(tmp_list, "coupled transport");
  }
  else if (pk_name == "flow and reactive transport") {
    submodel = (pk_model_["chemistry"] == "amanzi") ? "chemistry amanzi" : "chemistry alquimia";
    Teuchos::ParameterList& tmp_list = pk_tree.sublist("flow and reactive transport");
    tmp_list.set<std::string>("PK type", "flow reactive transport");
    tmp_list.sublist("reactive transport").set<std::string>("PK type", "reactive transport");
    tmp_list.sublist("reactive transport").sublist("transport").set<std::string>("PK type", "transport");
    tmp_list.sublist("reactive transport").sublist("chemistry").set<std::string>("PK type", submodel);
    tmp_list.sublist("flow").set<std::string>("PK type", pk_model_["flow"]);
  }
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

  // get the default time steps
  bool flag;
  node = GetUniqueElementByTagsString_("execution_controls, execution_control_defaults", flag);

  double t, dt_init_d, dt_max_d;
  std::map<double, double> init_dt, max_dt;

  dt_init_d = GetAttributeValueD_(node, "init_dt", TYPE_TIME, DVAL_MIN, DVAL_MAX, "s", false, RESTART_TIMESTEP);
  dt_max_d = GetAttributeValueD_(node, "max_dt", TYPE_TIME, DVAL_MIN, DVAL_MAX, "s", false, MAXIMUM_TIMESTEP);

  children = doc_->getElementsByTagName(mm.transcode("execution_control"));
  for (int i = 0; i < children->getLength(); ++i) {
    DOMNode* inode = children->item(i);
    if (inode->getNodeType() != DOMNode::ELEMENT_NODE) continue;

    t = GetAttributeValueD_(inode, "start", TYPE_TIME, DVAL_MIN, DVAL_MAX, "s");
    double dt_init = GetAttributeValueD_(inode, "init_dt", TYPE_TIME, DVAL_MIN, DVAL_MAX, "s", false, dt_init_d);
    double dt_max = GetAttributeValueD_(inode, "max_dt", TYPE_TIME, DVAL_MIN, DVAL_MAX, "s", false, dt_max_d);
    init_dt[t] = dt_init;
    max_dt[t] = dt_max;
  }

  // add start times of all boundary conditions to the list
  std::map<double, double> dt_init_map, dt_max_map;

  std::vector<std::string> bc_names = {
      "uniform_pressure", "linear_pressure",
      "hydrostatic", "linear_hydrostatic",
      "inward_mass_flux", "outward_mass_flux",
      "inward_volumetric_flux", "outward_volumetric_flux",
      "seepage_face", "aqueous_conc",
      "uniform_conc", "constraint",
      "uniform_temperature", "no_flow"};

  node_list = doc_->getElementsByTagName(mm.transcode("boundary_conditions"));
  if (node_list->getLength() > 0) {
    node = node_list->item(0);

    for (int n = 0; n < bc_names.size(); ++n) {
      children = static_cast<DOMElement*>(node)->getElementsByTagName(mm.transcode(bc_names[n].c_str()));
      int nchildren = children->getLength();

      for (int i = 0; i < nchildren; ++i) {
        DOMNode* inode = children->item(i);
        if (inode->getNodeType() != DOMNode::ELEMENT_NODE) continue;

        t = GetAttributeValueD_(inode, "start", TYPE_TIME, DVAL_MIN, DVAL_MAX, "s");
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
  src_names.push_back("perm_weighted");
  src_names.push_back("volume_weighted");
  src_names.push_back("uniform_conc");
  src_names.push_back("flow_weighted_conc");
  src_names.push_back("flow_mass_fraction_conc");
  src_names.push_back("diffusion_dominated_release");

  node_list = doc_->getElementsByTagName(mm.transcode("sources"));
  if (node_list->getLength() > 0) {
    node = node_list->item(0);

    for (int n = 0; n < src_names.size(); ++n) {
      children = static_cast<DOMElement*>(node)->getElementsByTagName(mm.transcode(src_names[n].c_str()));
      for (int i = 0; i < children->getLength(); ++i) {
        DOMNode* inode = children->item(i);
        if (inode->getNodeType() != DOMNode::ELEMENT_NODE) continue;

        t = GetAttributeValueD_(inode, "start", TYPE_TIME, DVAL_MIN, DVAL_MAX, "s");
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
Teuchos::ParameterList InputConverterU::TranslatePKs_(Teuchos::ParameterList& glist)
{
  Teuchos::ParameterList out_list;

  if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH)
      *vo_->os() << "Translating process kernels" << std::endl;
  Teuchos::OSTab tab = vo_->getOSTab();

  bool flag;
  DOMNode* node;

  // create PKs list
  Teuchos::ParameterList tp_list = glist.sublist("cycle driver").sublist("time periods");

  for (auto it = tp_list.begin(); it != tp_list.end(); ++it) {
    if ((it->second).isList()) {
      Teuchos::ParameterList& pk_tree = tp_list.sublist(it->first).sublist("PK tree");
      RegisterPKsList_(pk_tree, out_list);

      if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH) {
        std::string name = pk_tree.begin()->first;
        *vo_->os() << "TP: PK name=\"" << name << "\", factory: \"" 
                   << pk_tree.sublist(name).get<std::string>("PK type") << "\"" << std::endl;
      }
    }
  }

  // parse list of supported PKs
  for (auto it = out_list.begin(); it != out_list.end(); ++it) {
    if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH)
        *vo_->os() << "Next PK is \"" << it->first << "\"" << std::endl;

    if ((it->second).isList()) {
      // -- flow PKs
      if (it->first == "flow steady" || it->first == "flow matrix steady") {
        out_list.sublist(it->first) = TranslateFlow_("steady", "matrix");
      }
      if (it->first == "flow fracture steady") {
        out_list.sublist(it->first) = TranslateFlow_("steady", "fracture");
      }
      else if (it->first == "flow" || it->first == "flow matrix") {
        out_list.sublist(it->first) = TranslateFlow_("transient", "matrix");
      }
      else if (it->first == "flow fracture") {
        out_list.sublist(it->first) = TranslateFlow_("transient", "fracture");
      }
      // -- energy PKs
      else if (it->first == "energy") {
        out_list.sublist(it->first) = TranslateEnergy_();
      }
      // -- transport PKs
      else if (it->first == "transport" || it->first == "transport matrix") {
        out_list.sublist(it->first) = TranslateTransport_("matrix");
      }
      else if (it->first == "transport fracture") {
        out_list.sublist(it->first) = TranslateTransport_("fracture");
      }
      // -- chemistry PKs
      else if (it->first == "chemistry" || it->first == "chemistry matrix") {
        out_list.sublist(it->first) = TranslateChemistry_("matrix");
      }
      else if (it->first == "chemistry fracture") {
        out_list.sublist(it->first) = TranslateChemistry_("fracture");
      }
      // -- multiphase PKs
      else if (it->first == "multiphase") {
        auto& tmp = glist.sublist("state");
        out_list.sublist(it->first) = TranslateMultiphase_("matrix", tmp);
      }
      // -- coupled PKs
      else if (it->first == "coupled flow") {
        Teuchos::Array<std::string> pk_names;
        pk_names.push_back("flow matrix");
        pk_names.push_back("flow fracture");
        out_list.sublist(it->first).set<Teuchos::Array<std::string> >("PKs order", pk_names);
        out_list.sublist(it->first).set<int>("master PK index", 0);
      }
      else if (it->first == "coupled flow steady") {
        Teuchos::Array<std::string> pk_names;
        pk_names.push_back("flow matrix steady");
        pk_names.push_back("flow fracture steady");
        out_list.sublist(it->first).set<Teuchos::Array<std::string> >("PKs order", pk_names);
        out_list.sublist(it->first).set<int>("master PK index", 0);
      }
      else if (it->first == "coupled transport") {
        Teuchos::Array<std::string> pk_names;
        pk_names.push_back("transport matrix");
        pk_names.push_back("transport fracture");
        out_list.sublist(it->first).set<Teuchos::Array<std::string> >("PKs order", pk_names);
        out_list.sublist(it->first).set<int>("master PK index", 0);

        // implicit PK derived from a strongly coupled MPC needs a time integrator 
        if (transport_implicit_) {
          std::string err_options("residual"), mode("transient");
          std::string nonlinear_solver("nka"), tags_default("unstructured_controls, unstr_nonlinear_solver");

          node = GetUniqueElementByTagsString_(tags_default, flag);
          if (flag) nonlinear_solver = GetAttributeValueS_(node, "name", TYPE_NONE, false, "nka"); 

          out_list.sublist(it->first).sublist("time integrator") = TranslateTimeIntegrator_(
              err_options, nonlinear_solver, false, tags_default,
              dt_cut_[mode], dt_inc_[mode]);
        }
      }
      else if (it->first == "coupled chemistry") {
        Teuchos::Array<std::string> pk_names;
        pk_names.push_back("chemistry matrix");
        pk_names.push_back("chemistry fracture");
        out_list.sublist(it->first).set<Teuchos::Array<std::string> >("PKs order", pk_names);
        out_list.sublist(it->first).set<int>("master PK index", 0);
      }
      else if (it->first == "reactive transport") {
        Teuchos::Array<std::string> pk_names;
        pk_names.push_back("chemistry");
        pk_names.push_back("transport");
        out_list.sublist(it->first).set<Teuchos::Array<std::string> >("PKs order", pk_names);
      }
      else if (it->first == "flow and reactive transport") {
        Teuchos::Array<std::string> pk_names;
        pk_names.push_back("flow");
        pk_names.push_back("reactive transport");
        out_list.sublist(it->first).set<Teuchos::Array<std::string> >("PKs order", pk_names);
        out_list.sublist(it->first).set<int>("master PK index", 0);
      }
      else if (it->first == "flow and transport") {
        Teuchos::Array<std::string> pk_names;
        pk_names.push_back("flow");
        pk_names.push_back("transport");
        out_list.sublist(it->first).set<Teuchos::Array<std::string> >("PKs order", pk_names);
        out_list.sublist(it->first).set<int>("master PK index", 0);
      }
      else if (it->first == "flow and chemistry") {
        Teuchos::Array<std::string> pk_names;
        pk_names.push_back("flow");
        pk_names.push_back("chemistry");
        out_list.sublist(it->first).set<Teuchos::Array<std::string> >("PKs order", pk_names);
        out_list.sublist(it->first).set<int>("master PK index", 0);
      }
      else if (it->first == "flow and energy") {
        Teuchos::Array<std::string> pk_names;
        pk_names.push_back("flow");
        pk_names.push_back("energy");
        out_list.sublist(it->first).set<Teuchos::Array<std::string> >("PKs order", pk_names);
        out_list.sublist(it->first).set<int>("master PK index", 0);

        if (pk_master_.find("thermal richards") != pk_master_.end()) {
          // we use steady defaults so far
          out_list.sublist(it->first).sublist("time integrator") = TranslateTimeIntegrator_(
              "pressure, temperature", "nka", false,
              "unstructured_controls, unstr_thermal_richards_controls",
              TI_TS_REDUCTION_FACTOR, TI_TS_INCREASE_FACTOR);  
          out_list.sublist(it->first).sublist("verbose object") = verb_list_.sublist("verbose object");
        }
      }
      else if (it->first == "coupled flow and transport") {
        Teuchos::Array<std::string> pk_names;
        pk_names.push_back("coupled flow");
        pk_names.push_back("coupled transport");
        out_list.sublist(it->first).set<Teuchos::Array<std::string> >("PKs order", pk_names);
        out_list.sublist(it->first).set<int>("master PK index", 0);
      }
      else if (it->first == "coupled reactive transport") {
        Teuchos::Array<std::string> pk_names;
        pk_names.push_back("coupled chemistry");
        pk_names.push_back("coupled transport");
        out_list.sublist(it->first).set<Teuchos::Array<std::string> >("PKs order", pk_names);
        out_list.sublist(it->first).set<int>("master PK index", 0);
      }
    }
  }

  return out_list;
}


/* ******************************************************************
* Global changes to finalize MPC PKs
****************************************************************** */
void InputConverterU::FinalizeMPC_PKs_(Teuchos::ParameterList& glist)
{
  auto& pk_list = glist.sublist("PKs");
  auto& mesh_list = glist.sublist("mesh").sublist("unstructured");

  for (auto it = pk_list.begin(); it != pk_list.end(); ++it) {
    std::string name = it->first;
    if (name == "coupled flow" || name == "coupled flow steady") {
      std::string flow_m("flow matrix"), flow_f("flow fracture");
      if (name == "coupled flow steady") {
        flow_m += " steady";
        flow_f += " steady";
      }
      pk_list.sublist(name).sublist("time integrator") = 
          pk_list.sublist(flow_m).sublist("time integrator");
      
      pk_list.sublist(name).sublist("time integrator").sublist("BDF1").sublist("nka parameters")
         .set<std::string>("monitor", "monitor l2 residual");

      auto& tmp_m = pk_list.sublist(flow_m).sublist("time integrator");
      tmp_m.set<std::string>("time integration method", "none");
      tmp_m.remove("BDF1", false);
      tmp_m.remove("initialization", false);
      pk_list.sublist(flow_m).sublist("operators")
         .sublist("diffusion operator").sublist("matrix").set<Teuchos::Array<std::string> >("fracture", fracture_regions_);

      auto& tmp_f = pk_list.sublist(flow_f).sublist("time integrator");
      tmp_f.set<std::string>("time integration method", "none");
      tmp_f.remove("BDF1", false);
      tmp_f.remove("initialization", false);

      // Teuchos::Array<std::string> aux(1, CreateUniqueName_(fracture_regions_));
      Teuchos::Array<std::string> aux(1, "FRACTURE_NETWORK_INTERNAL");
      mesh_list.sublist("submesh").set<Teuchos::Array<std::string> >("regions", aux)
                                  .set<std::string>("extraction method", "manifold mesh");

      if (dim_ == 3) mesh_list.sublist("expert").set<bool>("request edges", true);
    }

    if (name == "coupled transport" && transport_implicit_) {
      pk_list.sublist("transport matrix").sublist("operators").sublist("advection operator")
             .sublist("matrix").set<Teuchos::Array<std::string> >("fracture", fracture_regions_);
    }
  }
}


/* ******************************************************************
* Empty
****************************************************************** */
void InputConverterU::RegisterPKsList_(
    Teuchos::ParameterList& pk_tree, Teuchos::ParameterList& pks_list)
{
  for (auto it = pk_tree.begin(); it !=pk_tree.end();++it) {
    if ((it->second).isList()) {
      pks_list.sublist(it->first);
      RegisterPKsList_(pk_tree.sublist(it->first), pks_list);
    }   
  }
}

}  // namespace AmanziInput
}  // namespace Amanzi
