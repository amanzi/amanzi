/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Daniil Svyatskiy (original version)
      Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Input Converter

*/

#include <algorithm>
#include <climits>
#include <set>
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
Teuchos::ParameterList
InputConverterU::TranslateCycleDriver_()
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
  char* tagname;
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

    dt_cut_[mode] = ConvertUnits_(
      GetAttributeValueS_(inode, "reduction_factor", TYPE_TIME, false, dt_cut_d), unit);
    dt_inc_[mode] = ConvertUnits_(
      GetAttributeValueS_(inode, "increase_factor", TYPE_TIME, false, dt_inc_d), unit);

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
    if (init_filename_.size() == 0)
      ThrowErrorIllformed_("execution_controls", "initialize", "filename");
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
      pk_model_["transport"] = "transport";
      pk_domain_["transport"] = "domain";
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
      node = GetUniqueElementByTagsString_("numerical_controls, unstructured_controls, "
                                           "unstr_steady-state_controls, unstr_initialization",
                                           flag);
      if (!flag) {
        msg << "Constant flow must have unstr_steady-state_controls->unstr_initialization list, "
               "unless state=off.\n";
        Exceptions::amanzi_throw(msg);
      }
    }

    Teuchos::ParameterList& tmp_list = out_list.sublist("time periods").sublist("TP 0");
    if (!coupled_flow_) {
      tmp_list.sublist("PK tree")
        .sublist("steady:flow")
        .set<std::string>("PK type", pk_model_["flow"]);
    } else {
      Teuchos::ParameterList& aux_list = tmp_list.sublist("PK tree").sublist("steady:coupled flow");
      aux_list.set<std::string>("PK type", "darcy matrix fracture");
      aux_list.sublist("steady:flow matrix").set<std::string>("PK type", "darcy");
      aux_list.sublist("steady:flow fracture").set<std::string>("PK type", "darcy");
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
    transport_implicit_ = (algorithm == "implicit" || algorithm == "implicit second-order");
  }

  std::string submodel;
  std::map<double, std::string>::iterator it = tp_mode.begin();
  while (it != tp_mode.end()) {
    switch (transient_model) {
    case 1:
      PopulatePKTree_(pk_tree_list, "transient:chemistry");
      break;
    case 2:
      if (!coupled_transport_)
        PopulatePKTree_(pk_tree_list, "transient:transport");
      else
        PopulatePKTree_(pk_tree_list, "transient:coupled transport");
      break;
    case 3:
      if (!coupled_transport_)
        PopulatePKTree_(pk_tree_list, "transient:reactive transport");
      else
        PopulatePKTree_(pk_tree_list, "transient:coupled reactive transport");
      break;
    case 4:
      if (!coupled_flow_)
        PopulatePKTree_(pk_tree_list, "transient:flow");
      else
        PopulatePKTree_(pk_tree_list, "transient:coupled flow");
      break;
    case 5:
      PopulatePKTree_(pk_tree_list, "transient:flow and chemistry");
      break;
    case 6:
      if (!coupled_flow_)
        PopulatePKTree_(pk_tree_list, "transient:flow and transport");
      else
        PopulatePKTree_(pk_tree_list, "transient:coupled flow and transport");
      break;
    case 7:
      if (!coupled_flow_)
        PopulatePKTree_(pk_tree_list, "transient:flow and reactive transport");
      else
        PopulatePKTree_(pk_tree_list, "transient:coupled flow and reactive transport");
      break;
    default:
      msg << "The model with id=" << transient_model << " is not supported by the MPC.\n";
      Exceptions::amanzi_throw(msg);
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
    out_list.set<Teuchos::Array<std::string>>("component names", comp_names_all_);
    out_list.set<int>("number of liquid components", phases_["water"].size());

    // available molar masses
    std::string name;
    std::vector<double> tmp(comp_names_all_.size(), 1.0);

    for (int i = 0; i < comp_names_all_.size(); ++i) {
      name = comp_names_all_[i];
      if (solute_molar_mass_.find(name) != solute_molar_mass_.end())
        tmp[i] = solute_molar_mass_[name];
    }
    out_list.set<Teuchos::Array<double>>("component molar masses", tmp);
  }

  out_list.sublist("time period control") = TranslateTimePeriodControls_();
  if (filename.size() > 0) {
    restart_ = true;
    out_list.sublist("restart").set<std::string>("file name", filename);
  }
  out_list.sublist("verbose object") = verb_list_.sublist("verbose object");

  return out_list;
}

} // namespace AmanziInput
} // namespace Amanzi
