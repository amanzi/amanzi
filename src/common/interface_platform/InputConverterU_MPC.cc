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

static const char delimiter = ':';

XERCES_CPP_NAMESPACE_USE

/* ******************************************************************
* Create new cycle driver list for this type lists:
<process_kernels>
  <pk mode="steady">
    <flow model="richards" state="on"/>
  </pk>
  <pk mode="transient1">
    <flow model="richards" state="on"/>
    <energy model="two-phase" state="on"/>
    <stronly_coupled name="mpc">flow,energy</strongly_coupled>
  </pk>
  <pk mode="transient2">
    <flow model="richards" state="on"/>
    <energy model="two-phase" state="on"/>
    <transport state="on"/>
    <stronly_coupled name="mpc">flow,energy</strongly_coupled>
    <sub_cycling>mpc,transport</sub_cycling>
  </pk>
  <pk mode="transient3">
    <flow model="richards" state="on"/>
    <energy model="two-phase" state="on"/>
    <transport state="on"/>
    <chemistry engine="none" process_model="none" state="on"/>
    <stronly_coupled name="mpc">flow,energy</strongly_coupled>
    <sub_cycling>mpc,transport,chemistry</sub_cycling>
  </pk>
  <pk mode="transient4">
    <flow model="richards" state="on"/>
    <energy model="two-phase" state="on"/>
    <transport state="on"/>
    <chemistry engine="none" process_model="none" state="on"/>
    <stronly_coupled name="mpc1">flow,energy</strongly_coupled>
    <weakly_coupled name="mpc2">transport,chemistry</weakly_coupled>
    <sub_cycling>mpc1,mpc2</sub_cycling>
  </pk>
</process_kernels>
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

  // global transport properties (should it be here?)
  node = GetUniqueElementByTagsString_(
      "numerical_controls, unstructured_controls, unstr_transport_controls, algorithm", flag);
  if (flag) {
    std::string algorithm = TrimString_(mm.transcode(node->getTextContent()));
    transport_implicit_ = (algorithm == "implicit");
  }

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

  GetUniqueElementByTagsString_("execution_controls", flag);
  GetUniqueElementByTagsString_("fracture_network", coupled_flow_);
  coupled_transport_ = coupled_flow_;
  coupled_energy_ = coupled_flow_;

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
    pk_domain_.clear();
    pk_region_.clear();

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
        pk_domain_["transport"] = "domain";

      } else if (strcmp(tagname, "energy") == 0) {
        model = GetAttributeValueS_(jnode, "model");
        pk_model_["energy"] = model;
        pk_master_["energy"] = true;
        GetAttributeValueS_(jnode, "state", "on");
        transient_model += 8;

      } else if (strcmp(tagname, "shallow_water") == 0) {
        model = GetAttributeValueS_(jnode, "model");
        pk_model_["shallow_water"] = model;

        std::string region = GetAttributeValueS_(jnode, "domain");
        surface_regions_.push_back(region);
        pk_region_["shallow_water"] = region;
        pk_domain_["shallow_water"] = "domain";

        GetAttributeValueS_(jnode, "state", "on");
        transient_model += 16;

      } else if (strcmp(tagname, "multiphase") == 0) {
        GetAttributeValueS_(jnode, "state", "on");
        model = GetAttributeValueS_(jnode, "model");
        pk_model_["multiphase"] = model;
        pk_master_["multiphase"] = true;
        transient_model += 32;
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
      PopulatePKTree_(pk_tree_list, Keys::merge(mode, "chemistry", delimiter));
      break;
    case 2:
      if (!coupled_transport_)
        PopulatePKTree_(pk_tree_list, Keys::merge(mode, "transport", delimiter));
      else
        PopulatePKTree_(pk_tree_list, Keys::merge(mode, "coupled transport", delimiter));
      break;
    case 3:
      if (!coupled_transport_)
        PopulatePKTree_(pk_tree_list, Keys::merge(mode, "reactive transport", delimiter));
      else
        PopulatePKTree_(pk_tree_list, Keys::merge(mode, "coupled reactive transport", delimiter));
      break;
    case 4:
      if (!coupled_flow_)
        PopulatePKTree_(pk_tree_list, Keys::merge(mode, "flow", delimiter));
      else
        PopulatePKTree_(pk_tree_list, Keys::merge(mode, "coupled flow", delimiter));
      break;
    case 5:
      PopulatePKTree_(pk_tree_list, Keys::merge(mode, "flow and chemistry", delimiter));
      break;
    case 6:
      if (!coupled_flow_)
        PopulatePKTree_(pk_tree_list, Keys::merge(mode, "flow and transport", delimiter));
      else
        PopulatePKTree_(pk_tree_list, Keys::merge(mode, "coupled flow and transport", delimiter));
      break;
    case 7:
      if (!coupled_flow_)
        PopulatePKTree_(pk_tree_list, Keys::merge(mode, "flow and reactive transport", delimiter));
      else 
        PopulatePKTree_(pk_tree_list, Keys::merge(mode, "coupled flow and reactive transport", delimiter));
      break;
    case 8:
      PopulatePKTree_(pk_tree_list, Keys::merge(mode, "energy", delimiter));
      break;
    case 12: 
      pk_master_["thermal flow"] = true;
      if (!coupled_flow_)
        PopulatePKTree_(pk_tree_list, Keys::merge(mode, "flow and energy", delimiter));
      else 
        PopulatePKTree_(pk_tree_list, Keys::merge(mode, "coupled flow and energy", delimiter));
      break;
    case 15: 
      pk_master_["thermal flow"] = true;
      if (!coupled_flow_)
        PopulatePKTree_(pk_tree_list, Keys::merge(mode, "thermal flow and reactive transport", delimiter));
      else 
        PopulatePKTree_(pk_tree_list, Keys::merge(mode, "coupled thermal flow and reactive transport", delimiter));
      break;
    case 16:
      PopulatePKTree_(pk_tree_list, Keys::merge(mode, "shallow water", delimiter));
      break;
    case 18:
      PopulatePKTree_(pk_tree_list, Keys::merge(mode, "shallow water and transport", delimiter));
      pk_domain_["shallow_water"] = "domain";
      pk_domain_["transport"] = "domain";
      break;
    case 20:
      PopulatePKTree_(pk_tree_list, Keys::merge(mode, "flow and shallow water", delimiter));
      pk_domain_["shallow_water"] = "surface";  // for integrated hydrology, SW is defined on manifold
      break;
    case 32: 
      pk_master_["multiphase"] = true;
      PopulatePKTree_(pk_tree_list, Keys::merge(mode, "multiphase", delimiter));
      break;
    default:
      Errors::Message msg;
      msg << "The model with id=" << transient_model << " is not supported by the MPC.\n";
      Exceptions::amanzi_throw(msg);
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

  Key prefix = Keys::split(pk_name, delimiter).first;
  Key basename = Keys::split(pk_name, delimiter).second;

  Teuchos::ParameterList& tmp = pk_tree.sublist(pk_name);

  if (basename == "flow") {
    tmp.set<std::string>("PK type", pk_model_["flow"]);    
  }
  else if (basename == "chemistry") {
    tmp.set<std::string>("PK type", "chemistry" + implicit);
  }
  else if (basename == "transport") {
    tmp.set<std::string>("PK type", "transport" + implicit);
  }
  else if (basename == "energy") {
    tmp.set<std::string>("PK type", pk_model_["energy"]);
  }
  else if (basename == "shallow water") {
    tmp.set<std::string>("PK type", "shallow water");
  }
  else if (basename == "multiphase") {
    tmp.set<std::string>("PK type", pk_model_["multiphase"]);    
  }
  else if (basename == "coupled flow") {
    tmp.set<std::string>("PK type", "darcy matrix fracture");
    tmp.sublist(Keys::merge(prefix, "flow matrix", delimiter))
       .set<std::string>("PK type", pk_model_["flow"]);
    tmp.sublist(Keys::merge(prefix, "flow fracture", delimiter))
       .set<std::string>("PK type", pk_model_["flow"]);
  }
  else if (basename == "coupled transport") {
    use_transport_dispersion_ = false;
    tmp.set<std::string>("PK type", "transport matrix fracture" + implicit);
    tmp.sublist(Keys::merge(prefix, "transport matrix", delimiter)).set<std::string>("PK type", "transport" + implicit);
    tmp.sublist(Keys::merge(prefix, "transport fracture", delimiter)).set<std::string>("PK type", "transport" + implicit);
  }
  else if (basename == "coupled chemistry") {
    submodel = (pk_model_["chemistry"] == "amanzi") ? "chemistry amanzi" : "chemistry alquimia";
    tmp.set<std::string>("PK type", "chemistry matrix fracture");
    tmp.sublist(Keys::merge(prefix, "chemistry matrix", delimiter)).set<std::string>("PK type", submodel);
    tmp.sublist(Keys::merge(prefix, "chemistry fracture", delimiter)).set<std::string>("PK type", submodel);
  }
  else if (basename == "reactive transport") {
    submodel = (pk_model_["chemistry"] == "amanzi") ? "chemistry amanzi" : "chemistry alquimia";
    tmp.set<std::string>("PK type", "reactive transport");
    tmp.sublist(Keys::merge(prefix, "transport", delimiter)).set<std::string>("PK type", "transport");
    tmp.sublist(Keys::merge(prefix, "chemistry", delimiter)).set<std::string>("PK type", submodel);  
  }
  else if (basename == "coupled reactive transport") {
    tmp.set<std::string>("PK type", "reactive transport matrix fracture");
    PopulatePKTree_(tmp, Keys::merge(prefix, "coupled chemistry", delimiter));
    PopulatePKTree_(tmp, Keys::merge(prefix, "coupled transport", delimiter));
  }
  else if (basename == "flow and chemistry") {
    submodel = (pk_model_["chemistry"] == "amanzi") ? "chemistry amanzi" : "chemistry alquimia";
    tmp.set<std::string>("PK type", "flow reactive transport");
    tmp.sublist(Keys::merge(prefix, "flow", delimiter)).set<std::string>("PK type", pk_model_["flow"]); 
    tmp.sublist(Keys::merge(prefix, "chemistry", delimiter)).set<std::string>("PK type", submodel);
  }
  else if (basename == "flow and transport") {
    tmp.set<std::string>("PK type", "flow reactive transport");
    tmp.sublist(Keys::merge(prefix, "flow", delimiter)).set<std::string>("PK type", pk_model_["flow"]);
    tmp.sublist(Keys::merge(prefix, "transport", delimiter)).set<std::string>("PK type", "transport");
  }
  else if (basename == "flow and energy") {
    tmp.set<std::string>("PK type", "thermal flow");
    tmp.sublist(Keys::merge(prefix, "flow", delimiter)).set<std::string>("PK type", pk_model_["flow"]);
    tmp.sublist(Keys::merge(prefix, "energy", delimiter)).set<std::string>("PK type", pk_model_["energy"]);
  }
  else if (basename == "flow and shallow water") {
    tmp.set<std::string>("PK type", "surface subsurface");
    tmp.sublist(Keys::merge(prefix, "flow", delimiter)).set<std::string>("PK type", pk_model_["flow"]);
    tmp.sublist(Keys::merge(prefix, "shallow water", delimiter)).set<std::string>("PK type", "shallow water");
  }
  else if (basename == "flow and energy fracture") {
    tmp.set<std::string>("PK type", "thermal flow");
    tmp.sublist(Keys::merge(prefix, "flow fracture", delimiter)).set<std::string>("PK type", pk_model_["flow"]);
    tmp.sublist(Keys::merge(prefix, "energy fracture", delimiter)).set<std::string>("PK type", pk_model_["energy"]);
  }
  else if (basename == "coupled flow and energy") {
    tmp.set<std::string>("PK type", "thermal flow matrix fracture");
    PopulatePKTree_(tmp, Keys::merge(prefix, "flow and energy", delimiter));
    PopulatePKTree_(tmp, Keys::merge(prefix, "flow and energy fracture", delimiter));
  }
  else if (basename == "coupled flow and transport") {
    tmp.set<std::string>("PK type", "flow reactive transport");
    PopulatePKTree_(tmp, Keys::merge(prefix, "coupled flow", delimiter));
    PopulatePKTree_(tmp, Keys::merge(prefix, "coupled transport", delimiter));
  }
  else if (basename == "flow and reactive transport") {
    submodel = (pk_model_["chemistry"] == "amanzi") ? "chemistry amanzi" : "chemistry alquimia";
    tmp.set<std::string>("PK type", "flow reactive transport");
    tmp.sublist(Keys::merge(prefix, "reactive transport", delimiter))
       .set<std::string>("PK type", "reactive transport");
    tmp.sublist(Keys::merge(prefix, "reactive transport", delimiter))
       .sublist(Keys::merge(prefix, "transport", delimiter)).set<std::string>("PK type", "transport");
    tmp.sublist(Keys::merge(prefix, "reactive transport", delimiter))
       .sublist(Keys::merge(prefix, "chemistry", delimiter)).set<std::string>("PK type", submodel);
    tmp.sublist(Keys::merge(prefix, "flow", delimiter)).set<std::string>("PK type", pk_model_["flow"]);
  }
  else if (basename == "coupled flow and reactive transport") {
    tmp.set<std::string>("PK type", "flow reactive transport");  // same as for single domain
    PopulatePKTree_(tmp, Keys::merge(prefix, "coupled flow", delimiter));
    PopulatePKTree_(tmp, Keys::merge(prefix, "coupled reactive transport", delimiter));
  }
  else if (basename == "coupled thermal flow and reactive transport") {
    tmp.set<std::string>("PK type", "flow reactive transport");
    PopulatePKTree_(tmp, Keys::merge(prefix, "coupled flow and energy", delimiter));
    PopulatePKTree_(tmp, Keys::merge(prefix, "coupled reactive transport", delimiter));
  }
  else if (basename == "shallow water and transport") {
    tmp.set<std::string>("PK type", "shallow water transport");
    tmp.sublist(Keys::merge(prefix, "shallow water", delimiter)).set<std::string>("PK type", "shallow water");
    tmp.sublist(Keys::merge(prefix, "transport", delimiter)).set<std::string>("PK type", "transport");
  }
  else {
    Errors::Message msg;
    msg << "Internal error: cannot add \"" << pk_name << "\" to the PK tree.\n";
    Exceptions::amanzi_throw(msg);
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

        t = GetAttributeValueD_(inode, "start", TYPE_TIME, DVAL_MIN, DVAL_MAX, "s", false, DVAL_MIN);
        if (t == DVAL_MIN) continue;

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

        t = GetAttributeValueD_(inode, "start", TYPE_TIME, DVAL_MIN, DVAL_MAX, "s", false, DVAL_MIN);
        if (t == DVAL_MIN) continue;

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
* Translate PKs list starting with lower-level PKs
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
        std::string pk = pk_tree.begin()->first;
        *vo_->os() << "TP: PK name=\"" << pk << "\", factory: \"" 
                   << pk_tree.sublist(pk).get<std::string>("PK type") << "\"" << std::endl;
      }
    }
  }

  // create list of PKs
  std::vector<std::string> pks;
  for (auto it = out_list.begin(); it != out_list.end(); ++it) {
    if ((it->second).isList()) pks.push_back(it->first);
  }

  // parse list of supported PKs in reverse order to allow MPC
  // contribute to its PKs.
  for (auto it = pks.crbegin(); it != pks.crend(); ++it) {
    std::string pk = *it;
    if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH)
        *vo_->os() << "Next PK is \"" << pk << "\"" << std::endl;

    std::string err_options("residual"), mode("transient");

    if (out_list.isSublist(pk)) {
      // -- flow PKs
      std::string prefix = Keys::split(pk, delimiter).first;
      std::string basename = Keys::split(pk, delimiter).second;

      if (prefix == "steady" && (basename == "flow" || basename == "flow matrix")) {
        out_list.sublist(pk) = TranslateFlow_("steady", "matrix");
      }
      else if (pk == "steady:flow fracture") {
        out_list.sublist(pk) = TranslateFlow_("steady", "fracture");
      }
      else if (prefix != "steady" && (basename == "flow" || basename == "flow matrix")) {
        out_list.sublist(pk) = TranslateFlow_("transient", "matrix");
      }
      else if (basename == "flow fracture") {
        out_list.sublist(pk) = TranslateFlow_("transient", "fracture");
      }
      // -- energy PKs
      else if (basename == "energy") {
        out_list.sublist(pk) = TranslateEnergy_("matrix");
      }
      else if (basename == "energy fracture") {
        out_list.sublist(pk) = TranslateEnergy_("fracture");
      }
      // -- transport PKs
      else if (basename == "transport") {
        out_list.sublist(pk) = TranslateTransport_(pk_domain_["transport"]);
      }
      else if (basename == "transport matrix") {
        out_list.sublist(pk) = TranslateTransport_("matrix");
      }
      else if (basename == "transport fracture") {
        out_list.sublist(pk) = TranslateTransport_("fracture");
      }
      // -- chemistry PKs
      else if (basename == "chemistry" || basename == "chemistry matrix") {
        out_list.sublist(pk) = TranslateChemistry_("matrix");
      }
      else if (basename == "chemistry fracture") {
        out_list.sublist(pk) = TranslateChemistry_("fracture");
      }
      // -- surface PKs
      else if (basename == "shallow water") {
        // temporarily, we run only stand-along SW
        out_list.sublist(pk) = TranslateShallowWater_(pk_domain_["shallow_water"]);
      }
      // -- multiphase PKs
      else if (basename == "multiphase") {
        auto& tmp = glist.sublist("state");
        out_list.sublist(pk) = TranslateMultiphase_("matrix", tmp);
      }
      // -- coupled PKs (matrix and fracture)
      else if (basename == "coupled flow") {
        Teuchos::Array<std::string> pk_names;
        pk_names.push_back(Keys::merge(prefix, "flow matrix", delimiter));
        pk_names.push_back(Keys::merge(prefix, "flow fracture", delimiter));
        out_list.sublist(pk).set<Teuchos::Array<std::string> >("PKs order", pk_names);
        out_list.sublist(pk).set<int>("master PK index", 0);
      }
      else if (pk == "steady:coupled flow") {
        Teuchos::Array<std::string> pk_names;
        pk_names.push_back("steady:flow matrix");
        pk_names.push_back("steady:flow fracture");
        out_list.sublist(pk).set<Teuchos::Array<std::string> >("PKs order", pk_names);
        out_list.sublist(pk).set<int>("master PK index", 0);
      }
      else if (basename == "coupled transport") {
        Teuchos::Array<std::string> pk_names;
        pk_names.push_back(Keys::merge(prefix, "transport matrix", delimiter));
        pk_names.push_back(Keys::merge(prefix, "transport fracture", delimiter));
        out_list.sublist(pk).set<Teuchos::Array<std::string> >("PKs order", pk_names);
        out_list.sublist(pk).set<int>("master PK index", 0);

        // implicit PK derived from a strongly coupled MPC needs a time integrator 
        if (transport_implicit_) {
          std::string nonlinear_solver("nka"), tags_default("unstructured_controls, unstr_nonlinear_solver");

          node = GetUniqueElementByTagsString_(tags_default, flag);
          if (flag) nonlinear_solver = GetAttributeValueS_(node, "name", TYPE_NONE, false, "nka"); 

          out_list.sublist(pk).sublist("time integrator") = TranslateTimeIntegrator_(
              err_options, nonlinear_solver, false, tags_default,
              dt_cut_[mode], dt_inc_[mode]);
        }
      }
      else if (basename == "coupled chemistry") {
        Teuchos::Array<std::string> pk_names;
        pk_names.push_back(Keys::merge(prefix, "chemistry matrix", delimiter));
        pk_names.push_back(Keys::merge(prefix, "chemistry fracture", delimiter));
        out_list.sublist(pk).set<Teuchos::Array<std::string> >("PKs order", pk_names);
        out_list.sublist(pk).set<int>("master PK index", 0);
      }
      // -- multiple PKs
      else if (basename == "reactive transport") {
        Teuchos::Array<std::string> pk_names;
        pk_names.push_back(Keys::merge(prefix, "chemistry", delimiter));
        pk_names.push_back(Keys::merge(prefix, "transport", delimiter));
        out_list.sublist(pk).set<Teuchos::Array<std::string> >("PKs order", pk_names);
      }
      else if (basename == "flow and reactive transport") {
        Teuchos::Array<std::string> pk_names;
        pk_names.push_back(Keys::merge(prefix, "flow", delimiter));
        pk_names.push_back(Keys::merge(prefix, "reactive transport", delimiter));
        out_list.sublist(pk).set<Teuchos::Array<std::string> >("PKs order", pk_names);
        out_list.sublist(pk).set<int>("master PK index", 0);
      }
      else if (basename == "flow and transport") {
        Teuchos::Array<std::string> pk_names;
        pk_names.push_back(Keys::merge(prefix, "flow", delimiter));
        pk_names.push_back(Keys::merge(prefix, "transport", delimiter));
        out_list.sublist(pk).set<Teuchos::Array<std::string> >("PKs order", pk_names);
        out_list.sublist(pk).set<int>("master PK index", 0);
      }
      else if (basename == "flow and chemistry") {
        Teuchos::Array<std::string> pk_names;
        pk_names.push_back(Keys::merge(prefix, "flow", delimiter));
        pk_names.push_back(Keys::merge(prefix, "chemistry", delimiter));
        out_list.sublist(pk).set<Teuchos::Array<std::string> >("PKs order", pk_names);
        out_list.sublist(pk).set<int>("master PK index", 0);
      }
      else if (basename == "flow and energy") {
        Teuchos::Array<std::string> pk_names;
        pk_names.push_back(Keys::merge(prefix, "flow", delimiter));
        pk_names.push_back(Keys::merge(prefix, "energy", delimiter));
        out_list.sublist(pk)
            .set<Teuchos::Array<std::string> >("PKs order", pk_names)
            .set<int>("master PK index", 0)
            .set<std::string>("domain name", "domain");

        auto& tmp = glist.sublist("state").sublist("evaluators");
        AddSecondaryFieldEvaluator_(tmp, 
            Keys::getKey("domain", "molar_density_liquid"), "molar density key",
            "eos", "liquid water 0-30C", "density");
        AddSecondaryFieldEvaluator_(tmp, 
            Keys::getKey("domain", "viscosity_liquid"), "viscosity key",
            "viscosity", "liquid water 0-30C", "viscosity");

        std::vector<std::string> controls({ "pressure" });
        out_list.sublist(pk_names[0]).sublist("time integrator")
            .set<Teuchos::Array<std::string> >("error control options", controls);

        out_list.sublist(pk).sublist("physical models and assumptions")
            .set<bool>("vapor diffusion", (pk_model_["energy"] == "two_phase energy"));

        err_options = "pressure, temperature";
      }
      else if (basename == "flow and shallow water") {
        Teuchos::Array<std::string> pk_names;
        pk_names.push_back(Keys::merge(prefix, "flow", delimiter));
        pk_names.push_back(Keys::merge(prefix, "shallow water", delimiter));
        out_list.sublist(pk)
            .set<Teuchos::Array<std::string> >("PKs order", pk_names)
            .set<int>("master PK index", 0);
      }
      else if (basename == "flow and energy fracture") {
        Teuchos::Array<std::string> pk_names;
        pk_names.push_back(Keys::merge(prefix, "flow fracture", delimiter));
        pk_names.push_back(Keys::merge(prefix, "energy fracture", delimiter));
        out_list.sublist(pk)
            .set<Teuchos::Array<std::string> >("PKs order", pk_names)
            .set<int>("master PK index", 0)
            .set<std::string>("domain name", "fracture");

        auto& tmp = glist.sublist("state").sublist("evaluators");
        AddSecondaryFieldEvaluator_(tmp, 
            Keys::getKey("fracture", "molar_density_liquid"), "molar density key",
            "eos", "liquid water 0-30C", "density");
        AddSecondaryFieldEvaluator_(tmp, 
            Keys::getKey("fracture", "viscosity_liquid"), "viscosity key",
            "viscosity", "liquid water 0-30C", "viscosity");

        std::vector<std::string> controls({ "pressure" });
        out_list.sublist(pk_names[0]).sublist("time integrator")
            .set<Teuchos::Array<std::string> >("error control options", controls);

        out_list.sublist(pk).sublist("physical models and assumptions")
            .set<bool>("vapor diffusion", (pk_model_["energy"] == "two_phase energy"));

        err_options = "pressure, temperature";
      }
      // -- multiple PKs (matrix and fracture)
      else if (basename == "coupled flow and transport") {
        Teuchos::Array<std::string> pk_names;
        pk_names.push_back(Keys::merge(prefix, "coupled flow", delimiter));
        pk_names.push_back(Keys::merge(prefix, "coupled transport", delimiter));
        out_list.sublist(pk).set<Teuchos::Array<std::string> >("PKs order", pk_names);
        out_list.sublist(pk).set<int>("master PK index", 0);
      }
      else if (basename == "coupled reactive transport") {
        Teuchos::Array<std::string> pk_names;
        pk_names.push_back(Keys::merge(prefix, "coupled chemistry", delimiter));
        pk_names.push_back(Keys::merge(prefix, "coupled transport", delimiter));
        out_list.sublist(pk).set<Teuchos::Array<std::string> >("PKs order", pk_names);
        out_list.sublist(pk).set<int>("master PK index", 0);
      }
      else if (basename == "coupled flow and energy") {
        Teuchos::Array<std::string> pk_names;
        pk_names.push_back(Keys::merge(prefix, "flow and energy", delimiter));
        pk_names.push_back(Keys::merge(prefix, "flow and energy fracture", delimiter));
        out_list.sublist(pk).set<Teuchos::Array<std::string> >("PKs order", pk_names);
        out_list.sublist(pk).set<int>("master PK index", 0);

        err_options = "pressure, temperature";
      }
      else if (basename == "coupled flow and reactive transport") {
        Teuchos::Array<std::string> pk_names;
        pk_names.push_back(Keys::merge(prefix, "coupled flow", delimiter));
        pk_names.push_back(Keys::merge(prefix, "coupled reactive transport", delimiter));
        out_list.sublist(pk).set<Teuchos::Array<std::string> >("PKs order", pk_names);
        out_list.sublist(pk).set<int>("master PK index", 0);
      }
      else if (basename == "coupled thermal flow and reactive transport") {
        Teuchos::Array<std::string> pk_names;
        pk_names.push_back(Keys::merge(prefix, "coupled flow and energy", delimiter));
        pk_names.push_back(Keys::merge(prefix, "coupled reactive transport", delimiter));
        out_list.sublist(pk)
            .set<Teuchos::Array<std::string> >("PKs order", pk_names)
            .set<int>("master PK index", 0);
      }
      else if (basename == "shallow water and transport") {
        Teuchos::Array<std::string> pk_names;
        pk_names.push_back(Keys::merge(prefix, "shallow water", delimiter));
        pk_names.push_back(Keys::merge(prefix, "transport", delimiter));
        out_list.sublist(pk)
            .set<Teuchos::Array<std::string> >("PKs order", pk_names)
            .set<std::string>("domain name", pk_domain_["shallow water"]);

        // ... coupling information 
        out_list.sublist(pk_names[1]).sublist("physical models and assumptions")
           .set<std::string>("darcy flux key", Keys::getKey(pk_domain_["shallow_water"], "riemann_flux"))
           .set<std::string>("saturation key", Keys::getKey(pk_domain_["shallow_water"], "ponded_depth"));
      }

      // add time integrator to PKs that have no transport and chemistry sub-PKs.
      if (pk.find("transport") == std::string::npos &&
          pk.find("chemistry") == std::string::npos) {
        if (!out_list.sublist(pk).isSublist("time integrator")) {
          out_list.sublist(pk).sublist("time integrator") = TranslateTimeIntegrator_(
              err_options, "nka", false,
              "unstructured_controls, unstr_transient_controls",
              dt_cut_[mode], dt_inc_[mode]);
          out_list.sublist(pk).sublist("verbose object") = verb_list_.sublist("verbose object");
        }
      }

      // lookup table
      if (eos_lookup_table_ != "") {
        out_list.sublist(pk).sublist("physical models and assumptions")
            .set<std::string>("eos lookup table", eos_lookup_table_);
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
    std::string prefix = Keys::split(name, delimiter).first;
    std::string basename = Keys::split(name, delimiter).second;

    if (basename == "coupled flow") {
      std::string flow_m = Keys::merge(prefix, "flow matrix", delimiter);
      std::string flow_f = Keys::merge(prefix, "flow fracture", delimiter);
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
                                  .set<std::string>("extraction method", "manifold mesh")
                                  .set<std::string>("domain name", "fracture");

      if (dim_ == 3) mesh_list.sublist("expert").set<bool>("request edges", true);
    }

    if (basename == "flow and shallow water") {
      mesh_list.sublist("expert").set<bool>("request edges", true);

      Teuchos::Array<std::string> aux(1, pk_region_["shallow_water"]);
      mesh_list.sublist("submesh").set<Teuchos::Array<std::string> >("regions", aux)
                                  .set<std::string>("extraction method", "manifold mesh")
                                  .set<std::string>("domain name", "surface");
    }

    if (basename == "coupled transport" && transport_implicit_) {
      std::string transport_m = Keys::merge(prefix, "transport matrix", delimiter);
      pk_list.sublist(transport_m).sublist("operators").sublist("advection operator")
             .sublist("matrix").set<Teuchos::Array<std::string> >("fracture", fracture_regions_);
    }

    if (basename == "coupled flow and energy") {
      std::vector<std::string> name_pks = { "flow and energy", "flow and energy fracture",
                                            "flow", "energy", "flow fracture", "energy fracture"};

      for (auto&& pk : name_pks) {
        std::string fullname = Keys::merge(prefix, pk, delimiter);
        auto& tmp = pk_list.sublist(fullname).sublist("time integrator");
        tmp.set<std::string>("time integration method", "none");
        tmp.remove("BDF1", false);
        tmp.remove("initialization", false);

        if (pk.find("fracture") == std::string::npos) {
          if (pk_list.sublist(fullname).isSublist("operators")) {
            auto& oplist = pk_list.sublist(fullname).sublist("operators").sublist("diffusion operator");
            oplist.sublist("matrix").set<Teuchos::Array<std::string> >("fracture", fracture_regions_);
            oplist.sublist("preconditioner").set<Teuchos::Array<std::string> >("fracture", fracture_regions_);
            oplist.sublist("vapor matrix").set<Teuchos::Array<std::string> >("fracture", fracture_regions_);
          }
        }

        if (pk == "flow" || pk == "flow fracture") {
          tmp.sublist("pressure-lambda constraints").set<std::string>("method", "projection");
        }

        if (pk == "energy") {
          auto& oplist = pk_list.sublist(fullname).sublist("operators").sublist("advection operator");
          oplist.set<Teuchos::Array<std::string> >("fracture", fracture_regions_);
        }
      }

      Teuchos::Array<std::string> aux(1, "FRACTURE_NETWORK_INTERNAL");
      mesh_list.sublist("submesh").set<Teuchos::Array<std::string> >("regions", aux)
                                  .set<std::string>("extraction method", "manifold mesh")
                                  .set<std::string>("domain name", "fracture");

      if (dim_ == 3) mesh_list.sublist("expert").set<bool>("request edges", true);
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
