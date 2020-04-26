/*
  Input Converter

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
* Create energy list.
****************************************************************** */
Teuchos::ParameterList InputConverterU::TranslateEnergy_()
{
  Teuchos::ParameterList out_list;

  if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH)
    *vo_->os() << "Translating energy" << std::endl;

  MemoryManager mm;
  DOMNode* node;
  DOMElement* element;

  // process expert parameters
  bool flag;
  node = GetUniqueElementByTagsString_("unstructured_controls, unstr_energy_controls", flag);

  // insert operator sublist
  std::string disc_method("mfd-optimized_for_sparsity");
  node = GetUniqueElementByTagsString_("unstructured_controls, unstr_energy_controls, discretization_method", flag);
  if (flag) disc_method = mm.transcode(node->getNodeName());

  std::string pc_method("linearized_operator");
  node = GetUniqueElementByTagsString_("unstructured_controls, unstr_energy_controls, preconditioning_strategy", flag);
  if (flag) pc_method = mm.transcode(node->getNodeName()); 

  std::string nonlinear_solver("nka");
  node = GetUniqueElementByTagsString_("unstructured_controls, unstr_nonlinear_solver", flag);
  element = static_cast<DOMElement*>(node);
  if (flag) nonlinear_solver = GetAttributeValueS_(element, "name", TYPE_NONE, false, "nka"); 

  bool modify_correction(false);
  node = GetUniqueElementByTagsString_("unstructured_controls, unstr_nonlinear_solver, modify_correction", flag);

  out_list.sublist("operators") = TranslateDiffusionOperator_(
      disc_method, pc_method, nonlinear_solver, "", "", false);

  // insert thermal conductivity evaluator with the default values (no 2.2 support yet)
  Teuchos::ParameterList& thermal = out_list.sublist("thermal conductivity evaluator")
                                            .sublist("thermal conductivity parameters");
  thermal.set<std::string>("thermal conductivity type", "two-phase Peters-Lidard");
  thermal.set<double>("thermal conductivity of rock", 0.2);
  thermal.set<double>("thermal conductivity of liquid", 0.1);
  thermal.set<double>("thermal conductivity of gas", 0.02);
  thermal.set<double>("unsaturated alpha", 1.0);
  thermal.set<double>("epsilon", 1.0e-10);

  // insert time integrator
  std::string err_options("energy"), unstr_controls("unstructured_controls, unstr_energy_controls");
  
  if (pk_master_.find("energy") != pk_master_.end()) {
    out_list.sublist("time integrator") = TranslateTimeIntegrator_(
        err_options, nonlinear_solver, modify_correction, unstr_controls,
        TI_TS_REDUCTION_FACTOR, TI_TS_INCREASE_FACTOR);  
  }

  // insert boundary conditions and source terms
  out_list.sublist("boundary conditions") = TranslateEnergyBCs_();

  // insert internal evaluators
  out_list.sublist("energy evaluator")
          .sublist("verbose object") = verb_list_.sublist("verbose object");
  out_list.sublist("enthalpy evaluator")
          .sublist("verbose object") = verb_list_.sublist("verbose object");

  out_list.sublist("verbose object") = verb_list_.sublist("verbose object");
  return out_list;
}


/* ******************************************************************
* Create list of energy BCs.
****************************************************************** */
Teuchos::ParameterList InputConverterU::TranslateEnergyBCs_()
{
  Teuchos::ParameterList out_list;

  MemoryManager mm;

  char *text;
  DOMNodeList *node_list, *children;
  DOMNode *node;
  DOMElement *element;

  node_list = doc_->getElementsByTagName(mm.transcode("boundary_conditions"));
  if (!node_list) return out_list;

  int ibc(0);
  children = node_list->item(0)->getChildNodes();
  int nchildren = children->getLength();

  for (int i = 0; i < nchildren; ++i) {
    DOMNode* inode = children->item(i);
    if (inode->getNodeType() != DOMNode::ELEMENT_NODE) continue;

    // read the assigned regions
    bool flag;
    node = GetUniqueElementByTagsString_(inode, "assigned_regions", flag);
    text = mm.transcode(node->getTextContent());
    std::vector<std::string> regions = CharToStrings_(text);

    vv_bc_regions_.insert(vv_bc_regions_.end(), regions.begin(), regions.end());

    node = GetUniqueElementByTagsString_(inode, "thermal_component", flag);
    if (!flag) continue;

    // process a group of similar elements defined by the first element
    std::string bctype;
    std::vector<DOMNode*> same_list = GetSameChildNodes_(node, bctype, flag, true);

    std::map<double, double> tp_values;
    std::map<double, std::string> tp_forms;

    for (int j = 0; j < same_list.size(); ++j) {
      DOMNode* jnode = same_list[j];
      element = static_cast<DOMElement*>(jnode);
      double t0 = GetAttributeValueD_(element, "start", TYPE_TIME, DVAL_MIN, DVAL_MAX, "y");

      tp_forms[t0] = GetAttributeValueS_(element, "function");
      tp_values[t0] = GetAttributeValueD_(element, "value", TYPE_NUMERICAL, 0.0, 1000.0, "K", false, 0.0);
    }

    // create vectors of values and forms
    std::vector<double> times, values;
    std::vector<std::string> forms;
    for (std::map<double, double>::iterator it = tp_values.begin(); it != tp_values.end(); ++it) {
      times.push_back(it->first);
      values.push_back(it->second);
      forms.push_back(tp_forms[it->first]);
    }
    forms.pop_back();

    // create names, modify data
    std::string bcname;
    if (bctype == "uniform_temperature") {
      bctype = "temperature";
      bcname = "boundary temperature";
    }
    std::stringstream ss;
    ss << "BC " << ibc++;

    // save in the XML files  
    Teuchos::ParameterList& tbc_list = out_list.sublist(bctype);
    Teuchos::ParameterList& bc = tbc_list.sublist(ss.str());
    bc.set<Teuchos::Array<std::string> >("regions", regions)
      .set<std::string>("spatial distribution method", "none");

    Teuchos::ParameterList& bcfn = bc.sublist(bcname);
    if (times.size() == 1) {
      bcfn.sublist("function-constant").set<double>("value", values[0]);
    } else {
      bcfn.sublist("function-tabular")
          .set<Teuchos::Array<double> >("x values", times)
          .set<Teuchos::Array<double> >("y values", values)
          .set<Teuchos::Array<std::string> >("forms", forms);
    }
  }

  return out_list;
}

}  // namespace AmanziInput
}  // namespace Amanzi


