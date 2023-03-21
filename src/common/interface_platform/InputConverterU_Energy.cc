/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Input Converter

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
Teuchos::ParameterList
InputConverterU::TranslateEnergy_(const std::string& domain)
{
  Teuchos::ParameterList out_list;

  if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH)
    *vo_->os() << "Translating energy, domain=" << domain << std::endl;

  MemoryManager mm;
  DOMNode* node;
  DOMElement* element;

  // process expert parameters
  bool flag;
  node = GetUniqueElementByTagsString_("unstructured_controls, unstr_energy_controls", flag);

  // create flow header
  out_list.set<std::string>("domain name", (domain == "matrix") ? "domain" : domain);

  // insert operator sublist
  std::string disc_method("mfd-optimized_for_sparsity");
  node = GetUniqueElementByTagsString_(
    "unstructured_controls, unstr_energy_controls, discretization_method", flag);
  if (flag) disc_method = mm.transcode(node->getNodeName());

  std::string pc_method("linearized_operator");
  node = GetUniqueElementByTagsString_(
    "unstructured_controls, unstr_energy_controls, preconditioning_strategy", flag);
  if (flag) pc_method = mm.transcode(node->getNodeName());

  std::string nonlinear_solver("nka");
  node = GetUniqueElementByTagsString_("unstructured_controls, unstr_nonlinear_solver", flag);
  element = static_cast<DOMElement*>(node);
  if (flag) nonlinear_solver = GetAttributeValueS_(element, "name", TYPE_NONE, false, "nka");

  bool modify_correction(false);
  node = GetUniqueElementByTagsString_(
    "unstructured_controls, unstr_nonlinear_solver, modify_correction", flag);

  out_list.sublist("operators") =
    TranslateDiffusionOperator_(disc_method, pc_method, nonlinear_solver, "", "", false, "energy");

  // insert thermal conductivity evaluator with the default values (no 2.2 support yet)
  Teuchos::ParameterList& thermal =
    out_list.sublist("thermal conductivity evaluator").sublist("thermal conductivity parameters");
  if (pk_model_["energy"] == "two-phase energy") {
    thermal.set<std::string>("thermal conductivity type", "two-phase Peters-Lidard");
    thermal.set<double>("thermal conductivity of gas", 0.02);
    thermal.set<double>("unsaturated alpha", 1.0);
    thermal.set<double>("epsilon", 1.0e-10);
  } else {
    thermal.set<std::string>("thermal conductivity type", "one-phase polynomial");
  }

  double cv_f(0.606), cv_r(0.2);
  node = GetUniqueElementByTagsString_("materials", flag);
  std::vector<DOMNode*> materials = GetChildren_(node, "material", flag);

  node = GetUniqueElementByTagsString_(materials[0], "thermal_properties, rock_conductivity", flag);
  if (flag) cv_f = GetTextContentD_(node, "W/m/K", true);

  node = GetUniqueElementByTagsString_(materials[0], "thermal_properties, rock_conductivity", flag);
  if (flag) cv_r = GetTextContentD_(node, "W/m/K", true);

  thermal.set<double>("thermal conductivity of liquid", cv_f);
  thermal.set<double>("thermal conductivity of rock", cv_r);
  thermal.set<double>("reference temperature", 298.15);

  // insert time integrator
  std::string err_options("energy"), unstr_controls("unstructured_controls, unstr_energy_controls");

  if (pk_master_.find("energy") != pk_master_.end()) {
    out_list.sublist("time integrator") = TranslateTimeIntegrator_(err_options,
                                                                   nonlinear_solver,
                                                                   modify_correction,
                                                                   unstr_controls,
                                                                   TI_SOLVER,
                                                                   TI_TS_REDUCTION_FACTOR,
                                                                   TI_TS_INCREASE_FACTOR);
  }

  // insert boundary conditions and source terms
  out_list.sublist("boundary conditions") = TranslateEnergyBCs_(domain);
  out_list.sublist("source terms") = TranslateSources_(domain, "energy");

  // insert internal evaluators
  out_list.sublist("energy evaluator").sublist("verbose object") =
    verb_list_.sublist("verbose object");
  out_list.sublist("enthalpy evaluator").sublist("verbose object") =
    verb_list_.sublist("verbose object");

  // cross coupling of PKs
  if (fracture_regions_.size() > 0 && domain == "fracture") {
    out_list.sublist("physical models and assumptions")
      .set<bool>("flow and transport in fractures", true);
  }

  out_list.sublist("verbose object") = verb_list_.sublist("verbose object");
  return out_list;
}


/* ******************************************************************
* Create list of energy BCs.
****************************************************************** */
Teuchos::ParameterList
InputConverterU::TranslateEnergyBCs_(const std::string& domain)
{
  Teuchos::ParameterList out_list;

  MemoryManager mm;

  char* text;
  DOMNodeList* children;
  DOMNode* node;
  DOMElement* element;

  // correct list of boundary conditions for given domain
  bool flag;
  if (domain == "matrix")
    node = GetUniqueElementByTagsString_("boundary_conditions", flag);
  else
    node = GetUniqueElementByTagsString_("fracture_network, boundary_conditions", flag);
  if (!flag) return out_list;

  int ibc(0);
  children = node->getChildNodes();
  int nchildren = children->getLength();

  for (int i = 0; i < nchildren; ++i) {
    DOMNode* inode = children->item(i);
    if (inode->getNodeType() != DOMNode::ELEMENT_NODE) continue;

    // read the assigned regions
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
    std::map<double, std::string> tp_forms, tp_formulas;

    for (int j = 0; j < same_list.size(); ++j) {
      DOMNode* jnode = same_list[j];
      element = static_cast<DOMElement*>(jnode);
      double t0 = GetAttributeValueD_(element, "start", TYPE_TIME, DVAL_MIN, DVAL_MAX, "y");

      tp_forms[t0] = GetAttributeValueS_(element, "function");
      tp_values[t0] =
        GetAttributeValueD_(element, "value", TYPE_NUMERICAL, 0.0, 1000.0, "K", false, 0.0);
      tp_formulas[t0] = GetAttributeValueS_(element, "formula", TYPE_NONE, false, "");
    }

    // create vectors of values and forms
    std::vector<double> times, values;
    std::vector<std::string> forms, formulas;
    for (std::map<double, double>::iterator it = tp_values.begin(); it != tp_values.end(); ++it) {
      times.push_back(it->first);
      values.push_back(it->second);
      forms.push_back(tp_forms[it->first]);
      formulas.push_back(tp_formulas[it->first]);
    }

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
    bc.set<Teuchos::Array<std::string>>("regions", regions)
      .set<std::string>("spatial distribution method", "none");

    Teuchos::ParameterList& bcfn = bc.sublist(bcname);
    TranslateGenericMath_(times, values, forms, formulas, bcfn);
  }

  return out_list;
}

} // namespace AmanziInput
} // namespace Amanzi
