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
InputConverterU::TranslateMechanics_(const std::string& domain)
{
  Teuchos::ParameterList out_list;

  if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH)
    *vo_->os() << "Translating mechanics, domain=" << domain << std::endl;

  MemoryManager mm;
  DOMNode *node;
  DOMElement* element;

  // process expert parameters
  bool flag;
  node = GetUniqueElementByTagsString_("unstructured_controls, unstr_mechanics_controls", flag);

  bool biot_model(false);
  node = GetUniqueElementByTagsString_(node, "biot_model", flag);
  if (flag) biot_model = GetTextContentL_(node, false);

  // create flow header
  out_list.set<std::string>("domain name", (domain == "matrix") ? "domain" : domain);

  // insert operator sublist
  std::string disc_method("mfd-default");
  std::string pc_method("linearized_operator");

  std::string nonlinear_solver("nka");
  node = GetUniqueElementByTagsString_("unstructured_controls, unstr_nonlinear_solver", flag);
  element = static_cast<DOMElement*>(node);
  if (flag) nonlinear_solver = GetAttributeValueS_(element, "name", TYPE_NONE, false, "nka");

  bool modify_correction(false);
  node = GetUniqueElementByTagsString_(
    "unstructured_controls, unstr_nonlinear_solver, modify_correction", flag);

  // insert time integrator
  std::string err_options("displacement"),
    unstr_controls("unstructured_controls, unstr_mechanics_controls");

  if (pk_master_.find("mechanics") != pk_master_.end()) {
    out_list.sublist("time integrator") = TranslateTimeIntegrator_(err_options,
                                                                   nonlinear_solver,
                                                                   modify_correction,
                                                                   unstr_controls,
                                                                   TI_SOLVER,
                                                                   TI_TS_REDUCTION_FACTOR,
                                                                   TI_TS_INCREASE_FACTOR);
  }

  // other parameters
  out_list.sublist("operators").sublist("elasticity operator")
    .set<std::string>("matrix type", "stiffness")
    .sublist("schema")
    .set<std::string>("base", "cell")
    .set<std::string>("method", "elasticity")
    .set<int>("method order", 1);

  out_list.sublist("physical models and assumptions")
    .set<bool>("use gravity", gravity_on_)
    .set<bool>("use biot model", biot_model);

  // insert boundary conditions and source terms
  out_list.sublist("boundary conditions") = TranslateMechanicsBCs_(domain);

  out_list.sublist("verbose object") = verb_list_.sublist("verbose object");
  return out_list;
}


/* ******************************************************************
* Create list of energy BCs.
****************************************************************** */
Teuchos::ParameterList
InputConverterU::TranslateMechanicsBCs_(const std::string& domain)
{
  Teuchos::ParameterList out_list;

  MemoryManager mm;

  char* text;
  DOMNodeList* children;
  DOMNode* node;

  // correct list of boundary conditions for given domain
  bool flag;
  node = GetUniqueElementByTagsString_("boundary_conditions", flag);
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

    node = GetUniqueElementByTagsString_(inode, "mechanics_component", flag);
    if (!flag) continue;

    // process a group of similar elements defined by the first element
    // create vectors of values and forms
    auto bcs = ParseCondList_(node, DVAL_MIN, DVAL_MAX, "m");

    std::stringstream ss;
    ss << "BC " << ibc++;

    // save in the XML files
    Teuchos::ParameterList& tbc_list = out_list.sublist(bcs.type);
    Teuchos::ParameterList& bc = tbc_list.sublist(ss.str());
    bc.set<Teuchos::Array<std::string>>("regions", regions)
      .set<std::string>("spatial distribution method", "none");

    if (bcs.type == "displacement") {
      Teuchos::ParameterList& bcfn = bc.sublist("no slip");
      bcfn.set<int>("number of dofs", dim_)
          .set<std::string>("function type", "composite function");
      for (int k = 0; k < dim_; ++k) {
        std::stringstream dof_str;
        dof_str << "dof " << k + 1 << " function";
        bcfn.sublist(dof_str.str()).sublist("function-constant").set<double>("value", bcs.vectors[0][k]);
      }
    } else if (bcs.type == "traction") {
      Teuchos::ParameterList& bcfn = bc.sublist("traction");
      bcfn.set<int>("number of dofs", dim_)
          .set<std::string>("function type", "composite function");
      for (int k = 0; k < dim_; ++k) {
        std::stringstream dof_str;
        dof_str << "dof " << k + 1 << " function";
        bcfn.sublist(dof_str.str()).sublist("function-constant").set<double>("value", bcs.vectors[0][k]);
      }
    } else if (bcs.type == "kinematic") {
      Teuchos::ParameterList& bcfn = bc.sublist("kinematic");
      bcfn.sublist("function-constant").set<double>("value", bcs.values[0]);
    }
  }

  return out_list;
}

} // namespace AmanziInput
} // namespace Amanzi
