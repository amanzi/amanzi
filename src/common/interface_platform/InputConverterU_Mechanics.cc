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
  DOMNode *node, *inode;
  DOMElement* element;

  // process expert parameters
  bool biot_undrained_split(false), biot_stress_split(false), use_fracture(false);
  std::string disc_method("elasticity");

  bool flag;
  node = GetUniqueElementByTagsString_("unstructured_controls, unstr_mechanics_controls", flag);
  if (flag) {
   // -- insert operator sublist
   inode = GetUniqueElementByTagsString_(node, "biot_model", flag);
    if (flag) {
      std::string method = GetTextContentS_(inode, "undrained_split, fixed_stress_split");
      biot_undrained_split = (method == "undrained_split");
      biot_stress_split = (method == "fixed_stress_split");
    }

    // -- discretization method
    inode = GetUniqueElementByTagsString_(node, "discretization_method", flag);
    if (flag) disc_method = GetTextContentS_(inode, "elasticity, BernardiRaugel");

    // -- insert fracture
    inode = GetUniqueElementByTagsString_(node, "use_fracture", flag);
    if (flag) use_fracture = GetTextContentL_(inode);
  }

  // create header
  out_list.set<std::string>("domain name", (domain == "matrix") ? "domain" : domain);

  std::string nonlinear_solver("nka");
  node = GetUniqueElementByTagsString_("unstructured_controls, unstr_nonlinear_solver", flag);
  element = static_cast<DOMElement*>(node);
  if (flag) nonlinear_solver = GetAttributeValueS_(element, "name", TYPE_NONE, false, "nka");

  bool modify_correction(false);
  node = GetUniqueElementByTagsString_(
    "unstructured_controls, unstr_nonlinear_solver, modify_correction", flag);

  // insert time integrator
  std::string err_options("displacement"),
    unstr_controls("unstructured_controls, unstr_transient_controls");

  if (pk_master_.find("mechanics") != pk_master_.end()) {
    out_list.sublist("time integrator") = TranslateTimeIntegrator_(err_options,
                                                                   nonlinear_solver,
                                                                   modify_correction,
                                                                   unstr_controls,
                                                                   "PCG for elasticity",
                                                                   TI_TS_REDUCTION_FACTOR,
                                                                   TI_TS_INCREASE_FACTOR);
  }

  // other parameters
  out_list.sublist("operators")
    .sublist("elasticity operator")
    .set<std::string>("matrix type", "stiffness")
    .sublist("schema")
    .set<std::string>("base", "cell")
    .set<std::string>("method", disc_method)
    .set<int>("method order", 1);

  if (fracture_network_ && use_fracture)
    out_list.sublist("operators")
      .sublist("elasticity operator")
      .sublist("schema")
      .set<Teuchos::Array<std::string>>("fracture", { "fracture" });

  out_list.sublist("physical models and assumptions")
    .set<bool>("use gravity", gravity_on_)
    .set<bool>("biot scheme: undrained split", biot_undrained_split)
    .set<bool>("biot scheme: fixed stress split", biot_stress_split);

  // small strain model
  auto tmp = TranslateMechanicsSSM_();
  if (tmp.numParams() > 0) out_list.sublist("small strain models") = tmp;

  // insert boundary conditions and source terms
  out_list.sublist("boundary conditions") = TranslateMechanicsBCs_(domain);

  out_list.sublist("verbose object") = verb_list_.sublist("verbose object");
  return out_list;
}


/* ******************************************************************
* Create list of permeability porosity models.
****************************************************************** */
Teuchos::ParameterList
InputConverterU::TranslateMechanicsSSM_()
{
  Teuchos::ParameterList out_list;

  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH)
    *vo_->os() << "Translating small strain models" << std::endl;

  MemoryManager mm;
  DOMNodeList* children;
  DOMNode* node;
  DOMElement* element;

  bool flag, found(false);

  node = GetUniqueElementByTagsString_("materials", flag);
  element = static_cast<DOMElement*>(node);
  children = element->getElementsByTagName(mm.transcode("material"));

  for (int i = 0; i < children->getLength(); ++i) {
    DOMNode* inode = children->item(i);

    bool flag;
    node = GetUniqueElementByTagsString_(inode, "assigned_regions", flag);
    std::vector<std::string> regions = CharToStrings_(mm.transcode(node->getTextContent()));

    // get optional compressibility
    node = GetUniqueElementByTagsString_(inode, "mechanical_properties, small_strain", flag);
    std::string model = GetAttributeValueS_(node, "model", TYPE_NONE, false, "");

    std::stringstream ss;
    ss << "SSM " << i;

    Teuchos::ParameterList& ssm_list = out_list.sublist(ss.str());
    ssm_list.set<Teuchos::Array<std::string>>("regions", regions);

    if (model == "hardin_drnevich") {
      found = true;
      double gamma = GetAttributeValueD_(node, "reference_shear_strain", TYPE_NUMERICAL, 0.0, DVAL_MAX, "Pa");
      double Gmax = GetAttributeValueD_(node, "maximum_shear_stress", TYPE_NUMERICAL, 0.0, DVAL_MAX, "Pa");

      ssm_list.set<std::string>("model", "Hardin Drnevich")
        .set<double>("reference shear strain", gamma)
        .set<double>("maximum shear stress", Gmax);
    }
  }

  if (!found) {
    Teuchos::ParameterList empty;
    out_list = empty;
  }

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
      bcfn.set<int>("number of dofs", dim_).set<std::string>("function type", "composite function");

      auto formulas = CharToStrings_(bcs.formulas[0].c_str());
      for (int k = 0; k < dim_; ++k) {
        std::stringstream dof_str;
        dof_str << "dof " << k + 1 << " function";
        if (formulas.size() == dim_) {
          bcfn.sublist(dof_str.str())
            .sublist("function-exprtk")
            .set<int>("number of arguments", dim_ + 1)
            .set<std::string>("formula", formulas[k]);
        } else {
          bcfn.sublist(dof_str.str())
            .sublist("function-constant")
            .set<double>("value", bcs.vectors[0][k]);
        }
      }
    } else if (bcs.type == "traction") {
      Teuchos::ParameterList& bcfn = bc.sublist("traction");
      bcfn.set<int>("number of dofs", dim_).set<std::string>("function type", "composite function");
      for (int k = 0; k < dim_; ++k) {
        std::stringstream dof_str;
        dof_str << "dof " << k + 1 << " function";
        bcfn.sublist(dof_str.str())
          .sublist("function-constant")
          .set<double>("value", bcs.vectors[0][k]);
      }
    } else if (bcs.type == "kinematic") {
      Teuchos::ParameterList& bcfn = bc.sublist("kinematic");
      bcfn.sublist("function-constant").set<double>("value", bcs.values[0]);
      if (bcs.kinematic != "") bc.set<std::string>("plane strain direction", bcs.kinematic);
    }
  }

  return out_list;
}

} // namespace AmanziInput
} // namespace Amanzi
