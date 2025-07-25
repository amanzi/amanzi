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
InputConverterU::TranslateShallowWater_(const std::string& domain)
{
  Teuchos::ParameterList out_list;

  if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH)
    *vo_->os() << "Translating shallow water, domain=" << domain << std::endl;

  MemoryManager mm;
  DOMNode* node;

  char* text;

  // set up default values for some expert parameters
  double cfl(0.5), limiter_cfl(0.5);

  // process expert parameters
  bool flag;
  std::string flux("Rusanov"), weight("constant");
  std::string controls("unstructured_controls, unstr_shallow_water_controls");
  node = GetUniqueElementByTagsString_(controls + ", cfl", flag);
  if (flag) {
    text = mm.transcode(node->getTextContent());
    cfl = strtod(text, NULL);
  }
  node = GetUniqueElementByTagsString_(controls + ", limiter_cfl", flag);
  if (flag) {
    text = mm.transcode(node->getTextContent());
    limiter_cfl = strtod(text, NULL);
  }

  node = GetUniqueElementByTagsString_(controls + ", reconstruction_weight", flag);
  if (flag) {
    weight = GetTextContentS_(node, "constant, inverse-distance");
    std::replace(weight.begin(), weight.end(), '-', ' ');
  }

  node = GetUniqueElementByTagsString_(controls + ", numerical_flux", flag);
  if (flag) {
    flux = GetTextContentS_(node, "Rusanov, central upwind");
  }

  out_list.set<std::string>("domain name", (domain == "matrix") ? "domain" : domain)
    .set<std::string>("numerical flux", flux)
    .set<int>("number of reduced cfl cycles", 10)
    .set<double>("cfl", cfl);

  out_list.sublist("reconstruction")
    .set<int>("polynomial order", 1)
    .set<std::string>("weight", weight)
    .set<std::string>("limiter", "Barth-Jespersen")
    .set<std::string>("limiter stencil", "cell to closest cells")
    .set<std::string>("limiter location", "face")
    .set<double>("limiter cfl", limiter_cfl);

  int nspace(1), ntime(1);
  std::string tags_default("unstructured_controls, unstr_shallow_water_controls");
  node = GetUniqueElementByTagsString_(tags_default + ", algorithm", flag);
  if (flag) {
    std::string order = GetTextContentS_(node, "explicit first-order, explicit second-order");
    if (order == "explicit first-order") {
      nspace = 1;
      ntime = 1;
    } else if (order == "explicit second-order") {
      nspace = 2;
      ntime = 2;
    }
  }
  out_list.set<int>("spatial discretization order", nspace);
  out_list.set<int>("temporal discretization order", ntime);

  // boundary conditions
  out_list.sublist("boundary conditions") = TranslateShallowWaterBCs_();

  // pipe extension
  if (strcmp(pk_model_["shallow_water"].begin()->c_str(), "pipe flow") == 0) {
    out_list.set<std::string>("domain name", "domain");
    out_list.set<std::string>("diameter key", "diameter")
      .set<std::string>("direction key", "direction");
  }

  out_list.sublist("verbose object") = verb_list_.sublist("verbose object");

  return out_list;
}


/* ******************************************************************
* Create list of BCs.
****************************************************************** */
Teuchos::ParameterList
InputConverterU::TranslateShallowWaterBCs_()
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

    node = GetUniqueElementByTagsString_(inode, "liquid_phase, liquid_component", flag);
    if (!flag) continue;

    // process a group of similar elements defined by the first element
    std::string bctype, bcname;
    std::vector<DOMNode*> same_list = GetSameChildNodes_(node, bctype, flag, true);
    if (bctype != "ponded_depth") continue;

    // create vectors of values and forms
    auto bcs = ParseCondList_(node, 0.0, 1000.0, "m");

    // .. modify data
    if (bctype == "ponded_depth") {
      bctype = "ponded depth";
      bcname = "ponded depth";
    }
    std::stringstream ss;
    ss << "BC " << ibc++;

    // ponded depth
    {
      Teuchos::ParameterList& tbc_list = out_list.sublist(bctype);
      Teuchos::ParameterList& bc = tbc_list.sublist(ss.str());
      bc.set<Teuchos::Array<std::string>>("regions", regions)
        .set<std::string>("spatial distribution method", "none");

      Teuchos::ParameterList& bcfn = bc.sublist(bcname);
      TranslateGenericMath_(bcs, bcfn);
    }

    // velocity FIXME
    {
      bctype = "velocity";
      bcname = "velocity";
      Teuchos::ParameterList& tbc_list = out_list.sublist(bctype);
      Teuchos::ParameterList& bc = tbc_list.sublist(ss.str());
      bc.set<Teuchos::Array<std::string>>("regions", regions)
        .set<std::string>("spatial distribution method", "none");

      Teuchos::ParameterList& bcfn = bc.sublist(bcname);
      bcfn.set<int>("number of dofs", 2).set<std::string>("function type", "composite function");

      bcfn.sublist("dof 1 function").sublist("function-constant").set<double>("value", 0.0);
      bcfn.sublist("dof 2 function").sublist("function-constant").set<double>("value", 0.0);
    }
  }

  return out_list;
}

} // namespace AmanziInput
} // namespace Amanzi
