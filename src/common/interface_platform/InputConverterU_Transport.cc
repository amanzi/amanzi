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
* Create transport list.
****************************************************************** */
Teuchos::ParameterList InputConverterU::TranslateTransport_()
{
  Teuchos::ParameterList out_list;

  if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH)
    *vo_->os() << "Translating transport" << std::endl;

  XString mm;

  char *text, *tagname;
  DOMNodeList *node_list, *children;
  DOMNode* node;

  // process CFL number
  bool flag;
  double cfl(1.0);

  node = getUniqueElementByTagNames_("unstructured_controls", "unstr_transport_controls", "cfl", flag);
  if (flag) {
    text = mm.transcode(node->getTextContent());
    cfl = strtod(text, NULL);
  }

  // set defaults for transport
  out_list.set<int>("spatial discretization order", 1);
  out_list.set<int>("temporal discretization order", 1);
  out_list.set<double>("cfl", cfl);
  out_list.set<std::string>("flow mode", "transient");

  out_list.set<std::string>("solver", "PCG with Hypre AMG");
  out_list.set<std::string>("enable internal tests", "no");
  out_list.set<bool>("transport subcycling", TRANSPORT_SUBCYCLING);

  // overwrite data from expert parameters  
  node = getUniqueElementByTagNames_("unstructured_controls", "unstr_transport_controls", "sub_cycling", flag);
  text = mm.transcode(node->getTextContent());
  if (flag) {
    out_list.set<bool>("transport subcycling", (strcmp(text, "on") == 0));
  }

  int poly_order(0);
  node = getUniqueElementByTagNames_("unstructured_controls", "unstr_transport_controls", "algorithm", flag);
  text = mm.transcode(node->getTextContent());
  if (strcmp(text, "explicit first-order") == 0) {
    out_list.set<int>("spatial discretization order", 1);
    out_list.set<int>("temporal discretization order", 1);
  } else if (strcmp(text, "explicit second-order") == 0) {
    out_list.set<int>("spatial discretization order", 2);
    out_list.set<int>("temporal discretization order", 2);
    poly_order = 1;
  } else {
    ThrowErrorMissattr_("unstructured_controls", "element", "explicit first-order", "unstr_transport_controls");
  }

  Teuchos::ParameterList& trp_lift = out_list.sublist("reconstruction");
  trp_lift.set<int>("polynomial order", poly_order);
  trp_lift.set<std::string>("limiter", "tensorial");
  trp_lift.set<bool>("limiter extension for transport", true);

  // check if we need to write a dispersivity sublist
  bool dispersion = doc_->getElementsByTagName(mm.transcode("dispersion_tensor"))->getLength() > 0;

  if (!dispersion) {
    dispersion = doc_->getElementsByTagName(mm.transcode("tortuosity"))->getLength() > 0;
  }

  // create dispersion list
  if (dispersion) {
    node_list = doc_->getElementsByTagName(mm.transcode("materials"));

    Teuchos::ParameterList& mat_list = out_list.sublist("material properties");
    mat_list.set<std::string>("numerical method", "two-point flux approximation");

    children = node_list->item(0)->getChildNodes();
    for (int i = 0; i < children->getLength(); ++i) {
      DOMNode* inode = children->item(i);
      if (inode->getNodeType() != DOMNode::ELEMENT_NODE) continue;

      tagname = mm.transcode(inode->getNodeName());
      if (strcmp(tagname, "material") != 0) continue;

      // -- regions
      node = getUniqueElementByTagNames_(inode, "assigned_regions", flag);
      std::vector<std::string> regions = CharToStrings_(mm.transcode(node->getTextContent()));

      Teuchos::ParameterList& tmp_list = mat_list.sublist(tagname);
      tmp_list.set<Teuchos::Array<std::string> >("regions", regions);

      // -- dispersion tensor
      node = getUniqueElementByTagNames_(inode, "mechanical_properties", "dispersion_tensor", flag);
      if (flag) {
        double al, alh, alv, at, ath, atv;
        std::string model = GetAttributeValueS_(static_cast<DOMElement*>(node), "type");
        if (strcmp(model.c_str(), "uniform_isotropic") == 0) { 
          tmp_list.set<std::string>("model", "Bear");

          al = GetAttributeValueD_(static_cast<DOMElement*>(node), "alpha_l");
          at = GetAttributeValueD_(static_cast<DOMElement*>(node), "alpha_t");

          tmp_list.sublist("parameters for Bear").set<double>("alphaL", al)
                                                 .set<double>("alphaT", at);
        } else if (strcmp(model.c_str(), "burnett_frind") == 0) {
          tmp_list.set<std::string>("model", "Burnett-Frind");

          al = GetAttributeValueD_(static_cast<DOMElement*>(node), "alpha_l");
          ath = GetAttributeValueD_(static_cast<DOMElement*>(node), "alpha_th");
          ath = GetAttributeValueD_(static_cast<DOMElement*>(node), "alpha_tv");

          tmp_list.sublist("parameters for Burnett-Frind")
              .set<double>("alphaL", al).set<double>("alphaTH", ath)
              .set<double>("alphaTH", atv);
        } else if (strcmp(model.c_str(), "lichtner_kelkar_robinson") == 0) {
          tmp_list.set<std::string>("model", "Lichtner-Kelkar-Robinson");

          alh = GetAttributeValueD_(static_cast<DOMElement*>(node), "alpha_lh");
          alv = GetAttributeValueD_(static_cast<DOMElement*>(node), "alpha_lv");
          ath = GetAttributeValueD_(static_cast<DOMElement*>(node), "alpha_th");
          ath = GetAttributeValueD_(static_cast<DOMElement*>(node), "alpha_tv");

          tmp_list.sublist("parameters for Lichtner-Kelkar-Robinson")
              .set<double>("alphaLH", alh).set<double>("alphaLV", alv)
              .set<double>("alphaTH", ath).set<double>("alphaTH", atv);
        } 
      }

      // -- tortousity
      node = getUniqueElementByTagNames_(inode, "mechanical_properties", "tortuosity", flag);
      if (flag) {
        double val = GetAttributeValueD_(static_cast<DOMElement*>(node), "value");
        tmp_list.set<double>("aqueous tortuosity", val);
      }
    }
  }

  // -- molecular diffusion
  //    check in phases->water list (other solutes are ignored)
  node = getUniqueElementByTagsString_("phases, liquid_phase, dissolved_components, solutes", flag);
  if (flag) {
    Teuchos::ParameterList& diff_list = out_list.sublist("molecular diffusion");
    std::vector<std::string> aqueous_names;
    std::vector<double> aqueous_values;

    children = node->getChildNodes();
    for (int i = 0; i < children->getLength(); ++i) {
      DOMNode* inode = children->item(i);
      if (inode->getNodeType() != DOMNode::ELEMENT_NODE) continue;

      tagname = mm.transcode(inode->getNodeName());
      if (strcmp(tagname, "solute") != 0) continue;

      double val = GetAttributeValueD_(static_cast<DOMElement*>(inode), "coefficient_of_diffusion");
      text = mm.transcode(inode->getTextContent());

      aqueous_names.push_back(TrimString_(text));
      aqueous_values.push_back(val);
    }

    diff_list.set<Teuchos::Array<std::string> >("aqueous names", aqueous_names);
    diff_list.set<Teuchos::Array<double> >("aqueous values", aqueous_values);
  }

  // -- molecular diffusion
  //    check in phases->air list (other solutes are ignored)
  node = getUniqueElementByTagsString_("phases, gas_phase, dissolved_components, solutes", flag);
  if (flag) {
    Teuchos::ParameterList& diff_list = out_list.sublist("molecular diffusion");
    std::vector<std::string> gaseous_names;
    std::vector<double> gaseous_values;

    children = node->getChildNodes();
    for (int i = 0; i < children->getLength(); ++i) {
      DOMNode* inode = children->item(i);
      if (inode->getNodeType() != DOMNode::ELEMENT_NODE) continue;

      tagname = mm.transcode(inode->getNodeName());
      if (strcmp(tagname, "solute") != 0) continue;

      double val = GetAttributeValueD_(static_cast<DOMElement*>(inode), "coefficient_of_diffusion");
      text = mm.transcode(inode->getTextContent());

      gaseous_names.push_back(TrimString_(text));
      gaseous_values.push_back(val);
    }

    diff_list.set<Teuchos::Array<std::string> >("gaseous names", gaseous_names);
    diff_list.set<Teuchos::Array<double> >("gaseous values", gaseous_values);
  }

  // create the sources and boundary conditions lists
  out_list.sublist("boundary conditions") = TranslateTransportBCs_();
  out_list.sublist("source terms") = TranslateTransportSources_();

  // remaining global parameters
  out_list.set<int>("number of aqueous components", phases_["water"].size());
  out_list.set<int>("number of gaseous components", phases_["gas"].size());

  out_list.sublist("VerboseObject") = verb_list_.sublist("VerboseObject");
  return out_list;
}


/* ******************************************************************
* Create list of transport sources.
****************************************************************** */
Teuchos::ParameterList InputConverterU::TranslateTransportBCs_()
{
  Teuchos::ParameterList out_list;

  XString mm;

  char *text, *tagname;
  DOMNodeList *node_list, *children;
  DOMNode *node, *phase;
  DOMElement* element;

  node_list = doc_->getElementsByTagName(mm.transcode("boundary_conditions"));
  if (!node_list) return out_list;

  children = node_list->item(0)->getChildNodes();

  for (int i = 0; i < children->getLength(); ++i) {
    DOMNode* inode = children->item(i);
    if (inode->getNodeType() != DOMNode::ELEMENT_NODE) continue;
    tagname = mm.transcode(inode->getNodeName());
    std::string bcname = GetAttributeValueS_(static_cast<DOMElement*>(inode), "name");

    // read the assigned regions
    bool flag;
    node = getUniqueElementByTagsString_(inode, "assigned_regions", flag);
    text = mm.transcode(node->getTextContent());
    std::vector<std::string> regions = CharToStrings_(text);

    vv_bc_regions_.insert(vv_bc_regions_.end(), regions.begin(), regions.end());

    phase = getUniqueElementByTagNames_(inode, "liquid_phase", flag);
    if (!flag) continue;

    // process solute elements
    // -- Dirichlet BCs for concentration
    std::string bctype, solute_name;

    element = static_cast<DOMElement*>(phase);
    DOMNodeList* solutes = element->getElementsByTagName(mm.transcode("solute_component"));

    if (solutes->getLength() > 0) {
      for (int n = 0; n < solutes->getLength(); ++n) {
        node = solutes->item(n);
        // process a group of similar elements defined by the first element
        std::vector<DOMNode*> same_list = getSameChildNodes_(node, bctype, flag, true);
        solute_name = GetAttributeValueS_(static_cast<DOMElement*>(same_list[0]), "name");

        std::map<double, double> tp_values;
        std::map<double, std::string> tp_forms;

        for (int j = 0; j < same_list.size(); ++j) {
          element = static_cast<DOMElement*>(same_list[j]);
          double t0 = GetAttributeValueD_(element, "start");
          tp_forms[t0] = GetAttributeValueS_(element, "function");
          tp_values[t0] = GetAttributeValueD_(element, "value");
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
     
        // save in the XML files  
        Teuchos::ParameterList& tbc_list = out_list.sublist("concentration");
        Teuchos::ParameterList& bc = tbc_list.sublist(solute_name).sublist(bcname);
        bc.set<Teuchos::Array<std::string> >("regions", regions);

        Teuchos::ParameterList& bcfn = bc.sublist("boundary concentration").sublist("function-tabular");
        bcfn.set<Teuchos::Array<double> >("x values", times);
        bcfn.set<Teuchos::Array<double> >("y values", values);
        bcfn.set<Teuchos::Array<std::string> >("forms", forms);
      }
    }

    // -- geochemistry BCs 
    node = getUniqueElementByTagNames_(phase, "geochemistry", flag);
    if (flag) {
      std::vector<DOMNode*> same_list = getSameChildNodes_(node, bctype, flag, true);
      std::map<double, std::string> tp_forms, tp_values;

      for (int j = 0; j < same_list.size(); ++j) {
        element = static_cast<DOMElement*>(same_list[j]);
        double t0 = GetAttributeValueD_(element, "start");
        tp_forms[t0] = GetAttributeValueS_(element, "function");
        tp_values[t0] = GetAttributeValueS_(element, "value");
      }

      // create vectors of values and forms
      std::vector<double> times;
      std::vector<std::string> forms, values;
      for (std::map<double, std::string>::iterator it = tp_values.begin(); it != tp_values.end(); ++it) {
        times.push_back(it->first);
        values.push_back(it->second);
        forms.push_back(tp_forms[it->first]);
      }
     
      // save in the XML files  
      Teuchos::ParameterList& tbc_list = out_list.sublist("geochemical conditions");
      Teuchos::ParameterList& bc = tbc_list.sublist(bcname);
      bc.set<Teuchos::Array<std::string> >("regions", regions);
      bc.set<Teuchos::Array<double> >("Times", times);
      bc.set<Teuchos::Array<std::string> >("Geochemical Conditions", values);
      bc.set<Teuchos::Array<std::string> >("Time Functions", forms);
    }
  }

  return out_list;
}


/* ******************************************************************
* Create list of transport sources.
****************************************************************** */
Teuchos::ParameterList InputConverterU::TranslateTransportSources_()
{
  Teuchos::ParameterList out_list;

  XString mm;

  char *text, *tagname;
  DOMNodeList *node_list, *children;
  DOMNode *node, *phase;
  DOMElement* element;

  node_list = doc_->getElementsByTagName(mm.transcode("sources"));
  if (!node_list) return out_list;

  children = node_list->item(0)->getChildNodes();

  for (int i = 0; i < children->getLength(); ++i) {
    DOMNode* inode = children->item(i);
    if (inode->getNodeType() != DOMNode::ELEMENT_NODE) continue;
    tagname = mm.transcode(inode->getNodeName());
    std::string bcname = GetAttributeValueS_(static_cast<DOMElement*>(inode), "name");

    // read the assigned regions
    bool flag;
    node = getUniqueElementByTagsString_(inode, "assigned_regions", flag);
    text = mm.transcode(node->getTextContent());
    std::vector<std::string> regions = CharToStrings_(text);

    vv_src_regions_.insert(vv_src_regions_.end(), regions.begin(), regions.end());

    phase = getUniqueElementByTagNames_(inode, "liquid_phase", flag);
    if (!flag) continue;

    // process solute elements
    // -- Dirichlet BCs for concentration
    std::string bctype, solute_name;

    element = static_cast<DOMElement*>(phase);
    DOMNodeList* solutes = element->getElementsByTagName(mm.transcode("solute_component"));

    if (solutes->getLength() > 0) {
      for (int n = 0; n < solutes->getLength(); ++n) {
        node = solutes->item(n);
        // get a group of similar elements defined by the first element
        std::vector<DOMNode*> same_list = getSameChildNodes_(node, bctype, flag, true);
        solute_name = GetAttributeValueS_(static_cast<DOMElement*>(same_list[0]), "name");

        // weighting method
        bool classical(true);
        std::string weight;
        text = mm.transcode(same_list[0]->getNodeName());
        if (strcmp(text, "volume_weighted") == 0) {
          weight = "volume";
        } else if (strcmp(text, "perm_weighted") == 0) {
          weight = "permeability";
        } else if (strcmp(text, "diffusion_dominated_release") == 0) {
          classical = false;
        } else {
          ThrowErrorIllformed_("sources", "element", text);
        } 

        if (classical) {
          std::map<double, double> tp_values;
          std::map<double, std::string> tp_forms;

          for (int j = 0; j < same_list.size(); ++j) {
            element = static_cast<DOMElement*>(same_list[j]);
            double t0 = GetAttributeValueD_(element, "start");
            tp_forms[t0] = GetAttributeValueS_(element, "function");
            tp_values[t0] = GetAttributeValueD_(element, "value");
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
     
          // save in the XML files  
          Teuchos::ParameterList& src_list = out_list.sublist("concentration");
          Teuchos::ParameterList& src = src_list.sublist(solute_name).sublist(bcname);
          src.set<Teuchos::Array<std::string> >("regions", regions);
          src.set<std::string>("spatial distribution method", weight);

          src.sublist("boundary concentration").sublist("function-tabular")
              .set<Teuchos::Array<double> >("x values", times)
              .set<Teuchos::Array<double> >("y values", values)
              .set<Teuchos::Array<std::string> >("forms", forms);
        } else {
          element = static_cast<DOMElement*>(same_list[0]);
          double total = GetAttributeValueD_(element, "inventory");
          double diff = GetAttributeValueD_(element, "diffusion_coeff");
          double length = GetAttributeValueD_(element, "mixing_length");
          std::vector<double> times;
          times.push_back(GetAttributeValueD_(element, "start"));
          times.push_back(GetAttributeValueD_(element, "end"));

          // save data in the XML
          Teuchos::ParameterList& src_list = out_list.sublist("concentration");
          Teuchos::ParameterList& src = src_list.sublist(solute_name).sublist(bcname);
          src.set<Teuchos::Array<std::string> >("regions", regions);

          std::vector<double> values(2, 0.0);
          std::vector<std::string> forms(1, "SQRT");
          double amplitude = 2 * total / length * std::pow(diff / M_PI, 0.5); 

          Teuchos::ParameterList& srcfn = src.sublist("boundary concentration").sublist("function-tabular");
          srcfn.set<Teuchos::Array<double> >("x values", times)
               .set<Teuchos::Array<double> >("y values", values)
               .set<Teuchos::Array<std::string> >("forms", forms);

          Teuchos::ParameterList& func = srcfn.sublist("SQRT").sublist("function-standard-math");
          func.set<std::string>("operator", "sqrt")
              .set<double>("parameter", 0.5)
              .set<double>("amplitude", amplitude)
              .set<double>("shift", times[0]);
        }
      }
    }
  }

  return out_list;
}

}  // namespace AmanziInput
}  // namespace Amanzi


