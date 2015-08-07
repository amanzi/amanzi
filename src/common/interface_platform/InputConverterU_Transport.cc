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
#include <xercesc/util/XMLString.hpp>

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
* Create flow list.
****************************************************************** */
Teuchos::ParameterList InputConverterU::TranslateTransport_()
{
  Teuchos::ParameterList out_list;

  if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH) {
    *vo_->os() << "Translating transport" << std::endl;
  }

  char *text_content, *tagname;
  XMLCh* xstr;
  DOMNodeList *node_list, *children;
  DOMNode* node;

  // process CFL number
  bool flag;
  double cfl(1.0);

  node = getUniqueElementByTagNames_("unstructured_controls", "unstr_transport_controls", "cfl", flag);
  if (flag) {
    text_content = XMLString::transcode(node->getTextContent());
    cfl = strtod(text_content, NULL);
    XMLString::release(&text_content);
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
  text_content = XMLString::transcode(node->getTextContent());
  if (flag) {
    out_list.set<bool>("transport subcycling", (strcmp(text_content, "on") == 0));
  }

  int poly_order(0);
  node = getUniqueElementByTagNames_("unstructured_controls", "unstr_transport_controls", "algorithm", flag);
  text_content = XMLString::transcode(node->getTextContent());
  if (strcmp(text_content, "explicit first-order") == 0) {
    out_list.set<int>("spatial discretization order", 1);
    out_list.set<int>("temporal discretization order", 1);
  } else if (strcmp(text_content, "explicit second-order") == 0) {
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
  xstr = XMLString::transcode("dispersion_tensor");
  bool dispersion = doc_->getElementsByTagName(xstr)->getLength() > 0;
  XMLString::release(&xstr);

  if (!dispersion) {
    xstr = XMLString::transcode("tortuosity");
    dispersion = doc_->getElementsByTagName(xstr)->getLength() > 0;
    XMLString::release(&xstr);
  }

  // create dispersion list
  if (dispersion) {
    xstr = XMLString::transcode("materials");
    node_list = doc_->getElementsByTagName(xstr);
    XMLString::release(&xstr);

    Teuchos::ParameterList& mat_list = out_list.sublist("material properties");
    mat_list.set<std::string>("numerical method", "two-point flux approximation");

    children = node_list->item(0)->getChildNodes();
    for (int i = 0; i < children->getLength(); ++i) {
      DOMNode* inode = children->item(i);
      if (inode->getNodeType() != DOMNode::ELEMENT_NODE) continue;

      tagname = XMLString::transcode(inode->getNodeName());
      if (strcmp(tagname, "material") != 0) continue;

      // -- regions
      node = getUniqueElementByTagNames_(inode, "assigned_regions", flag);
      std::vector<std::string> regions = GetElementVectorS_(static_cast<DOMElement*>(node));

      Teuchos::ParameterList& tmp_list = mat_list.sublist(tagname);
      tmp_list.set<Teuchos::Array<std::string> >("regions", regions);

      // -- dispersion tensor
      node = getUniqueElementByTagNames_(inode, "mechanical_properties", "dispersion_tensor", flag);
      if (flag) {
        double al, alh, alv, at, ath, atv;
        char* model = GetAttributeValueC_(static_cast<DOMElement*>(node), "type");
        if (strcmp(model, "uniform_isotropic") == 0) { 
          tmp_list.set<std::string>("model", "Bear");

          al = GetAttributeValueD_(static_cast<DOMElement*>(node), "alpha_l");
          at = GetAttributeValueD_(static_cast<DOMElement*>(node), "alpha_t");

          tmp_list.sublist("parameters for Bear").set<double>("alphaL", al)
                                                 .set<double>("alphaT", at);
        } else if (strcmp(model, "burnett_frind") == 0) {
          tmp_list.set<std::string>("model", "Burnett-Frind");

          al = GetAttributeValueD_(static_cast<DOMElement*>(node), "alpha_l");
          ath = GetAttributeValueD_(static_cast<DOMElement*>(node), "alpha_th");
          ath = GetAttributeValueD_(static_cast<DOMElement*>(node), "alpha_tv");

          tmp_list.sublist("parameters for Burnett-Frind")
              .set<double>("alphaL", al).set<double>("alphaTH", ath)
              .set<double>("alphaTH", atv);
        } else if (strcmp(model, "lichtner_kelkar_robinson") == 0) {
          tmp_list.set<std::string>("model", "Lichtner-Kelkar-Robinson");

          alh = GetAttributeValueD_(static_cast<DOMElement*>(node), "alpha_lh");
          alv = GetAttributeValueD_(static_cast<DOMElement*>(node), "alpha_lv");
          ath = GetAttributeValueD_(static_cast<DOMElement*>(node), "alpha_th");
          ath = GetAttributeValueD_(static_cast<DOMElement*>(node), "alpha_tv");

          tmp_list.sublist("parameters for Lichtner-Kelkar-Robinson")
              .set<double>("alphaLH", alh).set<double>("alphaLV", alv)
              .set<double>("alphaTH", ath).set<double>("alphaTH", atv);
        } 
        XMLString::release(&model);
      }

      // -- tortousity
      node = getUniqueElementByTagNames_(inode, "mechanical_properties", "tortuosity", flag);
      if (flag) {
        double val = GetAttributeValueD_(static_cast<DOMElement*>(node), "value");
        tmp_list.set<double>("aqueous tortuosity", val);
      }

      XMLString::release(&tagname);
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

      tagname = XMLString::transcode(inode->getNodeName());
      if (strcmp(tagname, "solute") != 0) continue;

      double val = GetAttributeValueD_(static_cast<DOMElement*>(inode), "coefficient_of_diffusion");
      text_content = XMLString::transcode(inode->getTextContent());

      aqueous_names.push_back(TrimString_(text_content));
      aqueous_values.push_back(val);

      XMLString::release(&text_content);
      XMLString::release(&tagname);
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

      tagname = XMLString::transcode(inode->getNodeName());
      if (strcmp(tagname, "solute") != 0) continue;

      double val = GetAttributeValueD_(static_cast<DOMElement*>(inode), "coefficient_of_diffusion");
      text_content = XMLString::transcode(inode->getTextContent());

      gaseous_names.push_back(TrimString_(text_content));
      gaseous_values.push_back(val);

      XMLString::release(&text_content);
      XMLString::release(&tagname);
    }

    diff_list.set<Teuchos::Array<std::string> >("gaseous names", gaseous_names);
    diff_list.set<Teuchos::Array<double> >("gaseous values", gaseous_values);
  }

    /*
      // now generate the source lists
      Teuchos::ParameterList src_list = CreateTransportSrcList_(plist);
      if (src_list.begin() != src_list.end()) { // the source lists are not empty
        out_list.sublist("source terms") = src_list;
      }

      // now generate the boundary conditions
      Teuchos::ParameterList& phase_list = plist->sublist("Phase Definitions");
      Teuchos::ParameterList& bc_list = plist->sublist("Boundary Conditions");

      for (Teuchos::ParameterList::ConstIterator i = bc_list.begin(); i != bc_list.end(); i++) {
        // read the assigned regions
        const std::string name(i->first);
        Teuchos::Array<std::string> regs = bc_list.sublist(name).get<Teuchos::Array<std::string> >("Assigned Regions");
        vv_bc_regions.insert(vv_bc_regions.end(), regs.begin(), regs.end());

        // only count sublists
        if (bc_list.isSublist(name)) {
          if (bc_list.sublist(name).isSublist("Solute BC")) {
            // read the solute bc stuff
            Teuchos::ParameterList& solbc = bc_list.sublist(name).sublist("Solute BC");
            for (Teuchos::ParameterList::ConstIterator ph = solbc.begin(); ph != solbc.end(); ph++) {
              std::string phase_name = ph->first; 
              int n;
              if (phase_name == "Aqueous") { 
                n = 0;
              } else if (phase_name == "Gaseous") {
                n = 1;
              } else {
                msg << "Boundary condition have unsupported phase " << phase_name;
                Exceptions::amanzi_throw(msg);
              }
              Teuchos::ParameterList& comps = solbc.sublist(phase_name).sublist(phases_[n].solute_name);

              for (std::vector<std::string>::iterator i = comp_names_all_.begin(); i != comp_names_all_.end(); i++) {
                if (comps.isSublist(*i)) {
                  std::stringstream compss;
                  compss << *i;
                  if (comps.sublist(*i).isSublist("BC: Uniform Concentration")) {
                    Teuchos::ParameterList& bcsub = comps.sublist(*i).sublist("BC: Uniform Concentration");

                    // If we only have one Geochemical condition, we don't need the Times and Time Functions
                    // entries, and make up entries for the 
                    if (bcsub.isParameter("Geochemical Condition")) { 
                      std::string cond_name = bcsub.get<std::string>("Geochemical Condition");

                      // Add an entry to Transport->boundary conditions->geochemical conditions.
                      Teuchos::ParameterList& gc_list = out_list.sublist("boundary conditions")
                          .sublist("geochemical conditions");
                      Teuchos::ParameterList& bc = gc_list.sublist(compss.str()).sublist(name);

                      // Fill it with made-up entries.
                      Teuchos::Array<double> times(2);
                      times[0] = -FLT_MAX;
                      times[1] = FLT_MAX;
                      Teuchos::Array<std::string> time_fns(1, "Constant");
                      Teuchos::Array<std::string> cond_names(2, cond_name);
                      bc.set<Teuchos::Array<double> >("Times", times);
                      bc.set<Teuchos::Array<std::string> >("Time Functions", time_fns);
                      bc.set<Teuchos::Array<std::string> >("regions", regs);
                      bc.set<Teuchos::Array<std::string> >("Geochemical Conditions", cond_names);
                    }

                    // Otherwise, we parse these entries.
                    else {
                      Teuchos::Array<double> times = bcsub.get<Teuchos::Array<double> >("Times");
                      Teuchos::Array<std::string> time_fns = bcsub.get<Teuchos::Array<std::string> >("Time Functions");

                      if (bcsub.isParameter("Geochemical Conditions")) { 
                        Teuchos::Array<std::string> cond_names = bcsub.get<Teuchos::Array<std::string> >("Geochemical Conditions");

                        // Add an entry to Transport->boundary conditions->geochemical conditions.
                        Teuchos::ParameterList& gc_list = out_list.sublist("boundary conditions")
                            .sublist("geochemical conditions");
                        Teuchos::ParameterList& bc = gc_list.sublist(compss.str()).sublist(name);

                        // Fill it with stuff.
                        bc.set<Teuchos::Array<double> >("Times", times);
                        bc.set<Teuchos::Array<std::string> >("Time Functions", time_fns);
                        bc.set<Teuchos::Array<std::string> >("Geochemical Conditions", cond_names);
                        bc.set<Teuchos::Array<std::string> >("regions", regs);
                      }
                      else { // ordinary Transport BCs.
                        Teuchos::ParameterList& tbc_list = out_list.sublist("boundary conditions").sublist("concentration");
                        Teuchos::ParameterList& bc = tbc_list.sublist(compss.str()).sublist(name);
                        bc.set<Teuchos::Array<std::string> >("regions",regs);
 
                        Teuchos::Array<double> values = bcsub.get<Teuchos::Array<double> >("Values");
                        Teuchos::ParameterList& bcfn = bc.sublist("boundary concentration").sublist("function-tabular");
                        bcfn.set<Teuchos::Array<double> >("y values", values);
                        bcfn.set<Teuchos::Array<double> >("x values", times);
                        bcfn.set<Teuchos::Array<std::string> >("forms", TranslateForms_(time_fns));
                      }
                    }
                  }
                }
              }
            }
          }
        }
       */

  // remaining global parameters
  out_list.set<int>("number of aqueous components", phases_["water"].size());
  out_list.set<int>("number of gaseous components", phases_["gas"].size());

  out_list.sublist("VerboseObject") = verb_list_.sublist("VerboseObject");
  return out_list;
}

}  // namespace AmanziInput
}  // namespace Amanzi


