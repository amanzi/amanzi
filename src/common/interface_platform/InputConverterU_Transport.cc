/*
  Input Converter

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Markus Berndt
           Erin Barker
           Jeffrey Johnson
           Konstantin Lipnikov (lipnikov@lanl.gov)
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
Teuchos::ParameterList InputConverterU::TranslateTransport_(const std::string& domain)
{
  Teuchos::ParameterList out_list, adv_list;

  if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH)
    *vo_->os() << "Translating transport, domain=" << domain << std::endl;

  MemoryManager mm;

  char *text, *tagname;
  DOMNodeList *node_list, *children;
  DOMNode* node;

  // create header
  out_list.set<std::string>("domain name", (domain == "matrix") ? "domain" : domain);

  // process CFL number
  bool flag;
  double cfl(1.0);

  node = GetUniqueElementByTagsString_("unstructured_controls, unstr_transport_controls, cfl", flag);
  if (flag) {
    text = mm.transcode(node->getTextContent());
    cfl = strtod(text, NULL);
  }

  // set defaults for transport
  out_list.set<int>("spatial discretization order", 1);
  out_list.set<int>("temporal discretization order", 1);
  out_list.set<double>("cfl", cfl);
  out_list.set<std::string>("flow mode", "transient");

  out_list.set<std::string>("solver", "Dispersion Solver");
  out_list.set<std::string>("preconditioner", LINEAR_SOLVER_PC);
  out_list.set<bool>("enable internal tests", false);
  out_list.set<bool>("transport subcycling", TRANSPORT_SUBCYCLING);

  // overwrite data from expert parameters  
  node = GetUniqueElementByTagsString_("unstructured_controls, unstr_transport_controls, sub_cycling", flag);
  if (flag) {
    text = mm.transcode(node->getTextContent());
    out_list.set<bool>("transport subcycling", (strcmp(text, "on") == 0));
  }

  int poly_order(0);
  std::string tags_default("unstructured_controls, unstr_transport_controls");
  node = GetUniqueElementByTagsString_(tags_default + ", algorithm", flag);
  if (flag) {
    std::string order = GetTextContentS_(node, "explicit first-order, explicit second-order, explicit, implicit");
    if (order == "explicit first-order") {
      out_list.set<int>("spatial discretization order", 1);
      out_list.set<int>("temporal discretization order", 1);
    } else if (order == "explicit second-order") {
      out_list.set<int>("spatial discretization order", 2);
      out_list.set<int>("temporal discretization order", 2);
      poly_order = 1;
    } else if (order == "explicit") {
      int nspace(-1), ntime(-1);
      node = GetUniqueElementByTagsString_(tags_default + ", spatial_order", flag);
      if (flag) nspace = std::strtol(mm.transcode(node->getTextContent()), NULL, 10);

      node = GetUniqueElementByTagsString_(tags_default + ", temporal_order", flag);
      if (flag) ntime = std::strtol(mm.transcode(node->getTextContent()), NULL, 10);

      if (nspace < 0 || nspace > 2 || ntime < 0 || ntime > 4) {
        Errors::Message msg;
        msg << "Transport algorithm=explicit requires spatial_order=1 or 2 "
            << "and temporal_order=1,2,3, or 4.\n";
        Exceptions::amanzi_throw(msg);
      }

      out_list.set<int>("spatial discretization order", nspace);
      out_list.set<int>("temporal discretization order", ntime);
      out_list.set<bool>("generic RK implementation", true);
      poly_order = 1;
    } else if (order == "implicit") {
      std::vector<std::string> dofs({"cell"});
      adv_list.sublist("matrix")
              .set<Teuchos::Array<std::string> >("schema", dofs)
              .set<int>("method order", 0)
              .set<std::string>("matrix type", "advection");
      poly_order = 1;
    }
  }

  // high-order transport
  // -- defaults
  Teuchos::ParameterList& trp_lift = out_list.sublist("reconstruction");
  trp_lift.set<int>("polynomial order", poly_order);
  trp_lift.set<std::string>("limiter", "tensorial");
  trp_lift.set<bool>("limiter extension for transport", true);
  trp_lift.set<std::string>("limiter stencil", "face to cells");

  // -- overwrite data from expert parameters  
  node = GetUniqueElementByTagsString_("unstructured_controls, unstr_transport_controls, limiter", flag);
  if (flag) {
    std::string limiter = GetTextContentS_(node, "tensorial, Kuzmin, Barth-Jespersen");
    trp_lift.set<std::string>("limiter", limiter);
  }

  node = GetUniqueElementByTagsString_("unstructured_controls, unstr_transport_controls, limiter_stencil", flag);
  if (flag) {
    std::string stencil = GetTextContentS_(node, "node-to-cells, face-to-cells, cell-to-closest-cells, cell-to-all-cells");
    std::replace(stencil.begin(), stencil.end(), '-', ' ');
    trp_lift.set<std::string>("limiter stencil", stencil);
  }

  // check if we need to write a dispersivity sublist
  node = doc_->getElementsByTagName(mm.transcode("materials"))->item(0);
  bool dispersion = (doc_->getElementsByTagName(mm.transcode("dispersion_tensor"))->getLength() > 0) ||
                    (doc_->getElementsByTagName(mm.transcode("tortuosity"))->getLength() > 0);

  // create dispersion list
  if (dispersion && domain == "matrix") {
    node_list = doc_->getElementsByTagName(mm.transcode("materials"));

    Teuchos::ParameterList& mat_list = out_list.sublist("material properties");

    children = node_list->item(0)->getChildNodes();
    for (int i = 0; i < children->getLength(); ++i) {
      DOMNode* inode = children->item(i);
      if (inode->getNodeType() != DOMNode::ELEMENT_NODE) continue;

      tagname = mm.transcode(inode->getNodeName());
      if (strcmp(tagname, "material") != 0) continue;

      // -- dispersion tensor
      Teuchos::ParameterList tmp_list;
      node = GetUniqueElementByTagsString_(inode, "mechanical_properties, dispersion_tensor", flag);
      if (flag) {
        double al, alh, alv, at, ath, atv;
        std::string model = GetAttributeValueS_(node, "type", "uniform_isotropic,burnett_frind,lichtner_kelkar_robinson");
        if (strcmp(model.c_str(), "uniform_isotropic") == 0) { 
          tmp_list.set<std::string>("model", "Bear");

          al = GetAttributeValueD_(node, "alpha_l", TYPE_NUMERICAL, 0.0, DVAL_MAX, "m");
          at = GetAttributeValueD_(node, "alpha_t", TYPE_NUMERICAL, 0.0, DVAL_MAX, "m");

          tmp_list.sublist("parameters for Bear").set<double>("alpha_l", al)
                                                 .set<double>("alpha_t", at);
        } else if (strcmp(model.c_str(), "burnett_frind") == 0) {
          tmp_list.set<std::string>("model", "Burnett-Frind");

          al = GetAttributeValueD_(node, "alpha_l", TYPE_NUMERICAL, 0.0, DVAL_MAX, "m");
          ath = GetAttributeValueD_(node, "alpha_th", TYPE_NUMERICAL, 0.0, DVAL_MAX, "m");
          atv = GetAttributeValueD_(node, "alpha_tv", TYPE_NUMERICAL, 0.0, DVAL_MAX, "m");

          tmp_list.sublist("parameters for Burnett-Frind")
              .set<double>("alpha_l", al).set<double>("alpha_th", ath)
              .set<double>("alpha_tv", atv);

          transport_permeability_ = true;
        } else if (strcmp(model.c_str(), "lichtner_kelkar_robinson") == 0) {
          tmp_list.set<std::string>("model", "Lichtner-Kelkar-Robinson");

          alh = GetAttributeValueD_(node, "alpha_lh", TYPE_NUMERICAL, 0.0, DVAL_MAX, "m");
          alv = GetAttributeValueD_(node, "alpha_lv", TYPE_NUMERICAL, 0.0, DVAL_MAX, "m");
          ath = GetAttributeValueD_(node, "alpha_th", TYPE_NUMERICAL, 0.0, DVAL_MAX, "m");
          atv = GetAttributeValueD_(node, "alpha_tv", TYPE_NUMERICAL, 0.0, DVAL_MAX, "m");

          tmp_list.sublist("parameters for Lichtner-Kelkar-Robinson")
              .set<double>("alpha_lh", alh).set<double>("alpha_lv", alv)
              .set<double>("alpha_th", ath).set<double>("alpha_tv", atv);

          transport_permeability_ = true;
        } 
      }

      // -- tortousity
      node = GetUniqueElementByTagsString_(inode, "mechanical_properties, tortuosity", flag);
      if (flag) {
        double val = GetAttributeValueD_(node, "value", TYPE_NUMERICAL, 0.0, DVAL_MAX);
        tmp_list.set<double>("aqueous tortuosity", val);
      }

      node = GetUniqueElementByTagsString_(inode, "mechanical_properties, tortuosity_gas", flag);
      if (flag) {
        double val = GetAttributeValueD_(node, "value", TYPE_NUMERICAL, 0.0, DVAL_MAX);
        tmp_list.set<double>("gaseous tortuosity", val);
      }

      if (tmp_list.numParams() > 0) { 
        node = GetUniqueElementByTagsString_(inode, "assigned_regions", flag);
        std::vector<std::string> regions = CharToStrings_(mm.transcode(node->getTextContent()));
        tmp_list.set<Teuchos::Array<std::string> >("regions", regions);

        std::string mat_name = GetAttributeValueS_(inode, "name");
        mat_list.sublist(mat_name) = tmp_list;
      }
    }
  }

  // -- molecular diffusion
  //    check in phases->water list (other solutes are ignored)
  //    check in phases->air list (other solutes are ignored)
  //    two names are supported (solutes and primaries)
  out_list.sublist("molecular diffusion") = TranslateMolecularDiffusion_();

  // add dispersion/diffusion operator
  node = GetUniqueElementByTagsString_(
      "unstructured_controls, unstr_transport_controls, dispersion_discretization_method", flag);
  std::string disc_methods;
  if (flag)
    disc_methods = mm.transcode(node->getTextContent());
  else
    disc_methods = (mesh_rectangular_) ? "mfd-monotone_for_hex" : "mfd-optimized_for_monotonicity";
  disc_methods.append(", mfd-two_point_flux_approximation");

  out_list.sublist("operators") = TranslateDiffusionOperator_(
      disc_methods, "diffusion_operator", "", "", "", false);

  // multiscale models sublist
  out_list.sublist("multiscale models") = TranslateTransportMSM_();
  if (out_list.sublist("multiscale models").numParams() > 0) {
    out_list.sublist("physical models and assumptions")
        .set<std::string>("multiscale model", "dual continuum discontinuous matrix");
  }

  // create the sources and boundary conditions lists
  out_list.sublist("boundary conditions") = TranslateTransportBCs_(domain);
  out_list.sublist("source terms") = TranslateTransportSources_();

  // remaining global parameters
  out_list.set<int>("number of aqueous components", phases_["water"].size());
  out_list.set<int>("number of gaseous components", phases_["air"].size());

  out_list.sublist("physical models and assumptions")
      .set<bool>("effective transport porosity", use_transport_porosity_);

  // cross coupling of PKs
  out_list.sublist("physical models and assumptions")
      .set<bool>("permeability field is required", transport_permeability_);

  if (fractures_ && domain != "domain") {
    out_list.sublist("physical models and assumptions").set<bool>("transport in fractures", true);
  }

  out_list.sublist("verbose object") = verb_list_.sublist("verbose object");

  // merging sublist
  if (adv_list.numParams() > 0)
    out_list.sublist("operators").sublist("advection operator") = adv_list;

  return out_list;
}


/* ******************************************************************
* Create list of molecular diffusion coefficients
****************************************************************** */
Teuchos::ParameterList InputConverterU::TranslateMolecularDiffusion_()
{
  Teuchos::ParameterList out_list;

  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH)
    *vo_->os() << "Translating molecular diffusion data" << std::endl;

  MemoryManager mm;

  bool flag;
  char *text, *tagname;
  DOMNodeList *children;
  DOMNode* node;

  std::vector<std::string> aqueous_names, gaseous_names;
  std::vector<double> aqueous_values, gaseous_values, molar_masses, henry_coef;

  // liquid phase
  std::string species("solute");
  node = GetUniqueElementByTagsString_("phases, liquid_phase, dissolved_components, solutes", flag);
  if (!flag) {
    node = GetUniqueElementByTagsString_("phases, liquid_phase, dissolved_components, primaries", flag);
    species = "primary";
  }

  if (flag) {
    children = node->getChildNodes();
    for (int i = 0; i < children->getLength(); ++i) {
      DOMNode* inode = children->item(i);
      if (inode->getNodeType() != DOMNode::ELEMENT_NODE) continue;

      tagname = mm.transcode(inode->getNodeName());
      if (strcmp(tagname, species.c_str()) != 0) continue;

      text = mm.transcode(inode->getTextContent());
      double val = GetAttributeValueD_(inode, "coefficient_of_diffusion", TYPE_NUMERICAL, 0.0, DVAL_MAX, "m^2/s", false);
      double mass = GetAttributeValueD_(inode, "molar_mass", TYPE_NUMERICAL, 0.0, DVAL_MAX, "kg/mol", false);

      aqueous_names.push_back(TrimString_(text));
      aqueous_values.push_back(val);
      molar_masses.push_back(mass);
    }
  }

  // gas phase
  species = "solute";
  node = GetUniqueElementByTagsString_("phases, gas_phase, dissolved_components, solutes", flag);
  if (!flag) {
    node = GetUniqueElementByTagsString_("phases, gas_phase, dissolved_components, primaries", flag);
    species = "primary";
  }

  if (flag) {
    children = node->getChildNodes();
    for (int i = 0; i < children->getLength(); ++i) {
      DOMNode* inode = children->item(i);
      if (inode->getNodeType() != DOMNode::ELEMENT_NODE) continue;

      tagname = mm.transcode(inode->getNodeName());
      if (strcmp(tagname, species.c_str()) != 0) continue;

      text = mm.transcode(inode->getTextContent());
      double val = GetAttributeValueD_(inode, "coefficient_of_diffusion", TYPE_NUMERICAL, 0.0, DVAL_MAX, "", false);
      double kh = GetAttributeValueD_(inode, "kh", TYPE_NUMERICAL, 0.0, DVAL_MAX);

      gaseous_names.push_back(TrimString_(text));
      gaseous_values.push_back(val);
      henry_coef.push_back(kh);
    }
  }

  out_list.set<Teuchos::Array<std::string> >("aqueous names", aqueous_names);
  out_list.set<Teuchos::Array<double> >("aqueous values", aqueous_values);
  out_list.set<Teuchos::Array<double> >("molar masses", molar_masses);
  if (gaseous_names.size() > 0) {
    out_list.set<Teuchos::Array<std::string> >("gaseous names", gaseous_names);
    out_list.set<Teuchos::Array<double> >("gaseous values", gaseous_values);
    out_list.set<Teuchos::Array<double> >("air-water partitioning coefficient", henry_coef);
    out_list.set<Teuchos::Array<double> >("Henry dimensionless constants", henry_coef);
  }

  return out_list;
}


/* ******************************************************************
* Create list of multiscale models.
****************************************************************** */
Teuchos::ParameterList InputConverterU::TranslateTransportMSM_()
{
  Teuchos::ParameterList out_list;

  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH)
    *vo_->os() << "Translating multiscale models" << std::endl;

  MemoryManager mm;
  DOMNodeList *node_list, *children;
  DOMNode *node;
  DOMElement* element;

  bool flag;
  std::string name, model, rel_perm;

  // check that list of models does exist
  node_list = doc_->getElementsByTagName(mm.transcode("materials_secondary_continuum"));
  if (node_list->getLength() == 0) return out_list;

  element = static_cast<DOMElement*>(node_list->item(0));
  children = element->getElementsByTagName(mm.transcode("material"));

  for (int i = 0; i < children->getLength(); ++i) {
    DOMNode* inode = children->item(i); 

    node = GetUniqueElementByTagsString_(inode, "multiscale_model", flag);
    if (!flag) ThrowErrorMissing_("materials", "element", "multiscale_model", "material");
    model = GetAttributeValueS_(node, "name");

    // verify material name
    name = GetAttributeValueS_(inode, "name");
    if (! FindNameInVector_(name, material_names_)) {
      Errors::Message msg("Materials names in primary and secondary continua do not match.\n");
      Exceptions::amanzi_throw(msg);
    } 

    // common for all models
    std::stringstream ss;
    ss << "MSM " << i;
    Teuchos::ParameterList& msm_plist = out_list.sublist(ss.str());
    Teuchos::ParameterList& msm_slist = msm_plist.sublist(model + " parameters");

    // -- regions
    node = GetUniqueElementByTagsString_(inode, "assigned_regions", flag);
    std::vector<std::string> regions = CharToStrings_(mm.transcode(node->getTextContent()));
    msm_plist.set<Teuchos::Array<std::string> >("regions", regions)
             .set<std::string>("multiscale model", model);

    // -- volume fraction
    node = GetUniqueElementByTagsString_(inode, "volume_fraction", flag);
    double vof = GetTextContentD_(node, "", true);

    // -- tortousity
    double tau(1.0);
    node = GetUniqueElementByTagsString_(inode, "mechanical_properties, tortuosity", flag);
    if (flag) tau = GetAttributeValueD_(node, "value", TYPE_NUMERICAL, 0.0, DVAL_MAX, "-");
    msm_slist.set<double>("matrix tortuosity", tau);

    // dual porosity model
    if (model == "dual porosity") {
      node = GetUniqueElementByTagsString_(inode, "multiscale_model, matrix", flag);
      if (!flag) ThrowErrorMissing_("materials", "element", "matrix", "multiscale_model");

      double depth = GetAttributeValueD_(node, "depth", TYPE_NUMERICAL, 0.0, DVAL_MAX, "m");
      double sigma = GetAttributeValueD_(node, "warren_root", TYPE_NUMERICAL, 0.0, 1000.0, "m^-2");

      msm_slist.set<double>("Warren Root parameter", sigma)
               .set<double>("matrix depth", depth)
               .set<double>("matrix volume fraction", vof);

    } else if (model == "generalized dual porosity") {
      node = GetUniqueElementByTagsString_(inode, "multiscale_model, matrix", flag);
      if (!flag) ThrowErrorMissing_("materials", "element", "matrix", "multiscale_model");

      int nnodes = GetAttributeValueL_(node, "number_of_nodes", TYPE_NUMERICAL, 0, INT_MAX, false, 1);
      double depth = GetAttributeValueD_(node, "depth", TYPE_NUMERICAL, 0.0, DVAL_MAX, "m");
    
      msm_slist.set<int>("number of matrix nodes", nnodes)
               .set<double>("matrix depth", depth)
               .set<double>("matrix volume fraction", vof);
    }
  }

  return out_list;
}


/* ******************************************************************
* Create list of transport BCs.
****************************************************************** */
Teuchos::ParameterList InputConverterU::TranslateTransportBCs_(const std::string& domain)
{
  Teuchos::ParameterList out_list;

  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH)
      *vo_->os() << "Translating boundary conditions" << std::endl;

  MemoryManager mm;

  char *text;
  DOMNodeList *children;
  DOMNode *node, *phase;
  DOMElement* element;

  bool flag;
  if (domain == "matrix")
    node = GetUniqueElementByTagsString_("boundary_conditions", flag);
  else
    node = GetUniqueElementByTagsString_("fracture_network, boundary_conditions", flag);
  if (!flag) return out_list;

  children = node->getChildNodes();
  int nchildren = children->getLength();
  if (nchildren == 0) return out_list;

  for (int i = 0; i < nchildren; ++i) {
    DOMNode* inode = children->item(i);
    if (inode->getNodeType() != DOMNode::ELEMENT_NODE) continue;
    std::string bcname = GetAttributeValueS_(inode, "name");

    // read the assigned regions
    node = GetUniqueElementByTagsString_(inode, "assigned_regions", flag);
    text = mm.transcode(node->getTextContent());
    std::vector<std::string> regions = CharToStrings_(text);

    vv_bc_regions_.insert(vv_bc_regions_.end(), regions.begin(), regions.end());

    // process different phases
    // -- liquid phase
    phase = GetUniqueElementByTagsString_(inode, "liquid_phase", flag);
    if (flag) {
      element = static_cast<DOMElement*>(phase);
      DOMNodeList* solutes = element->getElementsByTagName(mm.transcode("solute_component"));
      TranslateTransportBCsGroup_(bcname, regions, solutes, out_list);
    }

    // -- gas phase
    phase = GetUniqueElementByTagsString_(inode, "gas_phase", flag);
    if (flag) {
      element = static_cast<DOMElement*>(phase);
      DOMNodeList* solutes = element->getElementsByTagName(mm.transcode("solute_component"));
      TranslateTransportBCsGroup_(bcname, regions, solutes, out_list);
    }

    // geochemical BCs 
    node = GetUniqueElementByTagsString_(inode, "liquid_phase, geochemistry_component", flag);
    if (flag) {
      TranslateTransportGeochemistry_(node, bcname, regions, out_list);
    }
  }

  // backward compatibility: translate constraints for native chemistry
  TranslateTransportBCsAmanziGeochemistry_(out_list);

  return out_list;
}


/* ******************************************************************
* Create list of transport BCs for particular group of solutes.
* Solutes may have only one element, see schema for details.
****************************************************************** */
void InputConverterU::TranslateTransportBCsGroup_(
    std::string& bcname, std::vector<std::string>& regions,
    DOMNodeList* solutes, Teuchos::ParameterList& out_list)
{
  if (solutes->getLength() == 0) return;
 
  DOMNode* node = solutes->item(0);

  // get child nodes with the same tagname
  bool flag;
  std::string bctype, solute_name, tmp_name, unit("molar");
  std::vector<DOMNode*> same_list = GetSameChildNodes_(node, bctype, flag, true);

  while (same_list.size() > 0) {
    // process a group of elements named after the 0-th element
    solute_name = GetAttributeValueS_(same_list[0], "name");

    // Check for spatially dependent BCs. Only one is allowed (FIXME)
    bool space_bc, time_bc;
    DOMNode* knode = GetUniqueElementByTagsString_(same_list[0], "space", space_bc);
    DOMNode* lnode = GetUniqueElementByTagsString_(same_list[0], "time", time_bc);

    std::string filename, xheader, yheader;
    std::vector<double> times, values;
    std::vector<double> data, data_tmp;
    std::vector<std::string> forms, formulas;

    if (space_bc) {
      std::string tmp = GetAttributeValueS_(knode, "amplitude");
      data.push_back(ConvertUnits_(tmp, unit, solute_molar_mass_[solute_name]));

      data_tmp = GetAttributeVectorD_(knode, "center", -1, "m");
      data.insert(data.end(), data_tmp.begin(), data_tmp.end());
      data.push_back(GetAttributeValueD_(knode, "standard_deviation", TYPE_NUMERICAL, 0.0, DVAL_MAX, "m"));

      if (time_bc) {
        data[0] *= GetAttributeValueD_(lnode, "data", TYPE_NUMERICAL, DVAL_MIN, DVAL_MAX, "");
      }

      same_list.erase(same_list.begin());
    } else {
      DOMElement* element = static_cast<DOMElement*>(same_list[0]);
      filename = GetAttributeValueS_(element, "h5file", TYPE_NONE, false, "");
      if (filename != "") {
        xheader = GetAttributeValueS_(element, "times", TYPE_NONE);
        yheader = GetAttributeValueS_(element, "values", TYPE_NONE);

        same_list.erase(same_list.begin());
      } else {
        std::map<double, double> tp_values;
        std::map<double, std::string> tp_forms, tp_formulas;

        for (std::vector<DOMNode*>::iterator it = same_list.begin(); it != same_list.end(); ++it) {
          tmp_name = GetAttributeValueS_(*it, "name");

          if (tmp_name == solute_name) {
            double t0 = GetAttributeValueD_(*it, "start", TYPE_TIME, DVAL_MIN, DVAL_MAX, "s");
            tp_forms[t0] = GetAttributeValueS_(*it, "function");
            GetAttributeValueD_(*it, "value", TYPE_NUMERICAL, DVAL_MIN, DVAL_MAX, "molar");  // just a check
            tp_values[t0] = ConvertUnits_(GetAttributeValueS_(*it, "value"), unit, solute_molar_mass_[solute_name]);
            tp_formulas[t0] = GetAttributeValueS_(*it, "formula", TYPE_NONE, false, "");

            same_list.erase(it);
            it--;
          } 
        }

        // create vectors of values and forms
        for (std::map<double, double>::iterator it = tp_values.begin(); it != tp_values.end(); ++it) {
          times.push_back(it->first);
          values.push_back(it->second);
          forms.push_back(tp_forms[it->first]);
          formulas.push_back(tp_formulas[it->first]);
        }
      }
    }
     
    // save in the XML files  
    Teuchos::ParameterList& tbc_list = out_list.sublist("concentration");
    Teuchos::ParameterList& bc = tbc_list.sublist(solute_name).sublist(bcname);
    bc.set<Teuchos::Array<std::string> >("regions", regions);

    Teuchos::ParameterList& bcfn = bc.sublist("boundary concentration");
    if (space_bc) {
      TranslateFunctionGaussian_(data, bcfn);
    } else if (filename != "") {
      bcfn.sublist("function-tabular")
          .set<std::string>("file", filename)
          .set<std::string>("x header", xheader)
          .set<std::string>("y header", yheader);
    } else {
      TranslateGenericMath_(times, values, forms, formulas, bcfn);
    }

    // data distribution method
    bc.set<std::string>("spatial distribution method", "none");
    bc.set<bool>("use area fractions", WeightVolumeSubmodel_(regions));
  }
}


/* ******************************************************************
* Create list of transport BCs and sources for geochemistry.
****************************************************************** */
void InputConverterU::TranslateTransportGeochemistry_(
    DOMNode* node, std::string& bcname, std::vector<std::string>& regions,
    Teuchos::ParameterList& out_list)
{
  bool flag;
  std::string bctype;
  std::vector<DOMNode*> same_list = GetSameChildNodes_(node, bctype, flag, true);
  std::map<double, std::string> tp_forms, tp_values;

  for (int j = 0; j < same_list.size(); ++j) {
    DOMElement* element = static_cast<DOMElement*>(same_list[j]);
    double t0 = GetAttributeValueD_(element, "start", TYPE_TIME, DVAL_MIN, DVAL_MAX, "s");
    tp_values[t0] = GetAttributeValueS_(element, "name");
    tp_forms[t0] = GetAttributeValueS_(element, "function", TYPE_NONE, false, "constant");
    // no form -> use geochemistry engine
  }

  // create vectors of values and forms
  std::vector<double> times;
  std::vector<std::string> forms, values;
  for (std::map<double, std::string>::iterator it = tp_values.begin(); it != tp_values.end(); ++it) {
    times.push_back(it->first);
    values.push_back(it->second);
    forms.push_back(tp_forms[it->first]);
  }

  if (times.size() == 1) {
    times.push_back(times[0] + 1e+20);
    values.push_back(values[0]);
  } else {
    forms.pop_back();
  }

  // save in the XML files  
  Teuchos::ParameterList& tbc_list = out_list.sublist("geochemical");
  Teuchos::Array<std::string> solute_names;
  for (int i = 0; i < phases_["water"].size(); ++i) {
    solute_names.push_back(phases_["water"][i]);
  }
  Teuchos::ParameterList& bc = tbc_list.sublist(bcname);
  bc.set<Teuchos::Array<std::string> >("regions", regions)
    .set<Teuchos::Array<std::string> >("solutes", solute_names)
    .set<Teuchos::Array<double> >("times", times)
    .set<Teuchos::Array<std::string> >("geochemical conditions", values)
    .set<Teuchos::Array<std::string> >("time functions", forms);
}


/* ******************************************************************
* Create list of transport BCs for native chemistry.
* Solutes may have only one element, see schema for details.
****************************************************************** */
void InputConverterU::TranslateTransportBCsAmanziGeochemistry_(
    Teuchos::ParameterList& out_list)
{
  if (out_list.isSublist("geochemical") &&
      pk_model_["chemistry"] == "amanzi") {

    bool flag;
    std::string name, bc_name;
    DOMNode* node;
    DOMElement* element;

    node = GetUniqueElementByTagsString_("geochemistry, constraints", flag);

    Teuchos::ParameterList& bc_new = out_list.sublist("constraints");
    Teuchos::ParameterList& bc_old = out_list.sublist("geochemical");

    for (auto it = bc_old.begin(); it != bc_old.end(); ++it) {
      name = it->first;
 
      Teuchos::ParameterList& bco = bc_old.sublist(name);
      Teuchos::ParameterList& bcn = bc_new.sublist(name);

      std::vector<std::string> solutes = bco.get<Teuchos::Array<std::string> >("solutes").toVector();
      int nsolutes = solutes.size();
      Teuchos::Array<std::string> types(nsolutes);

      bcn.sublist("boundary constraints").set<int>("number of dofs", nsolutes);

      for (int n = 0; n < nsolutes; ++n) {
        std::stringstream dof;
        dof << "dof " << n + 1 << " function";
        Teuchos::ParameterList& fnc = bcn.sublist("boundary constraints").sublist(dof.str()).sublist("function-tabular");

        bcn.set("regions", bco.get<Teuchos::Array<std::string> >("regions"))
           .set<std::string>("spatial distribution method", "none")
           .set<bool>("use area fractions", false);
        fnc.set("x values", bco.get<Teuchos::Array<double> >("times"));
        fnc.set("forms", bco.get<Teuchos::Array<std::string> >("time functions"));
      
        // convert constraints to values
        Teuchos::Array<double> values;
        std::vector<std::string> constraints = 
            bco.get<Teuchos::Array<std::string> >("geochemical conditions").toVector();

        for (int i = 0; i < constraints.size(); ++i) {
          element = GetUniqueChildByAttribute_(node, "name", constraints[i], flag, true);
          element = GetUniqueChildByAttribute_(element, "name", solutes[n], flag, true);
          values.push_back(GetAttributeValueD_(element, "value"));
          types[n] = GetAttributeValueS_(element, "type");  // FIXME
        }
        fnc.set("y values", values);
      }
      bcn.set<Teuchos::Array<std::string> >("names", types);
    }

    out_list.remove("geochemical");
  }
}


/* ******************************************************************
* Create list of transport sources.
****************************************************************** */
Teuchos::ParameterList InputConverterU::TranslateTransportSources_()
{
  Teuchos::ParameterList out_list;

  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH)
      *vo_->os() << "Translating source terms" << std::endl;

  MemoryManager mm;

  char *text;
  DOMNodeList *node_list, *children;
  DOMNode *node;
  DOMElement* element;

  node_list = doc_->getElementsByTagName(mm.transcode("sources"));
  if (node_list->getLength() == 0) return out_list;

  children = node_list->item(0)->getChildNodes();

  for (int i = 0; i < children->getLength(); ++i) {
    DOMNode* inode = children->item(i);
    if (inode->getNodeType() != DOMNode::ELEMENT_NODE) continue;
    std::string srcname = GetAttributeValueS_(inode, "name");

    // read the assigned regions
    bool flag;
    node = GetUniqueElementByTagsString_(inode, "assigned_regions", flag);
    text = mm.transcode(node->getTextContent());
    std::vector<std::string> regions = CharToStrings_(text);

    vv_src_regions_.insert(vv_src_regions_.end(), regions.begin(), regions.end());

    // process different phases
    // -- liquid phase
    DOMNode* phase_l = GetUniqueElementByTagsString_(inode, "liquid_phase", flag);
    if (flag) {
      element = static_cast<DOMElement*>(phase_l);
      DOMNodeList* solutes = element->getElementsByTagName(mm.transcode("solute_component"));
      TranslateTransportSourcesGroup_(srcname, regions, solutes, phase_l, out_list);
    }

    // -- gas phase
    DOMNode* phase_g = GetUniqueElementByTagsString_(inode, "gas_phase", flag);
    if (flag) {
      element = static_cast<DOMElement*>(phase_g);
      DOMNodeList* solutes = element->getElementsByTagName(mm.transcode("solute_component"));
      TranslateTransportSourcesGroup_(srcname, regions, solutes, phase_l, out_list);
    }

    // geochemical sources
    node = GetUniqueElementByTagsString_(inode, "liquid_phase, geochemistry_component", flag);
    if (flag) {
      TranslateTransportGeochemistry_(node, srcname, regions, out_list);
    }
  }

  return out_list;
}


/* ******************************************************************
* Create list of transport sources.
****************************************************************** */
void InputConverterU::TranslateTransportSourcesGroup_(
    std::string& srcname, std::vector<std::string>& regions,
    DOMNodeList* solutes, DOMNode* phase_l, Teuchos::ParameterList& out_list)
{
  MemoryManager mm;
  DOMNodeList* node_list;
  DOMNode* node;
  DOMElement* element;

  for (int n = 0; n < solutes->getLength(); ++n) {
    node = solutes->item(n);

    // get a group of similar elements defined by the first element
    bool flag;
    std::string srctype, solute_name, srctype_flow, unit("mol/s"), unit_in;

    std::vector<DOMNode*> same_list = GetSameChildNodes_(node, srctype, flag, true);
    solute_name = GetAttributeValueS_(same_list[0], "name");

    // weighting method
    bool classical(true), mass_fraction(false);
    std::string weight("volume");

    char* text = mm.transcode(same_list[0]->getNodeName());
    if (strcmp(text, "perm_weighted") == 0) {
      weight = "permeability";
    } else if (strcmp(text, "volume_weighted") == 0) {
      weight = "volume";
    } else if (strcmp(text, "uniform_conc") == 0) {
      weight = "none";
      unit = "mol/s/m^3";
    } else if (strcmp(text, "flow_weighted_conc") == 0) {
      element = static_cast<DOMElement*>(phase_l);
      node_list = element->getElementsByTagName(mm.transcode("liquid_component")); 
      if (node_list->getLength() == 0)
        ThrowErrorIllformed_(srcname, "liquid_component", text);

      GetSameChildNodes_(node_list->item(0), srctype_flow, flag, true);
      weight = (srctype_flow == "volume_weighted") ? "volume" : "permeability";
    } else if (strcmp(text, "flow_mass_fraction_conc") == 0) {
      mass_fraction = true;
      unit = "-";
    } else if (strcmp(text, "diffusion_dominated_release") == 0) {
      classical = false;
    } else {
      ThrowErrorIllformed_(srcname, "element", text);
    } 
    if (weight == "permeability") transport_permeability_ = true;

    if (classical) {
      std::map<double, double> tp_values;
      std::map<double, std::string> tp_forms;

      for (int j = 0; j < same_list.size(); ++j) {
        element = static_cast<DOMElement*>(same_list[j]);
        double t0 = GetAttributeValueD_(element, "start", TYPE_TIME, DVAL_MIN, DVAL_MAX, "s");
        tp_forms[t0] = GetAttributeValueS_(element, "function");
        tp_values[t0] = ConvertUnits_(GetAttributeValueS_(element, "value"),
                                      unit_in, solute_molar_mass_[solute_name]);
        GetAttributeValueD_(element, "value", TYPE_NUMERICAL, DVAL_MIN, DVAL_MAX, unit);  // unit check

        // ugly correction when liquid/solute lists match
        if (mass_fraction) {
          element = static_cast<DOMElement*>(phase_l);
          node_list = element->getElementsByTagName(mm.transcode("liquid_component")); 
          std::vector<DOMNode*> tmp_list = GetSameChildNodes_(node_list->item(0), srctype_flow, flag, true);

          if (tmp_list.size() != same_list.size())
            ThrowErrorIllformed_(srcname, "liquid_component", text);

          double val = GetAttributeValueD_(tmp_list[j], "value");
          if (val > 0) {
            tp_values[t0] *= val / solute_molar_mass_[solute_name];
          } else {
            tp_values[t0] = val / rho_;
          }
        }
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

      // producer or injector?
      int well(0);
      for (int j = 0; j < values.size(); ++j) {
        if (values[j] < 0) well |= 1;
        if (values[j] > 0) well |= 2;
      }
      if (well == 3) {
        Errors::Message msg;
        msg << "Source \"" << srcname << "\" changes sign.\n";
        Exceptions::amanzi_throw(msg);
      }

      // save in the XML files  
      Teuchos::ParameterList& src_list = out_list.sublist("concentration");
      Teuchos::ParameterList& src = src_list.sublist(solute_name).sublist(srcname);
      src.set<Teuchos::Array<std::string> >("regions", regions);
      src.set<std::string>("spatial distribution method", weight);
      src.set<bool>("use volume fractions", WeightVolumeSubmodel_(regions));

      Teuchos::ParameterList& srcfn = src.sublist((well == 1) ? "producer" : "injector");
      if (times.size() == 1) {
        srcfn.sublist("function-constant").set<double>("value", values[0]);
      } else {
        srcfn.sublist("function-tabular")
            .set<Teuchos::Array<double> >("x values", times)
            .set<Teuchos::Array<double> >("y values", values)
            .set<Teuchos::Array<std::string> >("forms", forms);
      }
    } else {
      element = static_cast<DOMElement*>(same_list[0]);
      double total = GetAttributeValueD_(element, "total_inventory", TYPE_NUMERICAL, 0.0, DVAL_MAX, "mol");
      double diff = GetAttributeValueD_(element, "effective_diffusion_coefficient", TYPE_NUMERICAL, 0.0, DVAL_MAX);
      double length = GetAttributeValueD_(element, "mixing_length", TYPE_NUMERICAL, 0.0, DVAL_MAX, "m");
      std::vector<double> times;
      times.push_back(GetAttributeValueD_(element, "start", TYPE_TIME, DVAL_MIN, DVAL_MAX, "s"));

      element = static_cast<DOMElement*>(same_list[1]);
      times.push_back(GetAttributeValueD_(element, "start", TYPE_TIME, DVAL_MIN, DVAL_MAX, "s"));

      // save data in the XML
      Teuchos::ParameterList& src_list = out_list.sublist("concentration");
      Teuchos::ParameterList& src = src_list.sublist(solute_name).sublist(srcname);
      src.set<Teuchos::Array<std::string> >("regions", regions);
      src.set<std::string>("spatial distribution method", weight);
      src.set<bool>("use volume fractions", WeightVolumeSubmodel_(regions));

      std::vector<double> values(2, 0.0);
      std::vector<std::string> forms(1, "SQRT");
      double amplitude = 2 * total / length * std::pow(diff / M_PI, 0.5); 

      Teuchos::ParameterList& srcfn = src.sublist("well").sublist("function-tabular");
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


/* ******************************************************************
* Setup flag for volume fractions.
****************************************************************** */
bool InputConverterU::WeightVolumeSubmodel_(const std::vector<std::string>& regions)
{
  bool flag(false);
  for (int k = 0; k < regions.size(); ++k) {
    if (region_type_[regions[k]] == 1) flag = true;
  }
  return flag;
}

}  // namespace AmanziInput
}  // namespace Amanzi


