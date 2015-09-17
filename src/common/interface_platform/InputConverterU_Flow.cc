/*
  This is the input component of the Amanzi code. 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Markus Berndt (original version)
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
* Create flow list.
****************************************************************** */
Teuchos::ParameterList InputConverterU::TranslateFlow_(const std::string& mode)
{
  Teuchos::ParameterList out_list;
  Teuchos::ParameterList* flow_list;

  if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH)
    *vo_->os() << "Translating flow, mode=" << mode << std::endl;

  MemoryManager mm;
  DOMNode* node;

  // set up default values for some expert parameters
  double atm_pres(ATMOSPHERIC_PRESSURE);
  std::string rel_perm("upwind: amanzi");  
  std::string update_upwind("every timestep");

  // process expert parameters
  bool flag;
  node = GetUniqueElementByTagsString_("unstructured_controls, unstr_flow_controls, rel_perm_method", flag);
  if (flag) rel_perm = mm.transcode(node->getNodeName());

  // create flow header
  if (pk_model_["flow"] == "darcy") {
    Teuchos::ParameterList& darcy_list = out_list.sublist("Darcy problem");
    darcy_list.set<double>("atmospheric pressure", atm_pres);

    flow_list = &darcy_list;
    flow_single_phase_ = true;
  } else if (pk_model_["flow"] == "richards") {
    Teuchos::ParameterList& richards_list = out_list.sublist("Richards problem");
    Teuchos::ParameterList& upw_list = richards_list.sublist("upwind");
    upw_list.set<std::string>("relative permeability", rel_perm);
    upw_list.set<std::string>("upwind update", update_upwind);

    // "standard" is the most robust upwind method for variety of subsurface
    // scenarios. Note that "upwind: amanzi" requires "upwind method"="divk" 
    // to reproduce the same behavior on orthogonal meshes. 
    if (strcmp(rel_perm.c_str(), "upwind: amanzi") == 0) {
      upw_list.set<std::string>("upwind method", "divk");
      upw_list.sublist("upwind divk parameters").set<double>("tolerance", 1e-12);
    } else {
      upw_list.set<std::string>("upwind method", "standard");
      upw_list.sublist("upwind standard parameters").set<double>("tolerance", 1e-12);
    }
    richards_list.set<double>("atmospheric pressure", atm_pres);
    flow_list = &richards_list;

    richards_list.sublist("water retention models") = TranslateWRM_();
    richards_list.sublist("porosity models") = TranslatePOM_();
    if (richards_list.sublist("porosity models").numParams() > 0) {
      flow_list->sublist("physical models and assumptions")
          .set<std::string>("porosity model", "compressible: pressure function");
    }
    richards_list.sublist("multiscale models") = TranslateFlowMSM_();
    if (richards_list.sublist("multiscale models").numParams() > 0) {
      flow_list->sublist("physical models and assumptions")
          .set<std::string>("multiscale model", "dual porosity");
    }
  } else {
    Errors::Message msg;
    msg << "Internal error for flow model \"" << pk_model_["flow"] << "\".\n";
    Exceptions::amanzi_throw(msg);
  }

  // insert operator sublist
  std::string disc_method("mfd-optimized_for_sparsity");
  node = GetUniqueElementByTagsString_("unstructured_controls, unstr_flow_controls, discretization_method", flag);
  if (flag) disc_method = mm.transcode(node->getTextContent());

  std::string pc_method("linearized_operator");
  node = GetUniqueElementByTagsString_("unstructured_controls, unstr_flow_controls, preconditioning_strategy", flag);
  if (flag) pc_method = GetTextContentS_(node, "linearized_operator, diffusion_operator"); 

  std::string nonlinear_solver("nka");
  node = GetUniqueElementByTagsString_("unstructured_controls, unstr_nonlinear_solver", flag);
  if (flag) nonlinear_solver = GetAttributeValueS_(static_cast<DOMElement*>(node), "name", false, "nka"); 

  bool modify_correction(false);
  node = GetUniqueElementByTagsString_("unstructured_controls, unstr_nonlinear_solver, modify_correction", flag);

  // Newton method requires to overwrite other parameters.
  if (nonlinear_solver == "newton") {
    modify_correction = true;  // a hack
    if (disc_method != "fv-default" ||
        rel_perm != "upwind-darcy_velocity" ||
        update_upwind != "every nonlinear iteration" ||
        !modify_correction) {
      Errors::Message msg;
      msg << "Nonlinear solver \"newton\" requires \"upwind-darcy_velocity\"\n";
      msg << "\"modify_correction\"=true, and discretization method \"fv-default\".\n";
      Exceptions::amanzi_throw(msg);
    }
  }

  flow_list->sublist("operators") = TranslateDiffusionOperator_(
      disc_method, pc_method, nonlinear_solver, "vapor matrix");
  
  // insert time integrator
  std::string err_options, unstr_controls;
  if (mode == "steady") {
    err_options = "pressure";
    unstr_controls = "unstructured_controls, unstr_steady-state_controls";
  } else {
    err_options = "pressure, residual";
    unstr_controls = "unstructured_controls, unstr_transient_controls";

    // restart leads to a conflict
    node = GetUniqueElementByTagsString_(unstr_controls + ", unstr_initialization", flag); 
    if (flag && restart_) {
      Errors::Message msg;
      msg << "Parameters \"restart\" and \"unstr_transient_control->unstr_initialization\""
          << " are mutually exclusive.\n";
      Exceptions::amanzi_throw(msg);
    }
  } 
  
  if (pk_master_.find("flow") != pk_master_.end()) {
    flow_list->sublist("time integrator") = TranslateTimeIntegrator_(
        err_options, nonlinear_solver, modify_correction, unstr_controls,
        dt_cut_[mode], dt_inc_[mode]);
  }

  // insert boundary conditions and source terms
  flow_list->sublist("boundary conditions") = TranslateFlowBCs_();
  flow_list->sublist("source terms") = TranslateFlowSources_();

  // models and default assumptions. 
  // Note that MPC/PKs may overwrite these parameters
  flow_list->sublist("physical models and assumptions")
      .set<std::string>("water content model", "constant density");

  flow_list->sublist("VerboseObject") = verb_list_.sublist("VerboseObject");
  return out_list;
}


/* ******************************************************************
* Create list of water retention models.
****************************************************************** */
Teuchos::ParameterList InputConverterU::TranslateWRM_()
{
  Teuchos::ParameterList out_list;

  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH)
    *vo_->os() << "Translating water retension models" << std::endl;

  MemoryManager mm;
  DOMNodeList *node_list, *children;
  DOMNode* node;
  DOMElement* element;

  bool flag;
  std::string model, rel_perm;

  node_list = doc_->getElementsByTagName(mm.transcode("materials"));
  element = static_cast<DOMElement*>(node_list->item(0));
  children = element->getElementsByTagName(mm.transcode("material"));

  for (int i = 0; i < children->getLength(); ++i) {
    DOMNode* inode = children->item(i); 

    node = GetUniqueElementByTagsString_(inode, "cap_pressure", flag);
    model = GetAttributeValueS_(static_cast<DOMElement*>(node), "model", "van_genuchten, brooks_corey");
    DOMNode* nnode = GetUniqueElementByTagsString_(node, "parameters", flag);
    DOMElement* element_cp = static_cast<DOMElement*>(nnode);

    node = GetUniqueElementByTagsString_(inode, "rel_perm", flag);
    rel_perm = GetAttributeValueS_(static_cast<DOMElement*>(node), "model", "mualem, burdine");
    DOMNode* mnode = GetUniqueElementByTagsString_(node, "exp", flag);
    DOMElement* element_rp = (flag) ? static_cast<DOMElement*>(mnode) : NULL;

    // common stuff
    // -- assigned regions
    node = GetUniqueElementByTagsString_(inode, "assigned_regions", flag);
    std::vector<std::string> regions = CharToStrings_(mm.transcode(node->getTextContent()));

    // -- smoothing
    double krel_smooth = GetAttributeValueD_(element_cp, "optional_krel_smoothing_interval", false, 0.0);
    if (krel_smooth < 0.0) {
      Errors::Message msg;
      msg << "value of optional_krel_smoothing_interval must be non-negative.\n";
      Exceptions::amanzi_throw(msg);
    }

    // -- ell
    double ell, ell_d = (rel_perm == "mualem") ? ELL_MUALEM : ELL_BURDINE;
    ell = GetAttributeValueD_(element_rp, "value", false, ell_d);

    std::replace(rel_perm.begin(), rel_perm.begin() + 1, 'm', 'M');
    std::replace(rel_perm.begin(), rel_perm.begin() + 1, 'b', 'B');

    if (strcmp(model.c_str(), "van_genuchten") == 0) {
      double alpha = GetAttributeValueD_(element_cp, "alpha");
      double sr = GetAttributeValueD_(element_cp, "sr");
      double m = GetAttributeValueD_(element_cp, "m");

      for (std::vector<std::string>::const_iterator it = regions.begin(); it != regions.end(); ++it) {
        std::stringstream ss;
        ss << "WRM for " << *it;

        Teuchos::ParameterList& wrm_list = out_list.sublist(ss.str());

        wrm_list.set<std::string>("water retention model", "van Genuchten")
            .set<std::string>("region", *it)
            .set<double>("van Genuchten m", m)
            .set<double>("van Genuchten l", ell)
            .set<double>("van Genuchten alpha", alpha)
            .set<double>("residual saturation", sr)
            .set<double>("regularization interval", krel_smooth)
            .set<std::string>("relative permeability model", rel_perm);
      
        if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH) {
          Teuchos::ParameterList& file_list = wrm_list.sublist("output");
          std::stringstream name;
          name << *it << ".txt";
          file_list.set<std::string>("file", name.str());
          file_list.set<int>("number of points", 1000);

          *vo_->os() << "water retention curve file:" << name.str() << std::endl;
        }
      }
    } else if (strcmp(model.c_str(), "brooks_corey")) {
      double lambda = GetAttributeValueD_(element_cp, "lambda");
      double alpha = GetAttributeValueD_(element_cp, "alpha");
      double sr = GetAttributeValueD_(element_cp, "sr");

      for (std::vector<std::string>::const_iterator it = regions.begin(); it != regions.end(); ++it) {
        std::stringstream ss;
        ss << "WRM for " << *it;

        Teuchos::ParameterList& wrm_list = out_list.sublist(ss.str());

        wrm_list.set<std::string>("water retention model", "Brooks Corey")
            .set<std::string>("region", *it)
            .set<double>("Brooks Corey lambda", lambda)
            .set<double>("Brooks Corey alpha", alpha)
            .set<double>("Brooks Corey l", ell)
            .set<double>("residual saturation", sr)
            .set<double>("regularization interval", krel_smooth)
            .set<std::string>("relative permeability model", rel_perm);

        if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH) {
          Teuchos::ParameterList& file_list = wrm_list.sublist("output");
          std::stringstream name;
          name << *it << ".txt";
          file_list.set<std::string>("file", name.str());
          file_list.set<int>("number of points", 1000);

          *vo_->os() << "water retention curve file:" << name.str() << std::endl;
        }
      }
    }
  }

  return out_list;
}


/* ******************************************************************
* Create list of porosity models.
****************************************************************** */
Teuchos::ParameterList InputConverterU::TranslatePOM_() 
{
  Teuchos::ParameterList out_list;

  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH)
      *vo_->os() << "Translating porosity models" << std::endl;

  MemoryManager mm;
  DOMNodeList *node_list, *children;
  DOMNode* node;
  DOMElement* element;

  compressibility_ = false;

  node_list = doc_->getElementsByTagName(mm.transcode("materials"));
  element = static_cast<DOMElement*>(node_list->item(0));
  children = element->getElementsByTagName(mm.transcode("material"));

  for (int i = 0; i < children->getLength(); ++i) {
    DOMNode* inode = children->item(i); 
 
    // get assigned regions
    bool flag;
    node = GetUniqueElementByTagsString_(inode, "assigned_regions", flag);
    std::vector<std::string> regions = CharToStrings_(mm.transcode(node->getTextContent()));

    // get optional complessibility
    node = GetUniqueElementByTagsString_(inode, "mechanical_properties, porosity", flag);
    double phi = GetAttributeValueD_(static_cast<DOMElement*>(node), "value");
    double compres = GetAttributeValueD_(static_cast<DOMElement*>(node), "compressibility", false, 0.0);

    for (std::vector<std::string>::const_iterator it = regions.begin(); it != regions.end(); ++it) {
      std::stringstream ss;
      ss << "POM for " << *it;

      Teuchos::ParameterList& pom_list = out_list.sublist(ss.str());
      pom_list.set<std::string>("region", *it);

      // we can have either uniform of compressible rock
      if (compres == 0.0) {
        pom_list.set<std::string>("porosity model", "constant");
        pom_list.set<double>("value", phi);
      } else {
        pom_list.set<std::string>("porosity model", "compressible");
        pom_list.set<double>("undeformed soil porosity", phi);
        pom_list.set<double>("reference pressure", ATMOSPHERIC_PRESSURE);
        pom_list.set<double>("pore compressibility", compres);
        compressibility_ = true;
      }
    }
  }

  if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH)
      *vo_->os() << "compessibility models: " << compressibility_ << std::endl;

  if (!compressibility_) {
    Teuchos::ParameterList empty;
    out_list = empty;
  }
  return out_list;
}


/* ******************************************************************
* Create list of multiscale models.
****************************************************************** */
Teuchos::ParameterList InputConverterU::TranslateFlowMSM_()
{
  Teuchos::ParameterList out_list;

  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH)
    *vo_->os() << "Translating multiscale models" << std::endl;

  MemoryManager mm;
  DOMNodeList *node_list, *children;
  DOMNode* node;
  DOMElement* element;

  bool flag;
  std::string model, rel_perm;

  node_list = doc_->getElementsByTagName(mm.transcode("materials"));
  element = static_cast<DOMElement*>(node_list->item(0));
  children = element->getElementsByTagName(mm.transcode("material"));

  for (int i = 0; i < children->getLength(); ++i) {
    DOMNode* inode = children->item(i); 
    node = GetUniqueElementByTagsString_(inode, "multiscale_structure, cap_pressure", flag);
    if (!flag) continue;

    model = GetAttributeValueS_(static_cast<DOMElement*>(node), "model", "van_genuchten, brooks_corey");
    DOMNode* nnode = GetUniqueElementByTagsString_(node, "parameters", flag);
    DOMElement* element_cp = static_cast<DOMElement*>(nnode);

    node = GetUniqueElementByTagsString_(inode, "multiscale_structure, rel_perm", flag);
    rel_perm = GetAttributeValueS_(static_cast<DOMElement*>(node), "model", "mualem, burdine");
    DOMNode* mnode = GetUniqueElementByTagsString_(node, "exp", flag);
    DOMElement* element_rp = (flag) ? static_cast<DOMElement*>(mnode) : NULL;

    // common stuff
    // -- assigned regions
    node = GetUniqueElementByTagsString_(inode, "assigned_regions", flag);
    std::vector<std::string> regions = CharToStrings_(mm.transcode(node->getTextContent()));

    // -- mass transfer coefficient
    node = GetUniqueElementByTagsString_(inode, "multiscale_structure, mass_transfer_coefficient", flag);
    double alpha = std::strtod(mm.transcode(node->getTextContent()), NULL);
    
    for (std::vector<std::string>::const_iterator it = regions.begin(); it != regions.end(); ++it) {
      std::stringstream ss;
      ss << "MSM for " << *it;
      Teuchos::ParameterList& msm_list = out_list.sublist(ss.str());

      msm_list.set<std::string>("multiscale model", "dual porosity")
          .set<double>("mass transfer coefficient", alpha);
    }

    // -- porosity models
    node = GetUniqueElementByTagsString_(inode, "multiscale_structure, porosity", flag);
    double phi = GetAttributeValueD_(static_cast<DOMElement*>(node), "value");
    double compres = GetAttributeValueD_(static_cast<DOMElement*>(node), "compressibility", false, 0.0);

    for (std::vector<std::string>::const_iterator it = regions.begin(); it != regions.end(); ++it) {
      std::stringstream ss;
      ss << "MSM for " << *it;
      Teuchos::ParameterList& msm_list = out_list.sublist(ss.str());

      if (compres == 0.0) {
        msm_list.set<std::string>("porosity model", "constant");
        msm_list.set<double>("value", phi);
      } else {
        msm_list.set<std::string>("porosity model", "compressible");
        msm_list.set<double>("undeformed soil porosity", phi);
        msm_list.set<double>("reference pressure", ATMOSPHERIC_PRESSURE);
        msm_list.set<double>("pore compressibility", compres);
      }
    }

    // capillary pressure models
    // -- ell
    double ell, ell_d = (rel_perm == "mualem") ? ELL_MUALEM : ELL_BURDINE;
    ell = GetAttributeValueD_(element_rp, "value", false, ell_d);

    std::replace(rel_perm.begin(), rel_perm.begin() + 1, 'm', 'M');
    std::replace(rel_perm.begin(), rel_perm.begin() + 1, 'b', 'B');

    // -- van Genuchten or Brooks-Corey
    if (strcmp(model.c_str(), "van_genuchten") == 0) {
      double alpha = GetAttributeValueD_(element_cp, "alpha");
      double sr = GetAttributeValueD_(element_cp, "sr");
      double m = GetAttributeValueD_(element_cp, "m");

      for (std::vector<std::string>::const_iterator it = regions.begin(); it != regions.end(); ++it) {
        std::stringstream ss;
        ss << "MSM for " << *it;
        Teuchos::ParameterList& wrm_list = out_list.sublist(ss.str());

        wrm_list.set<std::string>("water retention model", "van Genuchten")
            .set<std::string>("region", *it)
            .set<double>("van Genuchten m", m)
            .set<double>("van Genuchten l", ell)
            .set<double>("van Genuchten alpha", alpha)
            .set<double>("residual saturation", sr)
            .set<std::string>("relative permeability model", rel_perm);
      }
    } else if (strcmp(model.c_str(), "brooks_corey")) {
      double lambda = GetAttributeValueD_(element_cp, "lambda");
      double alpha = GetAttributeValueD_(element_cp, "alpha");
      double sr = GetAttributeValueD_(element_cp, "sr");

      for (std::vector<std::string>::const_iterator it = regions.begin(); it != regions.end(); ++it) {
        std::stringstream ss;
        ss << "MSM for " << *it;

        Teuchos::ParameterList& wrm_list = out_list.sublist(ss.str());

        wrm_list.set<std::string>("water retention model", "Brooks Corey")
            .set<std::string>("region", *it)
            .set<double>("Brooks Corey lambda", lambda)
            .set<double>("Brooks Corey alpha", alpha)
            .set<double>("Brooks Corey l", ell)
            .set<double>("residual saturation", sr)
            .set<std::string>("relative permeability model", rel_perm);
      }
    }
  }

  return out_list;
}


/* ******************************************************************
* Create list of flow BCs.
****************************************************************** */
Teuchos::ParameterList InputConverterU::TranslateFlowBCs_()
{
  Teuchos::ParameterList out_list;

  MemoryManager mm;

  char *text, *tagname;
  DOMNodeList *node_list, *children;
  DOMNode *node;

  node_list = doc_->getElementsByTagName(mm.transcode("boundary_conditions"));
  if (!node_list) return out_list;

  int ibc(0);
  children = node_list->item(0)->getChildNodes();

  for (int i = 0; i < children->getLength(); ++i) {
    DOMNode* inode = children->item(i);
    if (inode->getNodeType() != DOMNode::ELEMENT_NODE) continue;
    tagname = mm.transcode(inode->getNodeName());

    // read the assigned regions
    bool flag;
    node = GetUniqueElementByTagsString_(inode, "assigned_regions", flag);
    text = mm.transcode(node->getTextContent());
    std::vector<std::string> regions = CharToStrings_(text);

    vv_bc_regions_.insert(vv_bc_regions_.end(), regions.begin(), regions.end());

    node = GetUniqueElementByTagsString_(inode, "liquid_phase, liquid_component", flag);
    if (!flag) continue;

    // process a group of similar elements defined by the first element
    std::string bctype;
    std::vector<DOMNode*> same_list = GetSameChildNodes_(node, bctype, flag, true);

    std::map<double, double> tp_values, tp_fluxes;
    std::map<double, std::string> tp_forms;

    for (int j = 0; j < same_list.size(); ++j) {
      DOMNode* jnode = same_list[j];
      double t0 = GetAttributeValueD_(static_cast<DOMElement*>(jnode), "start");

      tp_forms[t0] = GetAttributeValueS_(static_cast<DOMElement*>(jnode), "function");
      tp_values[t0] = GetAttributeValueD_(static_cast<DOMElement*>(jnode), "value", false, 0.0);
      tp_fluxes[t0] = GetAttributeValueD_(static_cast<DOMElement*>(jnode), "inward_mass_flux", false, 0.0);
    }

    // create vectors of values and forms
    std::vector<double> times, values, fluxes;
    std::vector<std::string> forms;
    for (std::map<double, double>::iterator it = tp_values.begin(); it != tp_values.end(); ++it) {
      times.push_back(it->first);
      values.push_back(it->second);
      fluxes.push_back(tp_fluxes[it->second]);
      forms.push_back(tp_forms[it->first]);
    }
    forms.pop_back();

    // create names, modify data
    std::string bcname;
    if (bctype == "inward_mass_flux") {
      bctype = "mass flux";
      bcname = "outward mass flux";
      for (int k = 0; k < values.size(); k++) values[k] *= -1;
    } else if (bctype == "outward_mass_flux") {
      bctype = "mass flux";
      bcname = "outward mass flux";
    } else if (bctype == "outward_volumetric_flux") {
      bctype = "mass flux";
      bcname = "outward mass flux";
      for (int k = 0; k < values.size(); k++) values[k] *= rho_;
    } else if (bctype == "inward_volumetric_flux") {
      bctype = "mass flux";
      bcname = "outward mass flux";
      for (int k = 0; k < values.size(); k++) values[k] *= -rho_;
    } else if (bctype == "uniform_pressure") {
      bctype = "pressure";
      bcname = "boundary pressure";
    } else if (bctype == "hydrostatic") {
      bctype = "static head";
      bcname = "water table elevation";
    } else if (bctype == "seepage_face") {
      bctype = "seepage face";
      bcname = "outward mass flux";
      values = fluxes;
      for (int k = 0; k < values.size(); k++) values[k] *= -1;
    }
    std::stringstream ss;
    ss << "BC " << ibc++;

    // save in the XML files  
    Teuchos::ParameterList& tbc_list = out_list.sublist(bctype);
    Teuchos::ParameterList& bc = tbc_list.sublist(ss.str());
    bc.set<Teuchos::Array<std::string> >("regions", regions);

    Teuchos::ParameterList& bcfn = bc.sublist(bcname);
    if (times.size() == 1) {
      bcfn.sublist("function-constant").set<double>("value", values[0]);
    } else {
      bcfn.sublist("function-tabular")
          .set<Teuchos::Array<double> >("x values", times)
          .set<Teuchos::Array<double> >("y values", values)
          .set<Teuchos::Array<std::string> >("forms", forms);
    }

    // special cases
    if (bctype == "mass flux") {
      bc.set<bool>("rainfall", false);
    } else if (bctype == "static head") {
      std::string tmp = GetAttributeValueS_(
          static_cast<DOMElement*>(same_list[0]), "coordinate_system", false, "absolute");
      bc.set<bool>("relative to top", (tmp == "relative to mesh top"));
    }
  }

  return out_list;
}


/* ******************************************************************
* Create list of flow sources.
****************************************************************** */
Teuchos::ParameterList InputConverterU::TranslateFlowSources_()
{
  Teuchos::ParameterList out_list;

  MemoryManager mm;

  char *text, *tagname;
  DOMNodeList *node_list, *children;
  DOMNode *node, *phase;
  DOMElement* element;

  node_list = doc_->getElementsByTagName(mm.transcode("sources"));
  if (node_list->getLength() == 0) return out_list;

  children = node_list->item(0)->getChildNodes();

  for (int i = 0; i < children->getLength(); ++i) {
    DOMNode* inode = children->item(i);
    if (inode->getNodeType() != DOMNode::ELEMENT_NODE) continue;
    tagname = mm.transcode(inode->getNodeName());
    std::string srcname = GetAttributeValueS_(static_cast<DOMElement*>(inode), "name");

    // read the assigned regions
    bool flag;
    node = GetUniqueElementByTagsString_(inode, "assigned_regions", flag);
    text = mm.transcode(node->getTextContent());
    std::vector<std::string> regions = CharToStrings_(text);

    vv_src_regions_.insert(vv_src_regions_.end(), regions.begin(), regions.end());

    // process flow sources for liquid saturation
    phase = GetUniqueElementByTagsString_(inode, "liquid_phase, liquid_component", flag);
    if (!flag) continue;

    // process a group of similar elements defined by the first element
    std::string srctype;
    std::vector<DOMNode*> same_list = GetSameChildNodes_(phase, srctype, flag, true);
    if (!flag || same_list.size() == 0) continue;

    std::map<double, double> tp_values;
    std::map<double, std::string> tp_forms;
 
    for (int j = 0; j < same_list.size(); ++j) {
       element = static_cast<DOMElement*>(same_list[j]);
       double t0 = GetAttributeValueD_(element, "start");
       tp_forms[t0] = GetAttributeValueS_(element, "function");
       tp_values[t0] = GetAttributeValueD_(element, "value");
    }

    std::string weight;
    if (srctype == "volume_weighted") {
      weight = "volume";
    } else if (srctype == "perm_weighted") {
      weight = "permeability";
    } else if (srctype == "uniform") {
      weight = "none";
    } else {
      ThrowErrorIllformed_("sources", "element", srctype);
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
    Teuchos::ParameterList& src = out_list.sublist(srcname);
    src.set<Teuchos::Array<std::string> >("regions", regions);
    src.set<std::string>("spatial distribution method", weight);

    Teuchos::ParameterList& srcfn = src.sublist("sink");
    if (times.size() == 1) {
      srcfn.sublist("function-constant").set<double>("value", values[0]);
    } else {
      srcfn.sublist("function-tabular")
          .set<Teuchos::Array<double> >("x values", times)
          .set<Teuchos::Array<double> >("y values", values)
          .set<Teuchos::Array<std::string> >("forms", forms);
    }
  }

  return out_list;
}

}  // namespace AmanziInput
}  // namespace Amanzi


