/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Markus Berndt (original version)
           Konstantin Lipnikov (lipnikov@lanl.gov)
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
#include "StringExt.hh"
#include "Units.hh"

#include "InputConverterU.hh"
#include "InputConverterU_Defs.hh"

namespace Amanzi {
namespace AmanziInput {

XERCES_CPP_NAMESPACE_USE

/* ******************************************************************
* Create flow list.
****************************************************************** */
Teuchos::ParameterList
InputConverterU::TranslateFlow_(const std::string& mode,
                                const std::string& domain,
                                const std::string& pk_model)
{
  Teuchos::ParameterList out_list;
  Teuchos::ParameterList* flow_list = nullptr;

  if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH)
    *vo_->os() << "Translating flow, mode=" << mode << " domain=" << domain << std::endl;

  MemoryManager mm;
  DOMNode* node;

  std::string msm;

  // set up default values for some expert parameters
  std::string rel_perm("upwind-darcy_velocity"), rel_perm_out;
  std::string update_upwind("every timestep");

  // process expert parameters
  bool flag;
  node = GetUniqueElementByTagsString_(
    "unstructured_controls, unstr_flow_controls, rel_perm_method", flag);
  if (flag) rel_perm = mm.transcode(node->getTextContent());

  rel_perm_out = rel_perm;
  Amanzi::replace_all(rel_perm_out, "-", ": ");
  std::replace(rel_perm_out.begin(), rel_perm_out.end(), '_', ' ');

  node = GetUniqueElementByTagsString_(
    "unstructured_controls, unstr_flow_controls, update_upwind_frequency", flag);
  if (flag) update_upwind = mm.transcode(node->getTextContent());
  replace(update_upwind.begin(), update_upwind.end(), '_', ' ');

  node = GetUniqueElementByTagsString_(
    "unstructured_controls, unstr_flow_controls, optional_fields", flag);
  if (flag) {
    auto fields = CharToStrings_(mm.transcode(node->getTextContent()));
    out_list.set<Teuchos::Array<std::string>>("optional fields", fields);
  }

  // create flow header
  out_list.set<std::string>("domain name", (domain == "matrix") ? "domain" : domain);

  out_list.sublist("fracture permeability models") = TranslateFlowFractures_(domain);
  if (out_list.sublist("fracture permeability models").numParams() > 0) {
    out_list.sublist("physical models and assumptions")
      .set<bool>("flow and transport in fractures", true);
  }

  if (pk_model == "darcy") {
    flow_list = &out_list;
    flow_single_phase_ = true;

    DOMNodeList* tmp = doc_->getElementsByTagName(mm.transcode("multiscale_model"));
    if (tmp->getLength() > 0) {
      msm = GetAttributeValueS_(tmp->item(0), "name");
      flow_list->sublist("physical models and assumptions")
        .set<std::string>("multiscale model", "dual continuum discontinuous matrix");
    }

  } else if (pk_model == "richards") {
    Teuchos::ParameterList& upw_list = out_list.sublist("relative permeability");
    upw_list.set<std::string>("upwind method", rel_perm_out);
    upw_list.set<std::string>("upwind frequency", update_upwind);
    upw_list.sublist("upwind parameters")
      .set<double>("tolerance", 1e-12)
      .set<std::string>("method", "cell-based")
      .set<int>("polynomial order", 1)
      .set<std::string>("limiter", "Barth-Jespersen");
    flow_list = &out_list;

    out_list.sublist("water retention models") = TranslateWRM_("flow");
    out_list.sublist("porosity models") = TranslatePOM_(domain);
    if (out_list.sublist("porosity models").numParams() > 0) {
      flow_list->sublist("physical models and assumptions")
        .set<std::string>("porosity model", "compressible: pressure function");
    }
    out_list.sublist("multiscale models") = TranslateFlowMSM_();
    if (out_list.sublist("multiscale models").numParams() > 0) {
      msm = "dual continuum discontinuous matrix";
      flow_list->sublist("physical models and assumptions")
        .set<std::string>("multiscale model", msm);
    }
    out_list.sublist("permeability porosity models") = TranslatePPM_(domain);
    if (out_list.sublist("permeability porosity models").numParams() > 0) {
      flow_list->sublist("physical models and assumptions")
        .set<bool>("permeability porosity model", true);
    }

  } else {
    Errors::Message msg;
    msg << "Internal error for flow model \"" << pk_model << "\".\n";
    Exceptions::amanzi_throw(msg);
  }

  // absolute permeability (default list so far)
  flow_list->sublist("absolute permeability").set<std::string>("coordinate system", "cartesian");

  // insert operator sublist
  std::string disc_method("mfd-optimized_for_sparsity");
  node = GetUniqueElementByTagsString_(
    "unstructured_controls, unstr_flow_controls, discretization_method", flag);
  if (flag) disc_method = mm.transcode(node->getTextContent());

  std::string pc_method("linearized_operator");
  if (pk_model == "darcy") pc_method = "diffusion_operator";
  node = GetUniqueElementByTagsString_(
    "unstructured_controls, unstr_flow_controls, preconditioning_strategy", flag);
  if (flag) pc_method = GetTextContentS_(node, "linearized_operator, diffusion_operator");

  std::string nonlinear_solver("nka");
  node = GetUniqueElementByTagsString_("unstructured_controls, unstr_nonlinear_solver", flag);
  if (flag) nonlinear_solver = GetAttributeValueS_(node, "name", TYPE_NONE, false, "nka");

  bool modify_correction(false);
  node = GetUniqueElementByTagsString_(
    "unstructured_controls, unstr_nonlinear_solver, modify_correction", flag);

  // Newton method requires to overwrite some parameters.
  if (nonlinear_solver == "newton") {
    modify_correction = true;
    out_list.sublist("relative permeability")
      .set<std::string>("upwind frequency", "every nonlinear iteration");

    if (disc_method != "fv-default" || rel_perm != "upwind-darcy_velocity" || !modify_correction) {
      Errors::Message msg;
      msg << "Nonlinear solver \"newton\" requires \"upwind-darcy_velocity\" (is \"" << rel_perm
          << "\")\n";
      msg << "\"modify_correction\"=true (is " << modify_correction
          << "), discretization method \"fv-default\" (is \"" << disc_method << "\").\n";
      Exceptions::amanzi_throw(msg);
    }
  }

  // Newton-Picard method requires to overwrite some parameters.
  if (nonlinear_solver == "newton-picard") {
    out_list.sublist("relative permeability")
      .set<std::string>("upwind frequency", "every nonlinear iteration");
  }

  std::string nonlinear_coef;
  if (pk_model == "richards") nonlinear_coef = rel_perm;
  flow_list->sublist("operators") = TranslateDiffusionOperator_(
    disc_method, pc_method, nonlinear_solver, nonlinear_coef, "vapor matrix", domain, true);

  // insert time integrator
  std::string err_options, unstr_controls;
  if (mode == "steady") {
    err_options = "pressure";
    unstr_controls = "numerical_controls, unstructured_controls, unstr_steady-state_controls";
  } else {
    err_options = "pressure, residual";
    unstr_controls = "numerical_controls, unstructured_controls, unstr_transient_controls";

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
    flow_list->sublist("time integrator") = TranslateTimeIntegrator_(err_options,
                                                                     nonlinear_solver,
                                                                     modify_correction,
                                                                     unstr_controls,
                                                                     TI_SOLVER,
                                                                     dt_cut_[mode],
                                                                     dt_inc_[mode]);
  }

  // insert boundary conditions and source terms
  flow_list->sublist("boundary conditions") = TranslateFlowBCs_(domain);
  flow_list->sublist("source terms") = TranslateSources_(domain, "flow");

  flow_list->sublist("verbose object") = verb_list_.sublist("verbose object");

  // statistics
  Teuchos::OSTab tab = vo_->getOSTab();
  if (msm.size() > 0 && vo_->getVerbLevel() >= Teuchos::VERB_HIGH)
    *vo_->os() << "found at least one multiscale model: \"" << msm << "\"\n";

  return out_list;
}


/* ******************************************************************
* Create list of water retention models.
****************************************************************** */
Teuchos::ParameterList
InputConverterU::TranslateWRM_(const std::string& pk_name)
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
    model =
      GetAttributeValueS_(node, "model", "van_genuchten, brooks_corey, saturated, linear, tabular");
    DOMNode* nnode = GetUniqueElementByTagsString_(node, "parameters", flag);
    DOMElement* element_cp = static_cast<DOMElement*>(nnode);

    node = GetUniqueElementByTagsString_(inode, "rel_perm", flag);
    DOMElement* element_rp = NULL;
    if (flag) {
      rel_perm = GetAttributeValueS_(node, "model", "mualem, burdine");
      DOMNode* mnode = GetUniqueElementByTagsString_(node, "exp", flag);
      element_rp = static_cast<DOMElement*>(mnode);
    }

    // common stuff
    // -- assigned regions
    node = GetUniqueElementByTagsString_(inode, "assigned_regions", flag);
    std::vector<std::string> regions = CharToStrings_(mm.transcode(node->getTextContent()));

    // -- smoothing
    double krel_smooth = GetAttributeValueD_(element_cp,
                                             "optional_krel_smoothing_interval",
                                             TYPE_NUMERICAL,
                                             0.0,
                                             DVAL_MAX,
                                             "",
                                             false,
                                             0.0);

    // -- ell
    double ell(0.0);
    if (element_rp != NULL) {
      double ell_d = (rel_perm == "mualem") ? ELL_MUALEM : ELL_BURDINE;
      ell = GetAttributeValueD_(element_rp, "value", TYPE_NUMERICAL, 0.0, 10.0, "", false, ell_d);
    }

    std::replace(rel_perm.begin(), rel_perm.begin() + 1, 'm', 'M');
    std::replace(rel_perm.begin(), rel_perm.begin() + 1, 'b', 'B');

    std::stringstream ss;
    ss << "WRM_" << i;
    Teuchos::ParameterList& wrm_list = out_list.sublist(ss.str());

    if (strcmp(model.c_str(), "van_genuchten") == 0) {
      double alpha =
        GetAttributeValueD_(element_cp, "alpha", TYPE_NUMERICAL, 0.0, DVAL_MAX, "Pa^-1");
      double sr = GetAttributeValueD_(element_cp, "sr", TYPE_NUMERICAL, 0.0, 1.0, "-");
      double m = GetAttributeValueD_(element_cp, "m", TYPE_NUMERICAL, 0.0, DVAL_MAX, "-");

      wrm_list.set<std::string>("water retention model", "van Genuchten")
        .set<Teuchos::Array<std::string>>("regions", regions)
        .set<double>("van Genuchten m", m)
        .set<double>("van Genuchten n", 1.0 / (1 - m))
        .set<double>("van Genuchten l", ell)
        .set<double>("van Genuchten alpha", alpha)
        .set<double>("residual saturation liquid", sr)
        .set<std::string>("relative permeability model", rel_perm);
      if (pk_name == "flow") {
        wrm_list.set<double>("regularization interval", krel_smooth);
      } else if (pk_name == "multiphase") {
        wrm_list.set<double>("regularization interval kr", krel_smooth)
          .set<double>("regularization interval pc", krel_smooth);
      }

      if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH) {
        Teuchos::ParameterList& file_list = wrm_list.sublist("output");
        std::string name = ss.str() + ".txt";
        file_list.set<std::string>("file", name);
        file_list.set<int>("number of points", 1000);

        *vo_->os() << "water retention curve file: wrm_" << name << std::endl;
      }
    } else if (strcmp(model.c_str(), "brooks_corey") == 0) {
      double lambda = GetAttributeValueD_(element_cp, "lambda", TYPE_NUMERICAL, 0.0, DVAL_MAX);
      double alpha =
        GetAttributeValueD_(element_cp, "alpha", TYPE_NUMERICAL, DVAL_MIN, DVAL_MAX, "Pa^-1");
      double sr = GetAttributeValueD_(element_cp, "sr", TYPE_NUMERICAL, 0.0, 1.0);

      wrm_list.set<std::string>("water retention model", "Brooks Corey")
        .set<Teuchos::Array<std::string>>("regions", regions)
        .set<double>("Brooks Corey lambda", lambda)
        .set<double>("Brooks Corey alpha", alpha)
        .set<double>("Brooks Corey l", ell)
        .set<double>("residual saturation liquid", sr)
        .set<double>("regularization interval", krel_smooth)
        .set<std::string>("relative permeability model", rel_perm);

      if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH) {
        Teuchos::ParameterList& file_list = wrm_list.sublist("output");
        std::string name = ss.str() + ".txt";
        file_list.set<std::string>("file", name);
        file_list.set<int>("number of points", 1000);

        *vo_->os() << "water retention curve file:" << name << std::endl;
      }
    } else if (strcmp(model.c_str(), "linear") == 0) {
      double alpha =
        GetAttributeValueD_(element_cp, "alpha", TYPE_NUMERICAL, DVAL_MIN, DVAL_MAX, "Pa^-1");
      double sr = GetAttributeValueD_(element_cp, "sr", TYPE_NUMERICAL, 0.0, 1.0);

      wrm_list.set<std::string>("water retention model", "linear")
        .set<Teuchos::Array<std::string>>("regions", regions)
        .set<double>("pc0", 1.0 / alpha)
        .set<double>("residual saturation liquid", sr);

    } else if (strcmp(model.c_str(), "saturated") == 0) {
      wrm_list.set<std::string>("water retention model", "saturated")
        .set<Teuchos::Array<std::string>>("regions", regions);

    } else if (strcmp(model.c_str(), "tabular") == 0) {
      auto cp = GetAttributeVectorD_(element_cp, "cp");
      auto kr = GetAttributeVectorD_(element_cp, "kr");
      auto sat = GetAttributeVectorD_(element_cp, "s");

      wrm_list.set<std::string>("water retention model", "tabular")
        .set<Teuchos::Array<std::string>>("regions", regions)
        .set<Teuchos::Array<double>>("cap pressure", cp)
        .set<Teuchos::Array<double>>("permeability", kr)
        .set<Teuchos::Array<double>>("saturation", sat);
    }
  }

  return out_list;
}


/* ******************************************************************
* Create list of porosity models.
****************************************************************** */
Teuchos::ParameterList
InputConverterU::TranslatePOM_(const std::string& domain)
{
  Teuchos::ParameterList out_list;

  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH)
    *vo_->os() << "Translating porosity models" << std::endl;

  MemoryManager mm;
  DOMNodeList* children;
  DOMNode* node;
  DOMElement* element;

  bool flag;

  compressibility_ = false;

  node = (domain == "fracture") ?
           GetUniqueElementByTagsString_("fracture_network, materials", flag) :
           GetUniqueElementByTagsString_("materials", flag);
  element = static_cast<DOMElement*>(node);
  children = element->getElementsByTagName(mm.transcode("material"));

  for (int i = 0; i < children->getLength(); ++i) {
    DOMNode* inode = children->item(i);

    // get assigned regions
    bool flag;
    node = GetUniqueElementByTagsString_(inode, "assigned_regions", flag);
    std::vector<std::string> regions = CharToStrings_(mm.transcode(node->getTextContent()));

    // get optional complessibility
    node = GetUniqueElementByTagsString_(inode, "mechanical_properties, porosity", flag);
    std::string type = GetAttributeValueS_(node, "type", TYPE_NONE, false, "");
    if (type == "h5file") {
      compressibility_ = false;
      break;
    }

    double phi = GetAttributeValueD_(node, "value", TYPE_NUMERICAL, 0.0, 1.0);
    double compres =
      GetAttributeValueD_(node, "compressibility", TYPE_NUMERICAL, 0.0, 1.0, "Pa^-1", false, 0.0);

    std::stringstream ss;
    ss << "POM " << i;

    Teuchos::ParameterList& pom_list = out_list.sublist(ss.str());
    pom_list.set<Teuchos::Array<std::string>>("regions", regions);

    // we can have either uniform of compressible rock
    if (compres == 0.0) {
      pom_list.set<std::string>("porosity model", "constant");
      pom_list.set<double>("value", phi);
    } else {
      pom_list.set<std::string>("porosity model", "compressible");
      pom_list.set<double>("undeformed soil porosity", phi);
      pom_list.set<double>("reference pressure", const_atm_pressure_);
      pom_list.set<double>("pore compressibility", compres);
      compressibility_ = true;
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
Teuchos::ParameterList
InputConverterU::TranslateFlowMSM_()
{
  Teuchos::ParameterList out_list;

  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH)
    *vo_->os() << "Translating multiscale models" << std::endl;

  MemoryManager mm;
  DOMNodeList *node_list, *children;
  DOMNode *node, *knode;
  DOMElement* element;

  bool flag;
  std::string name, dual, model, rel_perm;

  node_list = doc_->getElementsByTagName(mm.transcode("materials_secondary_continuum"));
  if (node_list->getLength() == 0) return out_list;

  element = static_cast<DOMElement*>(node_list->item(0));
  children = element->getElementsByTagName(mm.transcode("material"));

  for (int i = 0; i < children->getLength(); ++i) {
    DOMNode* inode = children->item(i);
    knode = GetUniqueElementByTagsString_(inode, "multiscale_model", flag);
    if (!flag) ThrowErrorMissing_("materials", "element", "multiscale_model", "material");
    dual = GetAttributeValueS_(knode, "name");

    // verify material name
    name = GetAttributeValueS_(inode, "name");
    if (!FindNameInVector_(name, material_names_)) {
      Errors::Message msg("Materials names in primary and secondary continua do not match.\n");
      Exceptions::amanzi_throw(msg);
    }

    // common stuff
    std::stringstream ss;
    ss << "MSM " << i;
    Teuchos::ParameterList& msm_list = out_list.sublist(ss.str());
    Teuchos::ParameterList& msm_slist = msm_list.sublist(dual + " parameters");

    // -- assigned regions
    node = GetUniqueElementByTagsString_(inode, "assigned_regions", flag);
    std::vector<std::string> regions = CharToStrings_(mm.transcode(node->getTextContent()));

    msm_list.set<Teuchos::Array<std::string>>("regions", regions)
      .set<std::string>("multiscale model", dual);

    // -- volume fraction
    node = GetUniqueElementByTagsString_(inode, "volume_fraction", flag);
    double vof = GetTextContentD_(node, "", true);

    // water retention model
    node = GetUniqueElementByTagsString_(inode, "cap_pressure", flag);

    model = GetAttributeValueS_(node, "model", "van_genuchten, brooks_corey");
    DOMNode* nnode = GetUniqueElementByTagsString_(node, "parameters", flag);
    DOMElement* element_cp = static_cast<DOMElement*>(nnode);

    node = GetUniqueElementByTagsString_(inode, "rel_perm", flag);
    rel_perm = GetAttributeValueS_(node, "model", "mualem, burdine");
    DOMNode* mnode = GetUniqueElementByTagsString_(node, "exp", flag);
    DOMElement* element_rp = (flag) ? static_cast<DOMElement*>(mnode) : NULL;

    // dual continuum models
    if (dual == "dual porosity") {
      node = GetUniqueElementByTagsString_(knode, "mass_transfer_coefficient", flag);
      double alpha = GetTextContentD_(node, "s^-1", true);

      msm_slist.set<double>("mass transfer coefficient", alpha);
    } else if (dual == "generalized dual porosity") {
      node = GetUniqueElementByTagsString_(knode, "matrix", flag);
      if (!flag) ThrowErrorMissing_("materials", "element", "matrix", "multiscale_model");

      int nnodes =
        GetAttributeValueL_(node, "number_of_nodes", TYPE_NUMERICAL, 0, INT_MAX, false, 1);
      double depth = GetAttributeValueD_(node, "depth", TYPE_NUMERICAL, 0.0, DVAL_MAX, "m");
      double perm = GetAttributeValueD_(node, "permeability", TYPE_NUMERICAL, 0.0, DVAL_MAX, "m^2");

      msm_slist.set<int>("number of matrix nodes", nnodes)
        .set<double>("matrix depth", depth)
        .set<double>("matrix volume fraction", vof)
        .set<double>("absolute permeability", perm);
    }

    // porosity models
    node = GetUniqueElementByTagsString_(inode, "mechanical_properties, porosity", flag);
    double phi = GetAttributeValueD_(node, "value", TYPE_NUMERICAL, 0.0, 1.0);
    double compres =
      GetAttributeValueD_(node, "compressibility", TYPE_NUMERICAL, 0.0, 1.0, "Pa^-1", false, 0.0);

    if (compres == 0.0) {
      msm_list.set<std::string>("porosity model", "constant");
      msm_list.set<double>("value", phi);
    } else {
      msm_list.set<std::string>("porosity model", "compressible");
      msm_list.set<double>("undeformed soil porosity", phi);
      msm_list.set<double>("reference pressure", const_atm_pressure_);
      msm_list.set<double>("pore compressibility", compres);
    }

    // capillary pressure models
    // -- ell
    double ell, ell_d = (rel_perm == "mualem") ? ELL_MUALEM : ELL_BURDINE;
    ell = GetTextContentD_(element_rp, "", false, ell_d);

    std::replace(rel_perm.begin(), rel_perm.begin() + 1, 'm', 'M');
    std::replace(rel_perm.begin(), rel_perm.begin() + 1, 'b', 'B');

    // -- van Genuchten or Brooks-Corey
    if (strcmp(model.c_str(), "van_genuchten") == 0) {
      double alpha =
        GetAttributeValueD_(element_cp, "alpha", TYPE_NUMERICAL, 0.0, DVAL_MAX, "Pa^-1");
      double sr = GetAttributeValueD_(element_cp, "sr", TYPE_NUMERICAL, 0.0, 1.0);
      double m = GetAttributeValueD_(element_cp, "m", TYPE_NUMERICAL, 0.0, DVAL_MAX);

      msm_list.set<std::string>("water retention model", "van Genuchten")
        .set<double>("van Genuchten m", m)
        .set<double>("van Genuchten l", ell)
        .set<double>("van Genuchten alpha", alpha)
        .set<double>("residual saturation liquid", sr)
        .set<std::string>("relative permeability model", rel_perm);
    } else if (strcmp(model.c_str(), "brooks_corey") == 0) {
      double lambda = GetAttributeValueD_(element_cp, "lambda", TYPE_NUMERICAL, 0.0, DVAL_MAX);
      double alpha =
        GetAttributeValueD_(element_cp, "alpha", TYPE_NUMERICAL, 0.0, DVAL_MAX, "Pa^-1");
      double sr = GetAttributeValueD_(element_cp, "sr", TYPE_NUMERICAL, 0.0, 1.0);

      msm_list.set<std::string>("water retention model", "Brooks Corey")
        .set<double>("Brooks Corey lambda", lambda)
        .set<double>("Brooks Corey alpha", alpha)
        .set<double>("Brooks Corey l", ell)
        .set<double>("residual saturation liquid", sr)
        .set<std::string>("relative permeability model", rel_perm);
    }
  }

  return out_list;
}


/* ******************************************************************
* Create list of permeability porosity models.
****************************************************************** */
Teuchos::ParameterList
InputConverterU::TranslatePPM_(const std::string& domain)
{
  Teuchos::ParameterList out_list;

  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH)
    *vo_->os() << "Translating permeability-porosity models" << std::endl;

  MemoryManager mm;
  DOMNodeList* children;
  DOMNode* node;
  DOMElement* element;

  bool flag, found(false);

  node = (domain == "fracture") ?
           GetUniqueElementByTagsString_("fracture_network, materials", flag) :
           GetUniqueElementByTagsString_("materials", flag);
  element = static_cast<DOMElement*>(node);
  children = element->getElementsByTagName(mm.transcode("material"));

  for (int i = 0; i < children->getLength(); ++i) {
    DOMNode* inode = children->item(i);

    // get assigned regions
    bool flag;
    node = GetUniqueElementByTagsString_(inode, "assigned_regions", flag);
    std::vector<std::string> regions = CharToStrings_(mm.transcode(node->getTextContent()));

    // get optional compressibility
    node = GetUniqueElementByTagsString_(inode, "permeability", flag);
    std::string model = GetAttributeValueS_(node, "model", TYPE_NONE, false, "");

    std::stringstream ss;
    ss << "PPM " << i;

    Teuchos::ParameterList& ppm_list = out_list.sublist(ss.str());
    ppm_list.set<Teuchos::Array<std::string>>("regions", regions);

    if (model == "power law") {
      found = true;
      double exp = GetAttributeValueD_(node, "exponent", TYPE_NUMERICAL, 0.0, 10.0, "-");
      double phi = GetAttributeValueD_(node, "reference_porosity", TYPE_NUMERICAL, 0.0, 1.0, "-");

      ppm_list.set<std::string>("permeability porosity model", model)
        .set<double>("undeformed soil porosity", phi)
        .set<double>("power law exponent", exp);
    } else if (model == "Kozeny-Carman") {
      found = true;
      double phi = GetAttributeValueD_(node, "reference_porosity", TYPE_NUMERICAL, 0.0, 1.0, "-");

      ppm_list.set<std::string>("permeability porosity model", model)
        .set<double>("undeformed soil porosity", phi);
    }
  }

  if (!found) {
    Teuchos::ParameterList empty;
    out_list = empty;
  }
  return out_list;
}


/* ******************************************************************
* Create list of flow BCs.
****************************************************************** */
Teuchos::ParameterList
InputConverterU::TranslateFlowBCs_(const std::string& domain)
{
  Teuchos::ParameterList out_list;

  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH)
    *vo_->os() << "Translating boundary conditions" << std::endl;

  MemoryManager mm;

  char* text;
  DOMNodeList* node_list;
  DOMNode* node;
  DOMElement* element;

  // correct list of boundary conditions for given domain
  bool flag;
  if (domain == "matrix")
    node = GetUniqueElementByTagsString_("boundary_conditions", flag);
  else
    node = GetUniqueElementByTagsString_("fracture_network, boundary_conditions", flag);
  if (!flag) return out_list;

  node_list = node->getChildNodes();

  // parse the list of BCs using non-intrusive first_child -> next_sibling loop
  std::set<std::string> active_bcs; // for statistics
  int ibc(0);

  DOMNode* inode = node_list->item(0);
  while (inode != NULL) {
    if (inode->getNodeType() != DOMNode::ELEMENT_NODE) {
      inode = inode->getNextSibling();
      continue;
    }

    // read the assigned regions
    node = GetUniqueElementByTagsString_(inode, "assigned_regions", flag);
    text = mm.transcode(node->getTextContent());
    std::vector<std::string> regions = CharToStrings_(text);

    node = GetUniqueElementByTagsString_(inode, "liquid_phase, liquid_component", flag);
    if (!flag) {
      inode = inode->getNextSibling();
      continue;
    }

    // process a group of similar elements defined by the first element
    // -- get BC type
    std::string bctype_in;
    std::vector<DOMNode*> same_list = GetSameChildNodes_(node, bctype_in, flag, true);

    // -- exceptions
    if (bctype_in == "no_flow") {
      Errors::Message msg;
      msg << "\"no-flow\" is the default BC. Use \"inward_mass_flux\" to bypass this error.\n";
      Exceptions::amanzi_throw(msg);
    }

    // -- identify a BC that do not require forms (the global BC)
    bool global_bc(false);
    if (bctype_in == "linear_pressure" || bctype_in == "linear_hydrostatic") { global_bc = true; }

    // -- identify a hard-coded BC that uses spatially dependent functions
    //    temporarily, we assume that it is also the global BC.
    bool space_bc, time_bc;
    DOMNode* knode = GetUniqueElementByTagsString_(same_list[0], "space", space_bc);
    DOMNode* lnode = GetUniqueElementByTagsString_(same_list[0], "time", time_bc);
    if (space_bc) global_bc = true;

    // -- define the expected unit
    std::string unit("kg/s/m^2");
    if (bctype_in == "outward_volumetric_flux" || bctype_in == "inward_volumetric_flux") {
      unit = "m/s";
    } else if (bctype_in == "uniform_pressure" || bctype_in == "linear_pressure") {
      unit = "Pa";
    } else if (bctype_in == "hydrostatic" || bctype_in == "linear_hydrostatic") {
      unit = "m";
    } else if (bctype_in == "inward_mass_flux" || bctype_in == "seepage_face") {
      unit = "kg/s/m^2";
    } else if (bctype_in == "inward_mass_flux_distributed") {
      unit = "kg/s";
    } else { // not flow BCs
      inode = inode->getNextSibling();
      continue;
    }

    // -- process global and local BC separately
    double refv;
    std::vector<double> grad, refc, data, data_tmp;
    BCs bcs;

    if (space_bc) {
      data.push_back(
        GetAttributeValueD_(knode, "amplitude", TYPE_NUMERICAL, DVAL_MIN, DVAL_MAX, "kg/m^2/s"));
      data_tmp = GetAttributeVectorD_(knode, "center", dim_, "m");
      data.insert(data.end(), data_tmp.begin(), data_tmp.end());
      data.push_back(
        GetAttributeValueD_(knode, "standard_deviation", TYPE_NUMERICAL, 0.0, DVAL_MAX, "m"));

      if (time_bc) {
        data[0] *= GetAttributeValueD_(lnode, "data", TYPE_NUMERICAL, DVAL_MIN, DVAL_MAX, "");
      }
    } else if (global_bc) {
      std::string unit_grad = unit + "/m";
      element = static_cast<DOMElement*>(same_list[0]);
      refv = GetAttributeValueD_(element, "value", TYPE_NUMERICAL, DVAL_MIN, DVAL_MAX, unit);
      grad = GetAttributeVectorD_(element, "gradient", dim_, unit_grad);
      refc = GetAttributeVectorD_(element, "reference_coord", dim_, "m");
    } else {
      bcs = ParseCondList_(node, DVAL_MIN, DVAL_MAX, unit);
    }

    // -- create BC names, modify input data
    std::string bcname, bctype(bctype_in);
    if (bctype_in == "inward_mass_flux") {
      bctype = "mass flux";
      bcname = "outward mass flux";
      for (int k = 0; k < bcs.values.size(); k++) bcs.values[k] *= -1;
    } else if (bctype_in == "inward_mass_flux_distributed") {
      bctype = "mass flux";
      bcname = "outward mass flux";
      for (int k = 0; k < bcs.values.size(); k++) bcs.values[k] *= -1;
      for (int k = 0; k < bcs.forms.size(); k++) bcs.forms[k] = "volume";
    } else if (bctype_in == "outward_mass_flux") {
      bctype = "mass flux";
      bcname = "outward mass flux";
    } else if (bctype_in == "outward_volumetric_flux") {
      bctype = "mass flux";
      bcname = "outward mass flux";
      for (int k = 0; k < bcs.values.size(); k++) bcs.values[k] *= rho_;
    } else if (bctype_in == "inward_volumetric_flux") {
      bctype = "mass flux";
      bcname = "outward mass flux";
      for (int k = 0; k < bcs.values.size(); k++) bcs.values[k] *= -rho_;
    } else if (bctype_in == "uniform_pressure" || bctype_in == "linear_pressure") {
      bctype = "pressure";
      bcname = "boundary pressure";
    } else if (bctype_in == "hydrostatic" || bctype_in == "linear_hydrostatic") {
      bctype = "static head";
      bcname = "water table elevation";
    } else if (bctype_in == "seepage_face") {
      bctype = "seepage face";
      bcname = "outward mass flux";
      bcs.values = bcs.fluxes;
      for (int k = 0; k < bcs.values.size(); k++) bcs.values[k] *= -1;
    } else {
      ThrowErrorIllformed_("boundary_conditions", "element", bctype_in);
    }
    active_bcs.insert(bctype);

    // save in the XML files
    std::stringstream ss;
    ss << "BC " << ibc++;

    Teuchos::ParameterList& tbc_list = out_list.sublist(bctype);
    Teuchos::ParameterList& bc = tbc_list.sublist(ss.str());
    bc.set<Teuchos::Array<std::string>>("regions", regions);

    // select one region for transport diagnostics (FIXME)
    if (bctype == "seepage face")
      transport_diagnostics_.insert(transport_diagnostics_.end(), regions.begin(), regions.end());

    Teuchos::ParameterList& bcfn = bc.sublist(bcname);
    if (space_bc) { // only one use case so far
      TranslateFunctionGaussian_(data, bcfn);
    } else if (global_bc) {
      TranslateFunctionGradient_(refv, grad, refc, bcfn);
    } else if (bcs.filename != "") {
      bcfn.sublist("function-tabular")
        .set<std::string>("file", bcs.filename)
        .set<std::string>("x header", bcs.xheader)
        .set<std::string>("y header", bcs.yheader);
    } else {
      TranslateGenericMath_(bcs, bcfn);
    }

    // data distribution method (only one distribution method at the moment - by volume)
    bc.set<bool>("use area fractions", WeightVolumeSubmodel_(regions));
    bc.set<std::string>("spatial distribution method", "none");

    if (bcs.forms.size() == 1 && bcs.forms[0] == "volume") {
      if (domain == "domain") {
        bc.set<std::string>("spatial distribution method", "volume");
      } else {
        bc.set<std::string>("spatial distribution method", "weight by field")
          .set<std::string>("field key", "fracture-aperture");
      }
    }

    // special cases and parameters without default values
    element = static_cast<DOMElement*>(same_list[0]);
    if (bctype == "mass flux") {
      bc.set<bool>("rainfall", false);
    } else if (bctype == "seepage face") {
      double tmp =
        GetAttributeValueD_(element, "flux_threshold", TYPE_NUMERICAL, 0.0, 0.1, "", false, 0.0);
      bc.set<bool>("rainfall", false)
        .set<std::string>("submodel", "PFloTran")
        .set<double>("seepage flux threshold", tmp);
    } else if (bctype == "static head") {
      std::string tmp =
        GetAttributeValueS_(element, "coordinate_system", TYPE_NONE, false, "absolute");
      bc.set<bool>("relative to top", (tmp == "relative to mesh top"));
      bc.set<bool>("relative to bottom", (tmp == "relative to mesh bottom"));

      Teuchos::ParameterList& bc_tmp = bc.sublist("static head").sublist("function-static-head");
      bc_tmp.set<int>("space dimension", dim_)
        .set<double>("density", rho_)
        .set<double>("gravity", (gravity_on_) ? const_gravity_ : 0.0)
        .set<double>("p0", const_atm_pressure_);
      bc_tmp.sublist(bcname) = bcfn;
      bc.remove(bcname);

      tmp = GetAttributeValueS_(element, "submodel", TYPE_NONE, false, "none");
      bc.set<bool>("no flow above water table", (tmp == "no_flow_above_water_table"));
    }

    vv_bc_regions_.insert(vv_bc_regions_.end(), regions.begin(), regions.end());

    inode = inode->getNextSibling();
  }

  // output statistics
  if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH) {
    Teuchos::OSTab tab2 = vo_->getOSTab();
    *vo_->os() << "active BCs: ";
    for (std::set<std::string>::iterator it = active_bcs.begin(); it != active_bcs.end(); ++it) {
      *vo_->os() << *it << ", ";
    }
    *vo_->os() << "no flow (default)" << std::endl;
  }

  return out_list;
}


/* ******************************************************************
* Create list of flow sources.
****************************************************************** */
Teuchos::ParameterList
InputConverterU::TranslateSources_(const std::string& domain, const std::string& pkname)
{
  Teuchos::ParameterList out_list;

  MemoryManager mm;

  bool flag;
  char* text;
  DOMNodeList* children;
  DOMNode *node, *phase;
  DOMElement* element;

  if (domain == "fracture")
    node = GetUniqueElementByTagsString_("fracture_network, sources", flag);
  else
    node = GetUniqueElementByTagsString_("sources", flag);

  if (!flag) return out_list;
  children = node->getChildNodes();

  for (int i = 0; i < children->getLength(); ++i) {
    DOMNode* inode = children->item(i);
    if (inode->getNodeType() != DOMNode::ELEMENT_NODE) continue;
    std::string srcname = GetAttributeValueS_(static_cast<DOMElement*>(inode), "name");

    // read the assigned regions
    node = GetUniqueElementByTagsString_(inode, "assigned_regions", flag);
    text = mm.transcode(node->getTextContent());
    std::vector<std::string> regions = CharToStrings_(text);

    vv_src_regions_.insert(vv_src_regions_.end(), regions.begin(), regions.end());

    // process flow sources for liquid saturation
    if (pkname == "flow") {
      phase = GetUniqueElementByTagsString_(inode, "liquid_phase, liquid_component", flag);
    } else if (pkname == "energy") {
      phase = GetUniqueElementByTagsString_(inode, "thermal_component", flag);
    } else {
      flag = false;
    }
    if (!flag) continue;

    // process a group of similar elements defined by the first element
    std::string srctype, weight, unit;
    std::vector<DOMNode*> same_list = GetSameChildNodes_(phase, srctype, flag, true);
    if (!flag || same_list.size() == 0) continue;

    if (srctype == "volume_weighted") {
      weight = "volume";
      if (pkname == "flow")
        unit = "kg/s";
      else if (pkname == "energy")
        unit = "J/s";
    } else if (srctype == "perm_weighted") {
      weight = "permeability";
      unit = "kg/s";
    } else if (srctype == "uniform") {
      weight = "none";
      unit = "kg/m^3/s";
    } else if (srctype == "peaceman_well") {
      weight = "simple well";
      unit = "Pa";
    } else {
      ThrowErrorIllformed_("sources", "element", srctype);
    }

    auto bcs = ParseCondList_(phase, DVAL_MIN, DVAL_MAX, unit);

    if (bcs.filename != "" && bcs.variable != "") {
      Teuchos::ParameterList& src = out_list.sublist("fields").sublist("SRC 0");
      src.set<Teuchos::Array<std::string>>("regions", regions)
        .set<std::string>("spatial distribution method", "field")
        .set<bool>("use volume fractions", false);
      src.sublist("field")
        .set<std::string>("field key", bcs.variable)
        .set<std::string>("component", "cell");

      auto& field_ev = glist_->sublist("state").sublist("evaluators").sublist(bcs.variable);
      field_ev.set<std::string>("evaluator type", "independent variable from file")
        .set<std::string>("filename", bcs.filename)
        .set<std::string>("domain name", (domain == "matrix") ? "domain" : domain)
        .set<std::string>("component name", "cell")
        .set<std::string>("mesh entity", "cell")
        .set<std::string>("variable name", bcs.variable)
        .set<int>("number of dofs", 1)
        .set<bool>("constant in time", true);
    } else {
      // additional options for submodels
      double peaceman_r, peaceman_d;
      if (srctype == "peaceman_well") {
        element = static_cast<DOMElement*>(same_list[0]);
        peaceman_r = GetAttributeValueD_(element, "radius", TYPE_NUMERICAL, 0.0, DVAL_MAX, "m");
        peaceman_d = GetAttributeValueD_(element, "depth", TYPE_NUMERICAL, DVAL_MIN, DVAL_MAX, "m");
      }

      // save in the XML files
      Teuchos::ParameterList& src =
        (pkname == "flow") ? out_list.sublist("wells").sublist(srcname) : out_list.sublist(srcname);
      src.set<Teuchos::Array<std::string>>("regions", regions);
      src.set<std::string>("spatial distribution method", weight);
      src.set<bool>("use volume fractions", WeightVolumeSubmodel_(regions));

      Teuchos::ParameterList* srcfn = &src.sublist("well");

      if (srctype == "peaceman_well") {
        srcfn->set<std::string>("submodel", "bhp");
        srcfn->set<double>("well radius", peaceman_r);
        srcfn->set<double>("depth", peaceman_d);
        srcfn = &srcfn->sublist("bhp");
      }

      TranslateGenericMath_(bcs, *srcfn);
    }
  }

  return out_list;
}


/* ******************************************************************
* Add optional fracture model to the list of physical models.
****************************************************************** */
Teuchos::ParameterList
InputConverterU::TranslateFlowFractures_(const std::string& domain)
{
  Teuchos::ParameterList out_list;

  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH)
    *vo_->os() << "Translating optional fracture models, domain=" << domain << std::endl;

  MemoryManager mm;
  DOMNodeList* children;
  DOMNode* node;
  DOMElement* element;

  bool flag;
  if (domain == "matrix") {
    return out_list;
  } else {
    node = GetUniqueElementByTagsString_("fracture_network, materials", flag);
    if (!flag) return out_list;
    element = static_cast<DOMElement*>(node);
  }

  children = element->getElementsByTagName(mm.transcode("material"));

  for (int i = 0; i < children->getLength(); ++i) {
    DOMNode* inode = children->item(i);

    // get assigned regions
    node = GetUniqueElementByTagsString_(inode, "assigned_regions", flag);
    std::vector<std::string> regions = CharToStrings_(mm.transcode(node->getTextContent()));

    // get permeability
    std::string type, model;
    node = GetUniqueElementByTagsString_(inode, "permeability", flag);
    if (flag) model = GetAttributeValueS_(node, "model", "cubic law, linear, constant");

    node = GetUniqueElementByTagsString_(inode, "aperture", flag);
    if (flag) {
      double aperture(0.0);
      type = GetAttributeValueS_(node, "model", TYPE_NONE, false, "");
      if (type == "")
        aperture =
          GetAttributeValueD_(node, "value", TYPE_NUMERICAL, 0.0, DVAL_MAX, "m^2", false, 0.0);

      for (auto it = regions.begin(); it != regions.end(); ++it) {
        std::stringstream ss;
        ss << "FPM for " << *it;

        Teuchos::ParameterList& fam_list = out_list.sublist(ss.str());
        fam_list.set<std::string>("region", *it)
          .set<std::string>("model", model)
          .set<double>("aperture", aperture);
      }
    }

    if (model == "constant") {
      double val = GetAttributeValueD_(node, "value", TYPE_NUMERICAL, 0.0, DVAL_MAX, "", true);

      for (auto it = regions.begin(); it != regions.end(); ++it) {
        std::stringstream ss;
        ss << "FPM for " << *it;
        out_list.sublist(ss.str()).set<double>("value", val);
      }
    }
  }

  return out_list;
}

} // namespace AmanziInput
} // namespace Amanzi
