
/*
  Input Converter

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

// Amanzi's
#include "errors.hh"

#include "InputConverterU.hh"
#include "InputConverterU_Defs.hh"

namespace Amanzi {
namespace AmanziInput {

/* ******************************************************************
* Create multiphase list.
****************************************************************** */
Teuchos::ParameterList InputConverterU::TranslateMultiphase_(const std::string& domain)
{
  Teuchos::ParameterList out_list;

  if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH)
    *vo_->os() << "Translating multiphase, domain=" << domain << std::endl;

  MemoryManager mm;

  bool flag;
  DOMNode* node;

  // solver data
  out_list.set<std::string>("Jacobian type", "analytic")
     .set<std::string>("linear solver", "AMESOS")
     .set<std::string>("preconditioner", "Euclid")
     .set<std::string>("NCP function", "min")
     .set<bool>("CPR enhancement", false);

  std::vector<int> blocks(1, 0);
  std::vector<std::string> pcs(1, "Hypre AMG");
  out_list.sublist("CPR parameters")
     .set<bool>("global solve", true)
     .set<Teuchos::Array<int> >("correction blocks", blocks)
     .set<Teuchos::Array<std::string> >("preconditioner", pcs);

  // chemical species
  out_list.sublist("molecular diffusion") = TranslateMolecularDiffusion_();

  out_list.set<int>("number of aqueous components", phases_["water"].size())
          .set<int>("number of gaseous components", phases_["air"].size())
          .set<double>("molar mass of water", 18.0e-3);

  // water retention models
  out_list.sublist("water retention models") = TranslateWRM_("multiphase");

  // time integrator
  std::string err_options("residual"), nonlinear_solver("newton");
  std::string unstr_controls("unstructured_controls, unstr_multiphase_controls");

  bool modify_correction(false);
  node = GetUniqueElementByTagsString_("unstructured_controls, unstr_nonlinear_solver, modify_correction", flag);

  out_list.sublist("time integrator") = TranslateTimeIntegrator_(
      err_options, nonlinear_solver, modify_correction, unstr_controls,
      TI_TS_REDUCTION_FACTOR, TI_TS_INCREASE_FACTOR);  

  // boundary and initial conditions
  out_list.sublist("boundary conditions") = TranslateMultiphaseBCs_();

  // operators
  std::string disc_method("fv-default, fv-default");
  out_list.sublist("operators") = TranslateDiffusionOperator_(
      disc_method, "", "", "upwind: face", "vapor matrix", true);

  out_list.sublist("operators").sublist("molecular diffusion operator") =
      out_list.sublist("operators").sublist("diffusion operator");

  out_list.sublist("operators").sublist("advection operator")
      .set<std::string>("discretization primary", "upwind")
      .set<int>("reconstruction order", 0);

  return out_list;
}


/* ******************************************************************
* Create list of multiphase BCs.
****************************************************************** */
Teuchos::ParameterList InputConverterU::TranslateMultiphaseBCs_()
{
  Teuchos::ParameterList out_list;

  MemoryManager mm;

  char *text, *tagname;
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
    tagname = mm.transcode(inode->getNodeName());

    // read the assigned regions
    bool flag;
    node = GetUniqueElementByTagsString_(inode, "assigned_regions", flag);
    text = mm.transcode(node->getTextContent());
    std::vector<std::string> regions = CharToStrings_(text);

    vv_bc_regions_.insert(vv_bc_regions_.end(), regions.begin(), regions.end());

    // try two components
    node = GetUniqueElementByTagsString_(inode, "liquid_phase, liquid_component", flag);
    if (!flag) 
      node = GetUniqueElementByTagsString_(inode, "liquid_phase, solute_component", flag);
    if (!flag) continue;

    // process a group of similar elements defined by the first element
    std::string bctype_in, bctype;
    std::vector<DOMNode*> same_list = GetSameChildNodes_(node, bctype_in, flag, true);

    while (same_list.size() > 0) {
      std::string solute_name = GetAttributeValueS_(same_list[0], "name", TYPE_NONE, false);

      std::map<double, double> tp_values;
      std::map<double, std::string> tp_forms;

      for (auto it = same_list.begin(); it != same_list.end(); ++it) {
        std::string tmp_name = GetAttributeValueS_(*it, "name", TYPE_NONE, false);

        if (tmp_name == solute_name) {
          double t0 = GetAttributeValueD_(*it, "start", TYPE_TIME, DVAL_MIN, DVAL_MAX, "s");
          tp_forms[t0] = GetAttributeValueS_(*it, "function");
          GetAttributeValueD_(*it, "value", TYPE_NUMERICAL, DVAL_MIN, DVAL_MAX);
          tp_values[t0] = GetAttributeValueD_(*it, "value", TYPE_NUMERICAL, DVAL_MIN, DVAL_MAX);

          same_list.erase(it);
          it--;
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

      // create names, modify data
      std::string bcname;
      if (bctype_in == "uniform_pressure") {
        bctype = "pressure liquid";
        bcname = "boundary pressure";
      }
      else if (bctype_in == "inward_volumetric_flux") {
        bctype = "mass flux total";
        bcname = "outward mass flux";
        for (int k = 0; k < values.size(); k++) values[k] *= -1;
      }
      else if (bctype_in == "saturation") {
        bctype = "saturation";
        bcname = "boundary saturation";
      }

      std::stringstream ss;
      ss << "BC " << ibc++;

      // save in the XML files  
      Teuchos::ParameterList& tbc_list = out_list.sublist(bctype);
      Teuchos::ParameterList& bc = tbc_list.sublist(ss.str());
      bc.set<Teuchos::Array<std::string> >("regions", regions)
        .set<std::string>("spatial distribution method", "none");
      if (solute_name != "")
        bc.set<std::string>("name", solute_name);

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
  }

  return out_list;
}

}  // namespace AmanziInput
}  // namespace Amanzi
