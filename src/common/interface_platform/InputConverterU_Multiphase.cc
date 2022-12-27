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

// Amanzi's
#include "errors.hh"

#include "InputConverterU.hh"
#include "InputConverterU_Defs.hh"

namespace Amanzi {
namespace AmanziInput {

/* ******************************************************************
* Create multiphase list.
****************************************************************** */
Teuchos::ParameterList
InputConverterU::TranslateMultiphase_(const std::string& domain, Teuchos::ParameterList& state_list)
{
  Teuchos::ParameterList out_list;
  multiphase_ = true;

  if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH)
    *vo_->os() << "Translating multiphase, domain=" << domain << std::endl;

  MemoryManager mm;

  bool flag;
  DOMNode* node;

  // solver data
  out_list.set<std::string>("Jacobian type", "analytic")
    .set<std::string>("linear solver", "GMRES for Newton-0")
    .set<std::string>("preconditioner", "Euclid")
    .set<std::string>("NCP function", "min")
    .set<bool>("CPR enhancement", false);

  std::vector<int> blocks(1, 0);
  std::vector<std::string> pcs(1, "Hypre AMG");
  out_list.sublist("CPR parameters")
    .set<bool>("global solve", true)
    .set<Teuchos::Array<int>>("correction blocks", blocks)
    .set<Teuchos::Array<std::string>>("preconditioner", pcs);

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
  node = GetUniqueElementByTagsString_(
    "unstructured_controls, unstr_nonlinear_solver, modify_correction", flag);

  out_list.sublist("time integrator") = TranslateTimeIntegrator_(err_options,
                                                                 nonlinear_solver,
                                                                 modify_correction,
                                                                 unstr_controls,
                                                                 TI_SOLVER,
                                                                 TI_TS_REDUCTION_FACTOR,
                                                                 TI_TS_INCREASE_FACTOR);

  // boundary and initial conditions
  out_list.sublist("boundary conditions") = TranslateMultiphaseBCs_();

  // operators
  std::string disc_method("fv-default, fv-default");
  out_list.sublist("operators") =
    TranslateDiffusionOperator_(disc_method, "", "", "upwind: face", "vapor matrix", true);

  out_list.sublist("operators").sublist("molecular diffusion operator") =
    out_list.sublist("operators").sublist("diffusion operator");

  out_list.sublist("operators")
    .sublist("advection operator")
    .set<std::string>("discretization primary", "upwind")
    .set<int>("reconstruction order", 0);

  // additional state list
  auto& fic = state_list.sublist("initial conditions");
  auto& fev = state_list.sublist("evaluators");

  // -- density
  auto& tmp = fev.sublist("molar_density_gas");
  tmp.set<std::string>("evaluator type", "eos")
    .set<std::string>("eos basis", "molar")
    .set<std::string>("molar density key", "molar_density_gas")
    .set<std::string>("pressure key", "pressure_gas");
  tmp.sublist("EOS parameters")
    .set<std::string>("eos type", "ideal gas")
    .set<double>("molar mass of gas", 28.9647e-03); // dry air (not used ?)

  AddIndependentFieldEvaluator_(
    fev, "molar_density_liquid", "All", "cell", rho_ / MOLAR_MASS_WATER);

  // -- viscosity
  double viscosity = fic.sublist("const_fluid_viscosity").template get<double>("value");
  AddIndependentFieldEvaluator_(fev, "viscosity_liquid", "All", "cell", viscosity);

  // -- diffusion
  auto diff = out_list.sublist("molecular diffusion")
                .get<Teuchos::Array<double>>("aqueous values")
                .toVector();
  AddIndependentFieldEvaluator_(fev, "molecular_diff_liquid", "All", "cell", diff[0]);

  diff = out_list.sublist("molecular diffusion")
           .get<Teuchos::Array<double>>("gaseous values")
           .toVector();
  AddIndependentFieldEvaluator_(fev, "molecular_diff_gas", "All", "cell", diff[0]);

  // -- pressure
  fic.sublist("pressure_liquid") = fic.sublist("pressure");
  fic.remove("pressure");

  // -- saturation
  fic.sublist("saturatiob_liquid") = fic.sublist("saturation");
  fic.remove("saturation");

  // -- temperature
  fev.sublist("temperature") = fic.sublist("temperature");
  fev.sublist("temperature").set<std::string>("evaluator type", "independent variable");
  fic.remove("temperature");

  // system of equations (fixed at the moment)
  // -- equations
  std::string pl("pressure_liquid"), xg("mole_fraction_gas");

  auto evals = Teuchos::Array<std::string>({ "ncp_g",
                                             "total_water_storage",
                                             "total_component_storage",
                                             "advection_water",
                                             "advection_liquid",
                                             "advection_gas",
                                             "diffusion_liquid",
                                             "diffusion_gas",
                                             "diffusion_vapor",
                                             "saturation_gas",
                                             "mole_fraction_liquid",
                                             "mole_fraction_vapor" });
  out_list.set<Teuchos::Array<std::string>>("evaluators", evals);

  Teuchos::ParameterList& peqn = out_list.sublist("system").sublist("pressure eqn");
  peqn.set<std::string>("primary unknown", pl)
    .set<Teuchos::Array<std::string>>("advection liquid",
                                      std::vector<std::string>({ "advection_water", pl }))
    .set<Teuchos::Array<double>>("advection factors", std::vector<double>({ 1.0, 1.0 }))
    .set<Teuchos::Array<std::string>>(
      "diffusion liquid", std::vector<std::string>({ "diffusion_vapor", "mole_fraction_vapor" }))
    .set<Teuchos::Array<double>>("diffusion factors", std::vector<double>({ 1.0, 1.0 }))
    .set<std::string>("accumulation", "total_water_storage");

  Teuchos::ParameterList& seqn = out_list.sublist("system").sublist("solute eqn");
  seqn.set<std::string>("primary unknown", xg)
    .set<Teuchos::Array<std::string>>("advection liquid",
                                      std::vector<std::string>({ "advection_liquid", pl }))
    .set<Teuchos::Array<std::string>>("advection gas",
                                      std::vector<std::string>({ "advection_gas", "pressure_gas" }))
    .set<Teuchos::Array<double>>("advection factors", std::vector<double>({ 1.0, 1.0 }))
    .set<Teuchos::Array<std::string>>(
      "diffusion liquid", std::vector<std::string>({ "diffusion_liquid", "mole_fraction_liquid" }))
    .set<Teuchos::Array<std::string>>(
      "diffusion gas", std::vector<std::string>({ "diffusion_gas", "mole_fraction_gas" }))
    .set<Teuchos::Array<double>>("diffusion factors", std::vector<double>({ 1.0, 1.0 }))
    .set<std::string>("accumulation", "total_component_storage");


  Teuchos::ParameterList& ceqn = out_list.sublist("system").sublist("constraint eqn");
  ceqn.set<std::string>("primary unknown", "saturation_liquid")
    .set<Teuchos::Array<std::string>>("ncp evaluators",
                                      std::vector<std::string>({ "ncp_f", "ncp_g" }));

  // -- evaluators
  fev.sublist("ncp_g")
    .set<std::string>("evaluator type", "ncp mole fraction gas")
    .set<std::string>("mole fraction vapor key", "mole_fraction_vapor")
    .set<std::string>("mole fraction gas key", "mole_fraction_gas")
    .set<std::string>("tag", "");

  fev.sublist("total_water_storage")
    .set<std::string>("evaluator type", "storage water")
    .set<std::string>("molar density liquid key", "molar_density_liquid")
    .set<std::string>("molar density gas key", "molar_density_gas")
    .set<std::string>("porosity key", "porosity")
    .set<std::string>("saturation liquid key", "saturation_liquid")
    .set<std::string>("mole fraction vapor key", "mole_fraction_vapor")
    .set<std::string>("tag", "");

  fev.sublist("total_component_storage")
    .set<std::string>("evaluator type", "storage component")
    .set<std::string>("saturation liquid key", "saturation_liquid")
    .set<std::string>("porosity key", "porosity")
    .set<std::string>("molar density liquid key", "molar_density_liquid")
    .set<std::string>("molar density gas key", "molar_density_gas")
    .set<std::string>("mole fraction liquid key", "mole_fraction_liquid")
    .set<std::string>("mole fraction gas key", "mole_fraction_gas")
    .set<std::string>("tag", "");

  fev.sublist("advection_water")
    .set<std::string>("evaluator type", "product")
    .set<Teuchos::Array<std::string>>(
      "dependencies",
      std::vector<std::string>(
        { "molar_density_liquid", "rel_permeability_liquid", "viscosity_liquid" }))
    .set<Teuchos::Array<int>>("powers", std::vector<int>({ 1, 1, -1 }))
    .set<std::string>("tag", "");

  fev.sublist("advection_liquid")
    .set<std::string>("evaluator type", "product")
    .set<Teuchos::Array<std::string>>("dependencies",
                                      std::vector<std::string>({ "molar_density_liquid",
                                                                 "mole_fraction_liquid",
                                                                 "rel_permeability_liquid",
                                                                 "viscosity_liquid" }))
    .set<Teuchos::Array<int>>("powers", std::vector<int>({ 1, 1, 1, -1 }))
    .set<std::string>("tag", "");

  fev.sublist("advection_gas")
    .set<std::string>("evaluator type", "product")
    .set<Teuchos::Array<std::string>>(
      "dependencies",
      std::vector<std::string>(
        { "molar_density_gas", "mole_fraction_gas", "rel_permeability_gas", "viscosity_gas" }))
    .set<Teuchos::Array<int>>("powers", std::vector<int>({ 1, 1, 1, -1 }))
    .set<std::string>("tag", "");

  fev.sublist("diffusion_liquid")
    .set<std::string>("evaluator type", "product")
    .set<Teuchos::Array<std::string>>(
      "dependencies",
      std::vector<std::string>(
        { "molecular_diff_liquid", "molar_density_liquid", "saturation_liquid" }))
    .set<Teuchos::Array<int>>("powers", std::vector<int>({ 1, 1, 1 }))
    .set<std::string>("tag", "");

  fev.sublist("diffusion_gas")
    .set<std::string>("evaluator type", "product")
    .set<Teuchos::Array<std::string>>(
      "dependencies",
      std::vector<std::string>({ "molecular_diff_gas", "molar_density_gas", "saturation_gas" }))
    .set<Teuchos::Array<int>>("powers", std::vector<int>({ 1, 1, 1 }))
    .set<std::string>("tag", "");

  fev.sublist("diffusion_vapor")
    .set<std::string>("evaluator type", "product")
    .set<Teuchos::Array<std::string>>(
      "dependencies",
      std::vector<std::string>(
        { "molecular_diff_gas", "molar_density_gas", "porosity", "saturation_gas" }))
    .set<Teuchos::Array<int>>("powers", std::vector<int>({ 1, 1, 1, 1 }))
    .set<std::string>("tag", "");

  fev.sublist("saturation_gas")
    .set<std::string>("evaluator type", "saturation gas")
    .set<std::string>("saturation liquid key", "saturation_liquid");

  fev.sublist("mole_fraction_liquid")
    .set<std::string>("evaluator type", "mole fraction liquid")
    .set<std::string>("pressure gas key", "pressure_gas")
    .set<std::string>("mole fraction gas key", "mole_fraction_gas")
    .set<std::string>("tag", "");

  fev.sublist("mole_fraction_vapor")
    .set<std::string>("evaluator type", "product")
    .set<Teuchos::Array<std::string>>(
      "dependencies", std::vector<std::string>({ "pressure_gas", "pressure_vapor" }))
    .set<Teuchos::Array<int>>("powers", std::vector<int>({ -1, 1 }))
    .set<std::string>("tag", "");

  return out_list;
}


/* ******************************************************************
* Create list of multiphase BCs.
****************************************************************** */
Teuchos::ParameterList
InputConverterU::TranslateMultiphaseBCs_()
{
  Teuchos::ParameterList out_list;

  MemoryManager mm;

  char *text, *tagname;
  DOMNodeList *node_list, *children;
  DOMNode* node;

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
    if (!flag) node = GetUniqueElementByTagsString_(inode, "liquid_phase, solute_component", flag);
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
      } else if (bctype_in == "inward_volumetric_flux") {
        bctype = "mass flux total";
        bcname = "outward mass flux";
        for (int k = 0; k < values.size(); k++) values[k] *= -1;
      } else if (bctype_in == "saturation") {
        bctype = "saturation";
        bcname = "boundary saturation";
      }

      std::stringstream ss;
      ss << "BC " << ibc++;

      // save in the XML files
      Teuchos::ParameterList& tbc_list = out_list.sublist(bctype);
      Teuchos::ParameterList& bc = tbc_list.sublist(ss.str());
      bc.set<Teuchos::Array<std::string>>("regions", regions)
        .set<std::string>("spatial distribution method", "none");
      if (solute_name != "") bc.set<std::string>("name", solute_name);

      Teuchos::ParameterList& bcfn = bc.sublist(bcname);
      if (times.size() == 1) {
        bcfn.sublist("function-constant").set<double>("value", values[0]);
      } else {
        bcfn.sublist("function-tabular")
          .set<Teuchos::Array<double>>("x values", times)
          .set<Teuchos::Array<double>>("y values", values)
          .set<Teuchos::Array<std::string>>("forms", forms);
      }
    }
  }

  return out_list;
}

} // namespace AmanziInput
} // namespace Amanzi
