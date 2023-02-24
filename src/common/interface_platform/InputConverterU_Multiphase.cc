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
  Errors::Message msg;
  Teuchos::ParameterList out_list;
  multiphase_ = true;

  if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH)
    *vo_->os() << "Translating multiphase, domain=" << domain << std::endl;

  MemoryManager mm;

  bool flag;
  DOMNode* node;

  // checks
  if (eos_model_ == "") {
    msg << "EOS model is required.\n";
    Exceptions::amanzi_throw(msg);
  }

  // header
  out_list.set<std::string>("domain name", (domain == "matrix") ? "domain" : domain);

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

  out_list.set<int>("number of aqueous components", phases_[LIQUID].dissolved.size())
    .set<int>("number of gaseous components", phases_[GAS].dissolved.size())
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
                                                                 dt_cut_["transient"],
                                                                 dt_inc_["transient"]);

  // boundary and initial conditions
  out_list.sublist("boundary conditions") = TranslateMultiphaseBCs_();

  // operators
  std::string disc_method("fv-default, fv-default");
  out_list.sublist("operators") =
    TranslateDiffusionOperator_(disc_method, "", "", "upwind: face", "vapor matrix", true);

  auto& tmp1 = out_list.sublist("operators").sublist("diffusion operator");
  auto& tmp2 = out_list.sublist("operators").sublist("molecular diffusion operator");
  tmp2.sublist("matrix") = tmp1.sublist("matrix");
  tmp2.sublist("preconditioner") = tmp1.sublist("preconditioner");

  tmp2.sublist("matrix").set<bool>("gravity", false);
  tmp2.sublist("preconditioner").set<bool>("gravity", false);

  out_list.sublist("operators")
    .sublist("advection operator")
    .set<std::string>("discretization primary", "upwind")
    .set<int>("reconstruction order", 0);

  // -- upwind 
  tmp1.sublist("upwind").set<std::string>("upwind method", "upwind: darcy velocity");
  tmp1.sublist("upwind").sublist("upwind standard parameters").set<double>("tolerance", 1e-12);

  // additional state list
  auto& fic = state_list.sublist("initial conditions");
  auto& fev = state_list.sublist("evaluators");

  // -- keys
  Key pressure_liquid_key = Keys::getKey(domain, "pressure_liquid");
  Key pressure_gas_key = Keys::getKey(domain, "pressure_gas");
  Key pressure_vapor_key = Keys::getKey(domain, "pressure_vapor");

  Key sat_liquid_key = Keys::getKey(domain, "saturation_liquid");
  Key sat_gas_key = Keys::getKey(domain, "saturation_gas");

  Key porosity_key = Keys::getKey(domain, "porosity");

  Key mol_density_liquid_key = Keys::getKey(domain, "molar_density_liquid");
  Key mol_density_gas_key = Keys::getKey(domain, "molar_density_gas");

  Key mass_density_liquid_key = Keys::getKey(domain, "mass_density_liquid");

  Key mole_xl_key = Keys::getKey(domain, "mole_fraction_liquid");
  Key mole_xg_key = Keys::getKey(domain, "mole_fraction_gas");
  Key mole_xv_key = Keys::getKey(domain, "mole_fraction_vapor");

  Key viscosity_liquid_key = Keys::getKey(domain, "viscosity_liquid");
  Key viscosity_gas_key = Keys::getKey(domain, "viscosity_gas");

  Key temperature_key = Keys::getKey(domain, "temperature");
  Key enthalpy_liquid_key = Keys::getKey(domain, "enthalpy_liquid");
  Key enthalpy_gas_key = Keys::getKey(domain, "enthalpy_gas");
  Key ie_liquid_key = Keys::getKey(domain, "internal_energy_liquid");
  Key ie_gas_key = Keys::getKey(domain, "internal_energy_gas");
  Key ie_rock_key = Keys::getKey(domain, "internal_energy_rock");

  Key diff_liquid_key = Keys::getKey(domain, "diffusion_liquid");
  Key diff_gas_key = Keys::getKey(domain, "diffusion_gas");
  Key diff_vapor_key = Keys::getKey(domain, "diffusion_vapor");
  Key conductivity_key = Keys::getKey(domain, "thermal_conductivity");

  Key adv_liquid_key = Keys::getKey(domain, "advection_liquid");
  Key adv_gas_key = Keys::getKey(domain, "advection_gas");
  Key adv_water_key = Keys::getKey(domain, "advection_water");
  Key adv_energy_liquid_key = Keys::getKey(domain, "adv_energy_liquid");
  Key adv_energy_gas_key = Keys::getKey(domain, "adv_energy_gas");

  Key mol_diff_liquid_key = Keys::getKey(domain, "molecular_diff_liquid");
  Key mol_diff_gas_key = Keys::getKey(domain, "molecular_diff_gas");

  Key storage_tcc_key = Keys::getKey(domain, "total_component_storage");
  Key storage_water_key = Keys::getKey(domain, "total_water_storage");
  Key energy_key = Keys::getKey(domain, "energy");

  Key ncp_f_key = Keys::getKey(domain, "ncp_f");
  Key ncp_g_key = Keys::getKey(domain, "ncp_g");

  Key relperm_liquid_key = Keys::getKey(domain, "rel_permeability_liquid");
  Key relperm_gas_key = Keys::getKey(domain, "rel_permeability_gas");

  // -- density
  auto& tmp = fev.sublist(mol_density_gas_key);
  tmp.set<std::string>("evaluator type", "eos")
    .set<std::string>("eos basis", "molar")
    .set<std::string>("molar density key", mol_density_gas_key)
    .set<std::string>("pressure key", pressure_gas_key);
  tmp.sublist("EOS parameters")
    .set<std::string>("eos type", "ideal gas")
    .set<double>("molar mass of gas", 28.9647e-03); // dry air (not used ?)

  fev.sublist(mol_density_liquid_key).set<std::string>("pressure key", pressure_liquid_key);

  // -- viscosity
  double viscosity = fic.sublist("const_fluid_viscosity").template get<double>("value");
  AddIndependentFieldEvaluator_(fev, viscosity_liquid_key, "All", "cell", viscosity);

  // -- diffusion
  auto diff_l = out_list.sublist("molecular diffusion")
      .get<Teuchos::Array<double>>("aqueous values").toVector();
  auto diff_g = out_list.sublist("molecular diffusion")
      .get<Teuchos::Array<double>>("gaseous values").toVector();

  if (diff_l.size() == 0 || diff_l.size() != diff_g.size()) {
    msg << "Incomplete definition of species: #aqueous=" << diff_l.size() 
        << ", #gaseous=" << diff_g.size() << ".\n";
    Exceptions::amanzi_throw(msg);
  }

  AddIndependentFieldEvaluator_(fev, mol_diff_liquid_key, "All", "cell", diff_l[0]);
  AddIndependentFieldEvaluator_(fev, mol_diff_gas_key, "All", "cell", diff_g[0]);

  // -- pressure (why do we use IC? FIXME)
  fic.sublist(pressure_liquid_key) = fic.sublist(Keys::getKey(domain, "pressure"));
  fic.remove(Keys::getKey(domain, "pressure"));

  // -- temperature
  if (isothermal_) {
    fev.sublist(temperature_key) = fic.sublist(temperature_key);
    fev.sublist(temperature_key).set<std::string>("evaluator type", "independent variable");
    fic.remove(temperature_key);
  }

  // system of equations (fixed at the moment)
  // -- equations
  std::string pl(pressure_liquid_key), pg(pressure_gas_key), xg(mole_xg_key);
  std::vector<double> ones({ 1.0, 1.0 });

  // liquid phase: primary component
  Teuchos::ParameterList& peqn = out_list.sublist("system").sublist("pressure eqn");
  peqn.set<std::string>("primary unknown", pl)
    .set<Teuchos::Array<std::string>>("advection liquid",
                                      std::vector<std::string>({ adv_water_key, pl }))
    .set<Teuchos::Array<double>>("advection factors", ones)
    .set<Teuchos::Array<double>>("diffusion factors", ones)
    .set<std::string>("accumulation", storage_water_key);
  if (phases_[GAS].model != "") {
    peqn.set<Teuchos::Array<std::string>>(
      "diffusion gas", std::vector<std::string>({ diff_vapor_key, mole_xv_key }));

    fev.sublist(storage_water_key)
      .set<std::string>("evaluator type", "storage water")
      .set<std::string>("molar density liquid key", mol_density_liquid_key)
      .set<std::string>("molar density gas key", mol_density_gas_key)
      .set<std::string>("porosity key", porosity_key)
      .set<std::string>("saturation liquid key", sat_liquid_key)
      .set<std::string>("mole fraction vapor key", mole_xv_key)
      .set<std::string>("tag", "");

    fev.sublist(mole_xv_key)
      .set<std::string>("evaluator type", "product")
      .set<Teuchos::Array<std::string>>(
        "dependencies", std::vector<std::string>({ pressure_gas_key, pressure_vapor_key }))
      .set<Teuchos::Array<int>>("powers", std::vector<int>({ -1, 1 }))
      .set<std::string>("tag", "");
  } else {
    fev.sublist(storage_water_key)
      .set<std::string>("evaluator type", "product")
      .set<Teuchos::Array<std::string>>(
        "dependencies", std::vector<std::string>( { mol_density_liquid_key, sat_liquid_key, porosity_key }))
      .set<Teuchos::Array<int>>("powers", std::vector<int>({ 1, 1, 1 }))
      .set<std::string>("tag", "");

    AddIndependentFieldEvaluator_(fev, mole_xv_key, "All", "cell", 0.0);
  }

  Teuchos::ParameterList& seqn = out_list.sublist("system").sublist("solute eqn");
  seqn.set<std::string>("primary unknown", xg)
    .set<Teuchos::Array<std::string>>("advection liquid",
                                      std::vector<std::string>({ adv_liquid_key, pl }))
    .set<Teuchos::Array<std::string>>("advection gas",
                                      std::vector<std::string>({ adv_gas_key, pg }))
    .set<Teuchos::Array<double>>("advection factors", ones)
    .set<Teuchos::Array<std::string>>("diffusion liquid",
                                      std::vector<std::string>({ diff_liquid_key, mole_xl_key }))
    .set<Teuchos::Array<std::string>>("diffusion gas",
                                      std::vector<std::string>({ diff_gas_key, mole_xg_key }))
    .set<Teuchos::Array<double>>("diffusion factors", ones)
    .set<std::string>("accumulation", storage_tcc_key);


  Teuchos::ParameterList& ceqn = out_list.sublist("system").sublist("constraint eqn");
  ceqn.set<std::string>("primary unknown", sat_liquid_key)
    .set<Teuchos::Array<std::string>>("ncp evaluators",
                                      std::vector<std::string>({ ncp_f_key, ncp_g_key }));

  if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH) {
    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "primary unknwons:" << pl << ", " << xg << ", " << sat_liquid_key << std::endl;
  }

  if (isothermal_ == false) {
    Teuchos::ParameterList& teqn = out_list.sublist("system").sublist("energy eqn");
    teqn.set<std::string>("primary unknown", temperature_key)
      .set<Teuchos::Array<std::string>>("advection liquid",
                                        std::vector<std::string>({ adv_energy_liquid_key, pl }))
      // .set<Teuchos::Array<std::string>>("advection gas",
      //                                   std::vector<std::string>({ adv_energy_gas_key, pg }))
      .set<Teuchos::Array<double>>("advection factors", ones)

      .set<Teuchos::Array<std::string>>("diffusion liquid",
                                        std::vector<std::string>({ conductivity_key, temperature_key }))
      .set<Teuchos::Array<double>>("diffusion factors", ones)
      .set<std::string>("accumulation", energy_key);

    Teuchos::ParameterList& thermal =
      out_list.sublist("thermal conductivity evaluator").sublist("thermal conductivity parameters");
    thermal.set<std::string>("thermal conductivity type", "two-phase Peters-Lidard");
    thermal.set<double>("thermal conductivity of gas", 0.02);
    thermal.set<double>("unsaturated alpha", 1.0);
    thermal.set<double>("epsilon", 1.0e-10);

    double cv_f(0.606), cv_r(0.2);
    node = GetUniqueElementByTagsString_("materials", flag);
    std::vector<DOMNode*> materials = GetChildren_(node, "material", flag);

    node = GetUniqueElementByTagsString_(materials[0], "thermal_properties, rock_conductivity", flag);
    if (flag) cv_f = GetTextContentD_(node, "W/m/K", true);

    node = GetUniqueElementByTagsString_(materials[0], "thermal_properties, rock_conductivity", flag);
    if (flag) cv_r = GetTextContentD_(node, "W/m/K", true);

    thermal.set<double>("thermal conductivity of liquid", cv_f);
    thermal.set<double>("thermal conductivity of rock", cv_r);
    thermal.set<double>("reference temperature", 298.15);

    if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH) {
      Teuchos::OSTab tab = vo_->getOSTab();
      *vo_->os() << "primary unknwons:" << temperature_key << std::endl;
    }
  }

  // -- evaluators
  auto evals = Teuchos::Array<std::string>({ ncp_g_key,
                                             storage_water_key,
                                             storage_tcc_key,
                                             adv_water_key,
                                             adv_liquid_key,
                                             adv_gas_key,
                                             diff_liquid_key,
                                             diff_gas_key,
                                             diff_vapor_key,
                                             sat_gas_key,
                                             mole_xl_key,
                                             mole_xv_key });
  if (isothermal_ == false) {
    evals.push_back(adv_energy_liquid_key);
    // evals.push_back(adv_energy_gas_key);
  }
  out_list.set<Teuchos::Array<std::string>>("evaluators", evals);

  fev.sublist(ncp_g_key)
    .set<std::string>("evaluator type", "ncp mole fraction gas")
    .set<std::string>("mole fraction vapor key", mole_xv_key)
    .set<std::string>("mole fraction gas key", mole_xg_key)
    .set<std::string>("tag", "");

  fev.sublist(storage_tcc_key)
    .set<std::string>("evaluator type", "storage component")
    .set<std::string>("saturation liquid key", sat_liquid_key)
    .set<std::string>("porosity key", porosity_key)
    .set<std::string>("molar density liquid key", mol_density_liquid_key)
    .set<std::string>("molar density gas key", mol_density_gas_key)
    .set<std::string>("mole fraction liquid key", mole_xl_key)
    .set<std::string>("mole fraction gas key", mole_xg_key)
    .set<std::string>("tag", "");

  fev.sublist(adv_water_key)
    .set<std::string>("evaluator type", "product")
    .set<Teuchos::Array<std::string>>(
      "dependencies",
      std::vector<std::string>(
        { mol_density_liquid_key, relperm_liquid_key, viscosity_liquid_key }))
    .set<Teuchos::Array<int>>("powers", std::vector<int>({ 1, 1, -1 }))
    .set<std::string>("tag", "");

  fev.sublist(adv_liquid_key)
    .set<std::string>("evaluator type", "product")
    .set<Teuchos::Array<std::string>>(
      "dependencies",
      std::vector<std::string>(
        { mol_density_liquid_key, mole_xl_key, relperm_liquid_key, viscosity_liquid_key }))
    .set<Teuchos::Array<int>>("powers", std::vector<int>({ 1, 1, 1, -1 }))
    .set<std::string>("tag", "");

  fev.sublist(adv_gas_key)
    .set<std::string>("evaluator type", "product")
    .set<Teuchos::Array<std::string>>(
      "dependencies",
      std::vector<std::string>(
        { mol_density_gas_key, mole_xg_key, relperm_gas_key, viscosity_gas_key }))
    .set<Teuchos::Array<int>>("powers", std::vector<int>({ 1, 1, 1, -1 }))
    .set<std::string>("tag", "");

  fev.sublist(diff_liquid_key)
    .set<std::string>("evaluator type", "product")
    .set<Teuchos::Array<std::string>>(
      "dependencies",
      std::vector<std::string>({ mol_diff_liquid_key, mol_density_liquid_key, porosity_key, sat_liquid_key }))
    .set<Teuchos::Array<int>>("powers", std::vector<int>({ 1, 1, 1, 1 }))
    .set<std::string>("tag", "");

  fev.sublist(diff_gas_key)
    .set<std::string>("evaluator type", "product")
    .set<Teuchos::Array<std::string>>(
      "dependencies",
      std::vector<std::string>({ mol_diff_gas_key, mol_density_gas_key, sat_gas_key }))
    .set<Teuchos::Array<int>>("powers", std::vector<int>({ 1, 1, 1 }))
    .set<std::string>("tag", "");

  fev.sublist(diff_vapor_key)
    .set<std::string>("evaluator type", "product")
    .set<Teuchos::Array<std::string>>(
      "dependencies",
      std::vector<std::string>(
        { mol_diff_gas_key, mol_density_gas_key, porosity_key, sat_gas_key }))
    .set<Teuchos::Array<int>>("powers", std::vector<int>({ 1, 1, 1, 1 }))
    .set<std::string>("tag", "");

  fev.sublist(sat_gas_key)
    .set<std::string>("evaluator type", "saturation gas")
    .set<std::string>("saturation liquid key", sat_liquid_key);

  fev.sublist(mole_xl_key)
    .set<std::string>("evaluator type", "mole fraction liquid")
    .set<std::string>("pressure gas key", pressure_gas_key)
    .set<std::string>("mole fraction gas key", mole_xg_key)
    .set<std::string>("tag", "");

  if (isothermal_ == false) {
    fev.sublist(adv_energy_liquid_key)
      .set<std::string>("evaluator type", "product")
      .set<Teuchos::Array<std::string>>("dependencies",
        std::vector<std::string>( { mol_density_liquid_key, relperm_liquid_key, 
                                    enthalpy_liquid_key, viscosity_liquid_key }))
      .set<Teuchos::Array<int>>("powers", std::vector<int>({ 1, 1, 1, -1 }))
      .set<std::string>("tag", "");

    fev.sublist(adv_energy_gas_key)
      .set<std::string>("evaluator type", "product")
      .set<Teuchos::Array<std::string>>("dependencies",
        std::vector<std::string>( { mol_density_gas_key, relperm_gas_key, 
                                    enthalpy_gas_key, viscosity_gas_key }))
      .set<Teuchos::Array<int>>("powers", std::vector<int>({ 1, 1, 1, -1 }))
      .set<std::string>("tag", "");

    fev.sublist(ie_liquid_key).set<std::string>("pressure key", pl);
    fev.sublist(ie_gas_key).set<std::string>("pressure key", pg);
    fev.sublist(ie_rock_key).set<std::string>("pressure key", pl);
  }

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
      if (bctype_in == "uniform_pressure" || bctype_in == "linear_pressure") {
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
      if (bctype_in == "linear_pressure") {
        double refv;
        std::vector<double> grad, refc;
        auto element = static_cast<DOMElement*>(same_list[0]);

        refv = GetAttributeValueD_(element, "value", TYPE_NUMERICAL, DVAL_MIN, DVAL_MAX, "Pa");
        grad = GetAttributeVectorD_(element, "gradient", dim_, "Pa/m");
        refc = GetAttributeVectorD_(element, "reference_coord", dim_, "m");
        grad.insert(grad.begin(), 0.0);
        refc.insert(refc.begin(), 0.0);

        bcfn.sublist("function-linear")
          .set<double>("y0", refv)
          .set<Teuchos::Array<double>>("x0", refc)
          .set<Teuchos::Array<double>>("gradient", grad);
      } else if (times.size() == 1) {
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
