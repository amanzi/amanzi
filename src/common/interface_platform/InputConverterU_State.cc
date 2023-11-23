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

#include <sstream>
#include <string>
#include <climits>

// TPLs
#include <xercesc/dom/DOM.hpp>

// Amanzi's
#include "errors.hh"
#include "exceptions.hh"
#include "dbc.hh"
#include "Key.hh"

#include "InputConverterU.hh"
#include "InputConverterU_Defs.hh"

namespace Amanzi {
namespace AmanziInput {

XERCES_CPP_NAMESPACE_USE

/* ******************************************************************
* STATE sublist
****************************************************************** */
Teuchos::ParameterList
InputConverterU::TranslateState_()
{
  Teuchos::ParameterList out_list;

  if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH) { *vo_->os() << "Translating state" << std::endl; }

  // first we write initial conditions for scalars and vectors, not region-specific
  Teuchos::ParameterList& out_ev = out_list.sublist("evaluators");
  Teuchos::ParameterList& out_ic = out_list.sublist("initial conditions");

  MemoryManager mm;
  Errors::Message msg;
  char* tagname;
  char* text_content;

  // --- eos lookup table
  bool flag;
  DOMNode* node = GetUniqueElementByTagsString_("phases, liquid_phase, eos", flag);
  if (flag) {
    eos_model_ = GetTextContentS_(node, "", false);
    if (eos_model_ == "false" || eos_model_ == "False") {
      eos_model_.clear();
    } else if (eos_model_ != "FEHM" && eos_model_ != "0-30C") {
      eos_lookup_table_ = eos_model_;
      eos_model_ = "tabular";
    }
  }

  // --- gravity
  Teuchos::Array<double> gravity(dim_);
  for (int i = 0; i != dim_ - 1; ++i) gravity[i] = 0.0;
  gravity[dim_ - 1] = -const_gravity_;
  out_ic.sublist("gravity").set<Teuchos::Array<double>>("value", gravity);

  double viscosity(0.0);
  rho_ = 1000.0;
  if (phases_[LIQUID].active) {
    // --- constant viscosities
    node = GetUniqueElementByTagsString_("phases, liquid_phase, viscosity", flag, true);
    viscosity = GetTextContentD_(node, "Pa*s");
    out_ic.sublist("const_fluid_viscosity").set<double>("value", viscosity);

    // --- constant density
    node = GetUniqueElementByTagsString_("phases, liquid_phase, density", flag, true);
    rho_ = GetTextContentD_(node, "kg/m^3");
    out_ic.sublist("const_fluid_density").set<double>("value", rho_);

    // --- constant compressibility
    node = GetUniqueElementByTagsString_("phases, liquid_phase, compressibility", flag, false);
    if (flag) {
      beta_ = GetTextContentD_(node, "Pa^-1");
      out_ic.sublist("const_fluid_compressibility").set<double>("value", beta_);
    }
  }

  if (eos_model_ == "") {
    AddIndependentFieldEvaluator_(out_ev, "mass_density_liquid", "All", "*", rho_);
    AddIndependentFieldEvaluator_(
      out_ev, "molar_density_liquid", "All", "*", rho_ / 0.0180153333333);
    AddIndependentFieldEvaluator_(out_ev, "viscosity_liquid", "All", "*", viscosity);
  }

  // --- region specific initial conditions from material properties
  std::map<std::string, int> reg2mat;
  int mat(0);

  // primary continuum
  TranslateCommonContinuumFields_("domain", out_ic, out_ev);

  DOMNodeList* node_list = doc_->getElementsByTagName(mm.transcode("materials"));
  DOMNodeList* children = node_list->item(0)->getChildNodes();

  for (int i = 0; i < children->getLength(); i++) {
    DOMNode* inode = children->item(i);
    if (DOMNode::ELEMENT_NODE == inode->getNodeType()) {
      std::string mat_name = GetAttributeValueS_(inode, "name");
      int mat_id = GetAttributeValueL_(inode, "id", TYPE_NUMERICAL, 0, INT_MAX, false, -1);

      node = GetUniqueElementByTagsString_(inode, "assigned_regions", flag);
      std::vector<std::string> regions = CharToStrings_(mm.transcode(node->getTextContent()));

      // record the material ID for each region that this material occupies
      for (int k = 0; k < regions.size(); k++) {
        if (reg2mat.find(regions[k]) == reg2mat.end()) {
          reg2mat[regions[k]] = mat++;
        } else {
          std::stringstream ss;
          ss << "There is more than one material assigned to region \"" << regions[k] << "\".";
          Exceptions::amanzi_throw(Errors::Message(ss.str().c_str()));
        }
      }

      // record user material ID
      for (int k = 0; k < regions.size(); k++) {
        material_regions_.push_back(regions[k]);
        material_names_.push_back(mat_name);
        material_ids_.push_back(mat_id);
      }

      // create regions string
      std::string reg_str = CreateNameFromVector_(regions);

      // -- transport porosity
      node =
        GetUniqueElementByTagsString_(inode, "mechanical_properties, transport_porosity", flag);
      if (flag) {
        use_transport_porosity_ = true;
        TranslateFieldEvaluator_(node, "transport_porosity", "-", reg_str, regions, out_ic, out_ev);
      } else if (use_transport_porosity_) {
        msg << "Transport porosity element must be specified for all materials or none.";
        Exceptions::amanzi_throw(msg);
      }

      // -- permeability. We parse matrix and fractures together.
      bool perm_err(false), conductivity(false);
      std::string perm_file, perm_format, unit("m^2");

      node = GetUniqueElementByTagsString_(inode, "permeability", flag);
      if (!flag) {
        conductivity = true;
        node = GetUniqueElementByTagsString_(inode, "hydraulic_conductivity", flag);
        unit = "m/s";
      }

      // First we get either permeability value or the file name
      int file(0);
      std::string file_name, model, format;
      std::vector<std::string> attr_names;
      double kx, ky, kz;

      kx = GetAttributeValueD_(node, "x", TYPE_NUMERICAL, 0.0, DVAL_MAX, unit, false, -1.0);
      ky = GetAttributeValueD_(node, "y", TYPE_NUMERICAL, 0.0, DVAL_MAX, unit, false, -1.0);
      kz = GetAttributeValueD_(node, "z", TYPE_NUMERICAL, 0.0, DVAL_MAX, unit, false, -1.0);

      model = GetAttributeValueS_(node, "model", TYPE_NONE, false, "");
      if (model == "file") file++;
      // format = GetAttributeValueS_(node, "format", TYPE_NONE, false, "");
      // if (format !="") file++;
      file_name = GetAttributeValueS_(node, "filename", TYPE_NONE, false, "");
      if (file_name != "") file++;
      attr_names = GetAttributeVectorS_(node, "attribute", false);
      if (attr_names.size() > 0) file++;

      if (conductivity) {
        kx *= viscosity / (rho_ * const_gravity_);
        ky *= viscosity / (rho_ * const_gravity_);
        kz *= viscosity / (rho_ * const_gravity_);
      }

      // Second, we copy collected data to XML file.
      // For now permeability is not dumped into checkpoint files.
      Teuchos::ParameterList& permeability_ic = out_ic.sublist("permeability");
      permeability_ic.set<bool>("write checkpoint", false);

      if (file == 2) {
        permeability_ic.set<std::string>("restart file", file_name);
        kx = ky = kz = 1.0;
      } else if (file == 3) {
        permeability_ic.sublist("exodus file initialization")
          .set<std::string>("file", file_name)
          .set<Teuchos::Array<std::string>>("attributes", attr_names);
        kx = ky = kz = 1.0;
      } else if (file == 0) {
        if (ky < 0) ky = kz; // x-z system was defined
        Teuchos::ParameterList& aux_list = permeability_ic.sublist("function")
                                             .sublist(reg_str)
                                             .set<Teuchos::Array<std::string>>("regions", regions)
                                             .set<std::string>("component", "cell")
                                             .sublist("function");
        aux_list.set<int>("number of dofs", dim_)
          .set<std::string>("function type", "composite function");
        aux_list.sublist("dof 1 function").sublist("function-constant").set<double>("value", kx);
        aux_list.sublist("dof 2 function").sublist("function-constant").set<double>("value", ky);
        if (dim_ == 3) {
          aux_list.sublist("dof 3 function").sublist("function-constant").set<double>("value", kz);
        } else {
          kz = 0.0;
        }
      } else {
        perm_err = true;
      }
      if (kx < 0.0 || ky < 0.0 || kz < 0.0 || perm_err) {
        ThrowErrorIllformed_(
          "materials", "permeability/hydraulic_conductivity", "file/filename/x/y/z");
      }

      // -- bulk modulus
      node = GetUniqueElementByTagsString_(inode, "mechanical_properties, bulk_modulus", flag);
      if (flag) {
        TranslateFieldEvaluator_(node, "bulk_modulus", "Pa", reg_str, regions, out_ic, out_ev);
      }

      // -- rock heat capacity
      node = GetUniqueElementByTagsString_(inode, "thermal_properties, rock_heat_capacity", flag);
      if (flag) {
        double cv =
          GetAttributeValueD_(node, "cv", TYPE_NUMERICAL, DVAL_MIN, DVAL_MAX, "m^2/s^2/K");
        std::string model = GetAttributeValueS_(node, "model", "linear");

        Teuchos::ParameterList& field_ev = out_ev.sublist("internal_energy_rock");
        field_ev.set<std::string>("evaluator type", "iem")
          .set<std::string>("internal energy key", "internal_energy_rock");

        field_ev.sublist("IEM parameters")
          .sublist(reg_str)
          .set<Teuchos::Array<std::string>>("regions", regions)
          .sublist("IEM parameters")
          .set<std::string>("iem type", model)
          .set<double>("heat capacity", cv);
      }
    }
  }

  // optional eos fields
  if (eos_model_ != "") {
    if (phases_[LIQUID].active) {
      AddSecondaryFieldEvaluator_(out_ev,
                                  Keys::getKey("domain", "molar_density_liquid"),
                                  "molar density key",
                                  "eos",
                                  "density");
    }
    if (phases_[GAS].active && phases_[GAS].dissolved.size() > 0) {
      AddSecondaryFieldEvaluator_(
        out_ev, Keys::getKey("domain", "molar_density_gas"), "molar density key", "eos", "density");
    }
  }

  // optional secondary continuum
  node_list = doc_->getElementsByTagName(mm.transcode("materials_secondary_continuum"));
  if (node_list->getLength() > 0) {
    children = node_list->item(0)->getChildNodes();

    for (int i = 0; i < children->getLength(); i++) {
      DOMNode* inode = children->item(i);
      if (DOMNode::ELEMENT_NODE != inode->getNodeType()) continue;

      node = GetUniqueElementByTagsString_(inode, "assigned_regions", flag);
      std::vector<std::string> regions = CharToStrings_(mm.transcode(node->getTextContent()));
      std::string reg_str = CreateNameFromVector_(regions);

      // -- dual porosity: matrix porosity
      node = GetUniqueElementByTagsString_(inode, "mechanical_properties, porosity", flag);
      if (flag) { TranslateFieldIC_(node, "porosity_msp", "-", reg_str, regions, out_ic); }
    }
  }

  // ----------------------------------------------------------------
  // optional fracture network
  // ----------------------------------------------------------------
  if (fracture_regions_.size() > 0) { TranslateCommonContinuumFields_("fracture", out_ic, out_ev); }

  if (fracture_regions_.size() > 0 && eos_model_ == "") {
    AddIndependentFieldEvaluator_(
      out_ev, "fracture-mass_density_liquid", "FRACTURE_NETWORK_INTERNAL", "*", rho_);
    AddIndependentFieldEvaluator_(out_ev,
                                  "fracture-molar_density_liquid",
                                  "FRACTURE_NETWORK_INTERNAL",
                                  "*",
                                  rho_ / 0.0180153333333);
    AddIndependentFieldEvaluator_(
      out_ev, "fracture-viscosity_liquid", "FRACTURE_NETWORK_INTERNAL", "*", viscosity);
  }

  node = GetUniqueElementByTagsString_("fracture_network, materials", flag);
  if (flag) {
    children = node->getChildNodes();

    for (int i = 0; i < children->getLength(); i++) {
      DOMNode* inode = children->item(i);
      if (DOMNode::ELEMENT_NODE != inode->getNodeType()) continue;

      node = GetUniqueElementByTagsString_(inode, "assigned_regions", flag);
      std::vector<std::string> regions = CharToStrings_(mm.transcode(node->getTextContent()));
      std::string reg_str = CreateNameFromVector_(regions);

      // material properties
      // -- aperture
      node = GetUniqueElementByTagsString_(inode, "aperture", flag);
      if (flag) {
        TranslateFieldEvaluator_(
          node, "fracture-aperture", "m", reg_str, regions, out_ic, out_ev, "value", "fracture");
      } else {
        msg << "Element \"aperture\" must be specified for all materials.";
        Exceptions::amanzi_throw(msg);
      }

      // -- diffusion to matrix (Darcy law)
      node = GetUniqueElementByTagsString_(inode, "flow_diffusion_to_matrix", flag);
      if (flag) {
        std::string model = GetAttributeValueS_(node, "model", TYPE_NONE, false, "");

        if (model == "constant") {
          TranslateFieldIC_(
            node, "fracture-diffusion_to_matrix", "s^-1", reg_str, regions, out_ic, "value");
        }
      } else {
        msg << "Element \"flow_diffusion_to_matrix\" must be specified for all materials.";
        Exceptions::amanzi_throw(msg);
      }

      // -- solute diffusion to matrix (Fick law)
      node = GetUniqueElementByTagsString_(inode, "solute_diffusion_to_matrix", flag);
      if (flag) {
        Teuchos::ParameterList& field_ev = out_ev.sublist("fracture-solute_diffusion_to_matrix");

        std::string model = GetAttributeValueS_(node, "model", TYPE_NONE, false, "");

        if (model == "constant") {
          double val = GetAttributeValueD_(node, "value", TYPE_NUMERICAL, 0.0, DVAL_MAX, "m/s");

          auto& tmp = field_ev.set<std::string>("evaluator type", "independent variable")
                        .sublist("function")
                        .sublist(reg_str)
                        .set<Teuchos::Array<std::string>>("regions", regions)
                        .set<std::string>("component", "cell")
                        .sublist("function");

          tmp.set<int>("number of dofs", 2).set<std::string>("function type", "composite function");
          tmp.sublist("dof 1 function").sublist("function-constant").set<double>("value", val);
          tmp.sublist("dof 2 function").sublist("function-constant").set<double>("value", val);
        }
        if (model == "standard") {
          field_ev.set<std::string>("evaluator type", "solute diffusion to matrix");

          field_ev.set<std::string>("porosity key", "porosity")
            .set<std::string>("tortuosity key", "tortuosity")
            .set<std::string>("saturation key", "saturation_liquid")
            .set<double>("molecular diffusion", 0.0);
        }
      }

      // -- heat diffusion to matrix (Fourier law)
      node = GetUniqueElementByTagsString_(inode, "heat_diffusion_to_matrix", flag);
      if (flag) {
        Teuchos::ParameterList& field_ev = out_ev.sublist("fracture-heat_diffusion_to_matrix");

        std::string model = GetAttributeValueS_(node, "model", TYPE_NONE, false, "");

        if (model == "constant") {
          double val = GetAttributeValueD_(node, "value", TYPE_NUMERICAL, 0.0, DVAL_MAX, "W/m^2/K");

          auto& tmp = field_ev.set<std::string>("evaluator type", "independent variable")
                        .sublist("function")
                        .sublist(reg_str)
                        .set<Teuchos::Array<std::string>>("regions", regions)
                        .set<std::string>("component", "cell")
                        .sublist("function");

          tmp.set<int>("number of dofs", 2).set<std::string>("function type", "composite function");
          tmp.sublist("dof 1 function").sublist("function-constant").set<double>("value", val);
          tmp.sublist("dof 2 function").sublist("function-constant").set<double>("value", val);
        }
        if (model == "standard") {
          field_ev.set<std::string>("evaluator type", "heat diffusion to matrix")
            .set<std::string>("thermal conductivity key", "fracture-thermal_conductivity")
            .set<std::string>("aperture key", "fracture-aperture");
        }
      }

      // -- fracture compliance
      node = GetUniqueElementByTagsString_(inode, "mechanical_properties, compliance", flag);
      if (flag) {
        TranslateFieldEvaluator_(node,
                                 "fracture-compliance",
                                 "m*Pa^-1",
                                 reg_str,
                                 regions,
                                 out_ic,
                                 out_ev,
                                 "value",
                                 "fracture");
      }

      // -- thermal conductivity
      node = GetUniqueElementByTagsString_(inode, "heat_flux_to_matrix", flag);
      if (flag) {
        TranslateFieldIC_(
          node, "fracture-heat_diffusion_to_matrix", "", reg_str, regions, out_ic, "normal");
      }

      // -- rock heat capacity
      node = GetUniqueElementByTagsString_(inode, "thermal_properties, rock_heat_capacity", flag);
      if (flag) {
        double cv =
          GetAttributeValueD_(node, "cv", TYPE_NUMERICAL, DVAL_MIN, DVAL_MAX, "m^2/s^2/K");
        std::string model = GetAttributeValueS_(node, "model", "linear");

        Teuchos::ParameterList& field_ev = out_ev.sublist("fracture-internal_energy_rock");
        field_ev.set<std::string>("evaluator type", "iem")
          .set<std::string>("internal energy key", "internal_energy_rock");

        field_ev.sublist("IEM parameters")
          .sublist(reg_str)
          .set<Teuchos::Array<std::string>>("regions", regions)
          .sublist("IEM parameters")
          .set<std::string>("iem type", model)
          .set<double>("heat capacity", cv);
      }
    }

    // -- eos
    if (eos_model_ != "") {
      AddSecondaryFieldEvaluator_(out_ev,
                                  Keys::getKey("fracture", "molar_density_liquid"),
                                  "molar density key",
                                  "eos",
                                  "density");
    }
  }

  // ---------------------------------------------------------
  // initialization of fields via the initial_conditions list.
  // We have to move most fields to evaluaton list
  // ---------------------------------------------------------
  node_list = doc_->getElementsByTagName(mm.transcode("initial_conditions"));
  int nchildren(0);
  if (node_list->getLength() != 0) {
    children = node_list->item(0)->getChildNodes();
    nchildren = children->getLength();
  }

  for (int i = 0; i < nchildren; i++) {
    DOMNode* inode = children->item(i);
    if (DOMNode::ELEMENT_NODE == inode->getNodeType()) {
      node = GetUniqueElementByTagsString_(inode, "assigned_regions", flag);
      text_content = mm.transcode(node->getTextContent());
      std::vector<std::string> regions = CharToStrings_(text_content);
      std::string reg_str = CreateNameFromVector_(regions);

      // ---------------------------------------------------------
      // liquid phase
      // ---------------------------------------------------------
      // -- uniform pressure
      node = GetUniqueElementByTagsString_(
        inode, "liquid_phase, liquid_component, uniform_pressure", flag);
      if (flag) {
        auto ics = ParseCondList_(node, "uniform_pressure", DVAL_MIN, DVAL_MAX, "Pa", false);

        Teuchos::ParameterList& icfn = out_ic.sublist("pressure")
                                         .sublist("function")
                                         .sublist(reg_str)
                                         .set<Teuchos::Array<std::string>>("regions", regions)
                                         .set<std::string>("component", "*")
                                         .sublist("function");

        TranslateGenericMath_(ics, icfn);
      }

      // -- linear pressure
      node = GetUniqueElementByTagsString_(
        inode, "liquid_phase, liquid_component, linear_pressure", flag);
      if (flag) {
        double p = GetAttributeValueD_(node, "value", TYPE_NUMERICAL, DVAL_MIN, DVAL_MAX, "Pa");
        std::vector<double> grad = GetAttributeVectorD_(node, "gradient", dim_, "Pa/m");
        std::vector<double> refc = GetAttributeVectorD_(node, "reference_coord", dim_, "m");

        Teuchos::ParameterList& icfn = out_ic.sublist("pressure")
                                         .sublist("function")
                                         .sublist(reg_str)
                                         .set<Teuchos::Array<std::string>>("regions", regions)
                                         .set<std::string>("component", "cell")
                                         .sublist("function");

        TranslateFunctionGradient_(p, grad, refc, icfn);
      }

      // -- uniform saturation
      node = GetUniqueElementByTagsString_(
        inode, "liquid_phase, liquid_component, uniform_saturation", flag);
      if (flag) {
        auto ics = ParseCondList_(node, "uniform_saturation", 0.0, 1.0, "-", false);

        Teuchos::ParameterList& icfn = out_ic.sublist("saturation_liquid")
                                         .sublist("function")
                                         .sublist(reg_str)
                                         .set<Teuchos::Array<std::string>>("regions", regions)
                                         .set<std::string>("component", "cell")
                                         .sublist("function");

        TranslateGenericMath_(ics, icfn);
      }

      // -- linear saturation
      node = GetUniqueElementByTagsString_(
        inode, "liquid_phase, liquid_component, linear_saturation", flag);
      if (flag) {
        double s = GetAttributeValueD_(node, "value", TYPE_NUMERICAL, 0.0, 1.0, "-");
        std::vector<double> grad = GetAttributeVectorD_(node, "gradient", dim_, "m^-1");
        std::vector<double> refc = GetAttributeVectorD_(node, "reference_coord", dim_, "m");

        Teuchos::ParameterList& icfn = out_ic.sublist("saturation_liquid")
                                         .sublist("function")
                                         .sublist(reg_str)
                                         .set<Teuchos::Array<std::string>>("regions", regions)
                                         .set<std::string>("component", "cell")
                                         .sublist("function");

        TranslateFunctionGradient_(s, grad, refc, icfn);
      }

      // -- darcy_flux, more precisely volumetric flow rate
      node = GetUniqueElementByTagsString_(inode, "liquid_phase, liquid_component, velocity", flag);
      if (flag) {
        std::vector<double> velocity;
        velocity.push_back(GetAttributeValueD_(node, coords_[0].c_str()));
        velocity.push_back(GetAttributeValueD_(node, coords_[1].c_str()));
        if (dim_ == 3) velocity.push_back(GetAttributeValueD_(node, coords_[2].c_str()));

        Teuchos::ParameterList& flowrate_ic = out_ic.sublist("volumetric_flow_rate");
        Teuchos::ParameterList& tmp_list =
          flowrate_ic.set<bool>("dot with normal", true)
            .sublist("function")
            .sublist(reg_str)
            .set<Teuchos::Array<std::string>>("regions", regions)
            .set<std::string>("component", "face")
            .sublist("function")
            .set<int>("number of dofs", dim_)
            .set<std::string>("function type", "composite function");

        for (int k = 0; k != dim_; ++k) {
          std::stringstream dof_str;
          dof_str << "dof " << k + 1 << " function";
          tmp_list.sublist(dof_str.str())
            .sublist("function-constant")
            .set<double>("value", velocity[k]);
        }
      }

      // -- solute concentration or fraction (liquid phase)
      int ncomp_l = phases_[LIQUID].dissolved.size();
      int ncomp_g = phases_[GAS].dissolved.size();
      int ncomp_all = ncomp_l + ncomp_g;

      node = GetUniqueElementByTagsString_(inode, "liquid_phase, solute_component", flag);
      if (flag && ncomp_all > 0) {
        std::vector<double> vals(ncomp_l, 0.0);

        DOMNodeList* children2 = node->getChildNodes();
        for (int j = 0; j < children2->getLength(); ++j) {
          DOMNode* jnode = children2->item(j);
          tagname = mm.transcode(jnode->getNodeName());

          if (strcmp(tagname, "uniform_conc") == 0) {
            std::string unit, text;
            text = GetAttributeValueS_(jnode, "name");
            int m = GetPosition_(phases_[LIQUID].dissolved, text);
            GetAttributeValueD_(
              jnode, "value", TYPE_NUMERICAL, 0.0, DVAL_MAX, "molar"); // just a check
            vals[m] =
              ConvertUnits_(GetAttributeValueS_(jnode, "value"), unit, solute_molar_mass_[text]);
          }
        }

        Teuchos::ParameterList& tcc_ic = out_ic.sublist("total_component_concentration");
        Teuchos::ParameterList& dof_list =
          tcc_ic.sublist("function")
            .sublist(reg_str)
            .set<Teuchos::Array<std::string>>("regions", regions)
            .set<std::string>("component", "cell")
            .sublist("function")
            .set<int>("number of dofs", ncomp_all)
            .set<std::string>("function type", "composite function");

        for (int k = 0; k < ncomp_l; k++) {
          std::string name = phases_[LIQUID].dissolved[k];
          std::stringstream dof_str;
          dof_str << "dof " << k + 1 << " function";
          dof_list.sublist(dof_str.str())
            .sublist("function-constant")
            .set<double>("value", vals[k]);
        }
      }

      // ---------------------------------------------------------
      // gas phase
      // ---------------------------------------------------------
      // -- solute concentation or fraction (gas phase)
      node = GetUniqueElementByTagsString_(inode, "gas_phase, solute_component", flag);
      if (flag) {
        int noffset;
        std::string field_name;
        std::vector<double> vals(ncomp_g, 0.0);

        DOMNodeList* children2 = node->getChildNodes();
        for (int j = 0; j < children2->getLength(); ++j) {
          DOMNode* jnode = children2->item(j);
          tagname = mm.transcode(jnode->getNodeName());

          if (strcmp(tagname, "uniform_conc") == 0) {
            field_name = "total_component_concentration";
            noffset = ncomp_l;
          } else if (strcmp(tagname, "uniform_mole_fraction") == 0) {
            field_name = "mole_fraction_gas";
            noffset = 0;
          } else {
            continue;
          }

          std::string text = GetAttributeValueS_(jnode, "name");
          int m = GetPosition_(phases_[GAS].dissolved, text);
          vals[m] = GetAttributeValueD_(jnode, "value", TYPE_NUMERICAL, 0.0, DVAL_MAX);
        }

        Teuchos::ParameterList& tcc_ic = out_ic.sublist(field_name);
        Teuchos::ParameterList& dof_list =
          tcc_ic.sublist("function")
            .sublist(reg_str)
            .set<Teuchos::Array<std::string>>("regions", regions)
            .set<std::string>("component", "cell")
            .sublist("function")
            .set<int>("number of dofs", noffset + ncomp_g)
            .set<std::string>("function type", "composite function");

        for (int k = 0; k < ncomp_g; k++) {
          std::string name = phases_[GAS].dissolved[k];
          std::stringstream dof_str;
          dof_str << "dof " << noffset + k + 1 << " function";
          dof_list.sublist(dof_str.str())
            .sublist("function-constant")
            .set<double>("value", vals[k]);
        }
      }

      // -- uniform temperature
      node = GetUniqueElementByTagsString_(inode, "uniform_temperature", flag);
      if (flag) {
        auto ics = ParseCondList_(node, "uniform_temperature", 0.0, 1000.0, "K", false);

        Teuchos::ParameterList& icfn = out_ic.sublist("temperature")
                                         .sublist("function")
                                         .sublist(reg_str)
                                         .set<Teuchos::Array<std::string>>("regions", regions)
                                         .set<std::string>("component", "*")
                                         .sublist("function");

        TranslateGenericMath_(ics, icfn);
      }

      // -- temperature
      node = GetUniqueElementByTagsString_(inode, "temperature", flag);
      if (flag) {
        TranslateFieldIC_(node, "temperature", "K", reg_str, regions, out_ic, "velue", { "*" });
      }

      // -- geochemical condition
      node = GetUniqueElementByTagsString_(
        inode, "liquid_phase, geochemistry_component, constraint", flag);
      if (flag) {
        std::string name = GetAttributeValueS_(node, "name");

        out_ic.sublist("geochemical conditions")
          .sublist(name)
          .set<Teuchos::Array<std::string>>("regions", regions);

        TranslateStateICsAmanziGeochemistry_(out_ic, name, regions, "domain");
      }

      // surface fields
      // -- ponded depth
      std::string domain = (dim_ == 2) ? "domain" : "surface";
      node =
        GetUniqueElementByTagsString_(inode, "liquid_phase, liquid_component, ponded_depth", flag);
      if (flag) {
        TranslateFieldIC_(
          node, Keys::getKey(domain, "ponded_depth"), "m", reg_str, regions, out_ic);
      }

      // -- bathymetry
      node = GetUniqueElementByTagsString_("geometric_model, bathymetry", flag);
      if (flag) {
        TranslateFieldIC_(node,
                          Keys::getKey(domain, "bathymetry"),
                          "",
                          reg_str,
                          regions,
                          out_ic,
                          "",
                          { "cell", "node" });
      }
    }
  }

  // initialization of fields via the initial_conditions in optional fracture network
  node = GetUniqueElementByTagsString_("fracture_network, initial_conditions", flag);
  if (flag) {
    children = node->getChildNodes();

    for (int i = 0; i < children->getLength(); i++) {
      DOMNode* inode = children->item(i);
      if (DOMNode::ELEMENT_NODE != inode->getNodeType()) continue;

      node = GetUniqueElementByTagsString_(inode, "assigned_regions", flag);
      std::vector<std::string> regions = CharToStrings_(mm.transcode(node->getTextContent()));
      std::string reg_str = CreateNameFromVector_(regions);

      // -- uniform fracture pressure
      node = GetUniqueElementByTagsString_(
        inode, "liquid_phase, liquid_component, uniform_pressure", flag);
      if (flag) {
        double p = GetAttributeValueD_(node, "value", TYPE_NUMERICAL, DVAL_MIN, DVAL_MAX);

        Teuchos::ParameterList& pressure_ic = out_ic.sublist("fracture-pressure");
        pressure_ic.sublist("function")
          .sublist(reg_str)
          .set<Teuchos::Array<std::string>>("regions", regions)
          .set<std::string>("component", "*")
          .sublist("function")
          .sublist("function-constant")
          .set<double>("value", p);
      }

      // -- linear fracture pressure
      node = GetUniqueElementByTagsString_(
        inode, "liquid_phase, liquid_component, linear_pressure", flag);
      if (flag) {
        double p = GetAttributeValueD_(node, "value", TYPE_NUMERICAL, DVAL_MIN, DVAL_MAX, "Pa");
        std::vector<double> grad = GetAttributeVectorD_(node, "gradient", dim_, "Pa/m");
        std::vector<double> refc = GetAttributeVectorD_(node, "reference_coord", dim_, "m");

        Teuchos::ParameterList& icfn = out_ic.sublist("fracture-pressure")
                                         .sublist("function")
                                         .sublist(reg_str)
                                         .set<Teuchos::Array<std::string>>("regions", regions)
                                         .set<std::string>("component", "*")
                                         .sublist("function");

        TranslateFunctionGradient_(p, grad, refc, icfn);
      }

      // -- fracture volumetric flow rate
      node = GetUniqueElementByTagsString_(inode, "liquid_phase, liquid_component, velocity", flag);
      if (flag) {
        std::vector<double> velocity;
        velocity.push_back(GetAttributeValueD_(node, coords_[0].c_str()));
        velocity.push_back(GetAttributeValueD_(node, coords_[1].c_str()));
        velocity.push_back(GetAttributeValueD_(node, coords_[2].c_str()));

        Teuchos::ParameterList& flowrate_ic = out_ic.sublist("fracture-volumetric_flow_rate");
        Teuchos::ParameterList& tmp_list =
          flowrate_ic.set<bool>("dot with normal", true)
            .sublist("function")
            .sublist(reg_str)
            .set<Teuchos::Array<std::string>>("regions", regions)
            .set<std::string>("component", "face")
            .sublist("function")
            .set<int>("number of dofs", dim_)
            .set<std::string>("function type", "composite function");

        for (int k = 0; k != dim_; ++k) {
          std::stringstream dof_str;
          dof_str << "dof " << k + 1 << " function";
          tmp_list.sublist(dof_str.str())
            .sublist("function-constant")
            .set<double>("value", velocity[k]);
        }
      }

      // -- total_component_concentration (liquid phase)
      int ncomp_l = phases_[LIQUID].dissolved.size();

      node = GetUniqueElementByTagsString_(inode, "liquid_phase, solute_component", flag);
      if (flag && ncomp_l > 0) {
        std::vector<double> vals(ncomp_l, 0.0);

        DOMNodeList* children2 = node->getChildNodes();
        for (int j = 0; j < children2->getLength(); ++j) {
          DOMNode* jnode = children2->item(j);
          tagname = mm.transcode(jnode->getNodeName());

          if (strcmp(tagname, "uniform_conc") == 0) {
            std::string unit, text;
            text = GetAttributeValueS_(jnode, "name");
            int m = GetPosition_(phases_[LIQUID].dissolved, text);
            GetAttributeValueD_(
              jnode, "value", TYPE_NUMERICAL, 0.0, DVAL_MAX, "molar"); // just a check
            vals[m] =
              ConvertUnits_(GetAttributeValueS_(jnode, "value"), unit, solute_molar_mass_[text]);
          }
        }

        Teuchos::ParameterList& tcc_ic = out_ic.sublist("fracture-total_component_concentration");
        Teuchos::ParameterList& dof_list =
          tcc_ic.sublist("function")
            .sublist(reg_str)
            .set<Teuchos::Array<std::string>>("regions", regions)
            .set<std::string>("component", "cell")
            .sublist("function")
            .set<int>("number of dofs", ncomp_l)
            .set<std::string>("function type", "composite function");

        for (int k = 0; k < ncomp_l; k++) {
          std::string name = phases_[LIQUID].dissolved[k];
          std::stringstream dof_str;
          dof_str << "dof " << k + 1 << " function";
          dof_list.sublist(dof_str.str())
            .sublist("function-constant")
            .set<double>("value", vals[k]);
        }
      }

      // -- geochemical condition
      node = GetUniqueElementByTagsString_(
        inode, "liquid_phase, geochemistry_component, constraint", flag);
      if (flag) {
        std::string name = GetAttributeValueS_(node, "name");

        out_ic.sublist("geochemical conditions")
          .sublist(name)
          .set<Teuchos::Array<std::string>>("regions", regions);

        TranslateStateICsAmanziGeochemistry_(out_ic, name, regions, "fracture");
      }

      // -- uniform temperature
      node = GetUniqueElementByTagsString_(inode, "uniform_temperature", flag);
      if (flag) {
        double val = GetAttributeValueD_(node, "value", TYPE_NUMERICAL, 0.0, 1000.0, "K");

        Teuchos::ParameterList& temperature_ic = out_ic.sublist("fracture-temperature");
        temperature_ic.sublist("function")
          .sublist(reg_str)
          .set<Teuchos::Array<std::string>>("regions", regions)
          .set<std::string>("component", "*")
          .sublist("function")
          .sublist("function-constant")
          .set<double>("value", val);
      }

      // -- temperature
      node = GetUniqueElementByTagsString_(inode, "temperature", flag);
      if (flag)
        TranslateFieldIC_(
          node, "fracture-temperature", "K", reg_str, regions, out_ic, "value", { "*" });
    }
  }

  // atmospheric pressure
  out_ic.sublist("atmospheric_pressure").set<double>("value", const_atm_pressure_);

  // add mesh partitions to the state list
  out_list.sublist("mesh partitions") = TranslateMaterialsPartition_();

  // visualization blacklist and whitelist
  node = GetUniqueElementByTagsString_("output, vis, blacklist", flag);
  if (flag) {
    text_content = mm.transcode(node->getTextContent());
    out_list.set<Teuchos::Array<std::string>>("blacklist", CharToStrings_(text_content));
  }

  node = GetUniqueElementByTagsString_("output, vis, whitelist", flag);
  if (flag) {
    text_content = mm.transcode(node->getTextContent());
    out_list.set<Teuchos::Array<std::string>>("whitelist", CharToStrings_(text_content));
  }

  out_list.sublist("verbose object") = verb_list_.sublist("verbose object");

  return out_list;
}


/* ******************************************************************
* Select proper evaluator based on the list of input parameters.
****************************************************************** */
void
InputConverterU::TranslateCommonContinuumFields_(const std::string& domain,
                                                 Teuchos::ParameterList& out_ic,
                                                 Teuchos::ParameterList& out_ev)
{
  MemoryManager mm;
  Errors::Message msg;

  bool flag;
  DOMNode* node;
  DOMNodeList* children;

  // material independent fields
  // -- viscosity
  node = GetUniqueElementByTagsString_("phases, gas_phase, viscosity", flag);
  if (flag) {
    double viscosity = GetTextContentD_(node, "Pa*s");
    out_ic.sublist("const_gas_viscosity").set<double>("value", viscosity);
    AddIndependentFieldEvaluator_(
      out_ev, Keys::getKey(domain, "viscosity_gas"), "All", "cell", viscosity);
  }

  // -- molar heat capacity
  node = GetUniqueElementByTagsString_("phases, gas_phase, molar_heat_capacity", flag);
  if (flag) {
    double cv = GetTextContentD_(node, "J/mol/K");

    auto key = Keys::getKey(domain, "internal_energy_gas");
    Teuchos::ParameterList& field_ev = out_ev.sublist(key);
    field_ev.set<std::string>("evaluator type", "iem").set<std::string>("internal energy key", key);

    field_ev.sublist("IEM parameters")
      .sublist("All")
      .set<Teuchos::Array<std::string>>("regions", std::vector<std::string>({ "All" }))
      .sublist("IEM parameters")
      .set<std::string>("iem type", "linear")
      .set<double>("heat capacity", cv);
  }

  if (domain == "domain") {
    DOMNodeList* node_list = doc_->getElementsByTagName(mm.transcode("materials"));
    children = node_list->item(0)->getChildNodes();
  } else {
    node = GetUniqueElementByTagsString_("fracture_network, materials", flag);
    children = node->getChildNodes();
  }

  for (int i = 0; i < children->getLength(); i++) {
    DOMNode* inode = children->item(i);
    if (DOMNode::ELEMENT_NODE == inode->getNodeType()) {
      node = GetUniqueElementByTagsString_(inode, "assigned_regions", flag);
      std::vector<std::string> regions = CharToStrings_(mm.transcode(node->getTextContent()));
      std::string reg_str = CreateNameFromVector_(regions);

      // porosity
      node = GetUniqueElementByTagsString_(inode, "mechanical_properties, porosity", flag);
      if (flag) {
        TranslateFieldEvaluator_(
          node, Keys::getKey(domain, "porosity"), "-", reg_str, regions, out_ic, out_ev);
      } else {
        msg << "Porosity element must be specified under mechanical_properties";
        Exceptions::amanzi_throw(msg);
      }

      // specific_yield
      node = GetUniqueElementByTagsString_(inode, "mechanical_properties, specific_yield", flag);
      if (flag)
        TranslateFieldIC_(
          node, Keys::getKey(domain, "specific_yield"), "-", reg_str, regions, out_ic);

      // specific storage
      node = GetUniqueElementByTagsString_(inode, "mechanical_properties, specific_storage", flag);
      if (flag) {
        auto key = Keys::getKey(domain, "specific_storage");
        Teuchos::ParameterList& field_ev = out_ev.sublist(key);

        std::string model = GetAttributeValueS_(node, "model", TYPE_NONE, false, "");
        field_ev.set<std::string>("evaluator type", "specific storage");

        Teuchos::ParameterList& params =
          field_ev.sublist("specific storage parameters").sublist(reg_str);
        params.set<Teuchos::Array<std::string>>("regions", regions);

        double val1, val2;
        if (model == "constant") {
          val1 = GetAttributeValueD_(node, "value", TYPE_NUMERICAL, 0.0, DVAL_MAX, "m^-1");
          params.set<std::string>("model", "constant").set<double>("value", val1);
        }
        if (model == "standard") {
          std::vector<std::string> deps({ Keys::getKey(domain, "porosity") });
          field_ev.set<Teuchos::Array<std::string>>("dependencies", deps);

          val1 = GetAttributeValueD_(
            node, "fluid_compressibility", TYPE_NUMERICAL, 0.0, DVAL_MAX, "Pa^-1");
          val2 = GetAttributeValueD_(
            node, "matrix_compressibility", TYPE_NUMERICAL, 0.0, DVAL_MAX, "Pa^-1");
          params.set<std::string>("model", "standard")
            .set<double>("fluid compressibility", val1)
            .set<double>("matrix compressibility", val2)
            .set<double>("gravity", const_gravity_)
            .set<double>("fluid density", rho_);
        }
      }

      // particle density
      node = GetUniqueElementByTagsString_(inode, "mechanical_properties, particle_density", flag);
      if (flag)
        TranslateFieldEvaluator_(node,
                                 Keys::getKey(domain, "particle_density"),
                                 "kg*m^-3",
                                 reg_str,
                                 regions,
                                 out_ic,
                                 out_ev,
                                 "value",
                                 domain);

      // internal energy for liquid
      node = GetUniqueElementByTagsString_(inode, "thermal_properties, liquid_heat_capacity", flag);
      if (flag) {
        double cv =
          GetAttributeValueD_(node, "cv", TYPE_NUMERICAL, DVAL_MIN, DVAL_MAX, "kg*m^2/s^2/mol/K");
        std::string model = GetAttributeValueS_(node, "model", "linear");

        auto key = Keys::getKey(domain, "internal_energy_liquid");
        Teuchos::ParameterList& field_ev = out_ev.sublist(key);
        field_ev.set<std::string>("evaluator type", "iem")
          .set<std::string>("internal energy key", key);

        field_ev.sublist("IEM parameters")
          .sublist(reg_str)
          .set<Teuchos::Array<std::string>>("regions", regions)
          .sublist("IEM parameters")
          .set<std::string>("iem type", model)
          .set<double>("heat capacity", cv);

        if (eos_lookup_table_.size() > 0) {
          field_ev.sublist("IEM parameters")
            .sublist(reg_str)
            .set<Teuchos::Array<std::string>>("regions", regions)
            .sublist("IEM parameters")
            .set<std::string>("iem type", "tabular")
            .set<std::string>("table name", eos_lookup_table_)
            .set<std::string>("field name", "internal_energy");
        }
      }
    }
  }
}


/* ******************************************************************
* Select proper evaluator based on the list of input parameters.
****************************************************************** */
void
InputConverterU::TranslateFieldEvaluator_(DOMNode* node,
                                          const std::string& field,
                                          const std::string& unit,
                                          const std::string& reg_str,
                                          const std::vector<std::string>& regions,
                                          Teuchos::ParameterList& out_ic,
                                          Teuchos::ParameterList& out_ev,
                                          std::string data_key,
                                          std::string domain)
{
  MemoryManager mm;

  std::string model = GetAttributeValueS_(node, "model", TYPE_NONE, false, "");
  if (model == "file") { // Amanzi restart file
    std::string filename = GetAttributeValueS_(node, "filename");
    Teuchos::ParameterList& field_ic = out_ic.sublist(field);
    field_ic.set<std::string>("restart file", filename);

    Teuchos::ParameterList& field_ev = out_ev.sublist(field);
    field_ev.set<std::string>("evaluator type", "constant variable");
  } else if (model == "h5file") { // regular h5 file
    std::string filename = GetAttributeValueS_(node, "filename");
    bool temporal =
      GetAttributeValueS_(node, "constant_in_time", TYPE_NUMERICAL, false, "true") == "true";

    Teuchos::ParameterList& field_ev = out_ev.sublist(field);
    field_ev.set<std::string>("evaluator type", "independent variable from file")
      .set<std::string>("filename", filename)
      .set<std::string>("domain name", domain)
      .set<std::string>("component name", "cell")
      .set<std::string>("mesh entity", "cell")
      .set<std::string>("variable name", field)
      .set<int>("number of dofs", 1)
      .set<bool>("constant in time", temporal);
  } else if (model == "") {
    Teuchos::ParameterList& field_ev = out_ev.sublist(field);

    if (static_cast<DOMElement*>(node)->hasAttribute(mm.transcode("formula"))) {
      std::string formula = GetAttributeValueS_(node, "formula");

      field_ev.sublist("function")
        .sublist(reg_str)
        .set<Teuchos::Array<std::string>>("regions", regions)
        .set<std::string>("component", "cell")
        .sublist("function")
        .sublist("function-exprtk")
        .set<int>("number of arguments", dim_ + 1)
        .set<std::string>("formula", formula);
      field_ev.set<std::string>("evaluator type", "independent variable");

    } else {
      double val =
        GetAttributeValueD_(node, data_key.c_str(), TYPE_NUMERICAL, DVAL_MIN, DVAL_MAX, unit);

      field_ev.sublist("function")
        .sublist(reg_str)
        .set<Teuchos::Array<std::string>>("regions", regions)
        .set<std::string>("component", "cell")
        .sublist("function")
        .sublist("function-constant")
        .set<double>("value", val);
      field_ev.set<std::string>("evaluator type", "independent variable")
        .set<bool>("constant in time", true);
    }
  }
}


/* ******************************************************************
* Select proper IC based on the list of input parameters.
****************************************************************** */
void
InputConverterU::TranslateFieldIC_(DOMNode* node,
                                   std::string field,
                                   std::string unit,
                                   const std::string& reg_str,
                                   const std::vector<std::string>& regions,
                                   Teuchos::ParameterList& out_ic,
                                   std::string data_key,
                                   const std::vector<std::string>& components)
{
  MemoryManager mm;

  std::string type = GetAttributeValueS_(node, "type", TYPE_NONE, false, "");
  if (type == "file") {
    std::string filename = GetAttributeValueS_(node, "filename");
    Teuchos::ParameterList& field_ic = out_ic.sublist(field);
    field_ic.set<std::string>("restart file", filename);
  } else {
    Teuchos::ParameterList& field_ic = out_ic.sublist(field);

    if (static_cast<DOMElement*>(node)->hasAttribute(mm.transcode("formula"))) {
      std::string formula = GetAttributeValueS_(node, "formula");

      field_ic.sublist("function")
        .sublist(reg_str)
        .set<Teuchos::Array<std::string>>("regions", regions)
        .sublist("function")
        .sublist("function-exprtk")
        .set<int>("number of arguments", dim_ + 1)
        .set<std::string>("formula", formula);

      if (components.size() == 1) {
        field_ic.sublist("function")
          .sublist(reg_str)
          .set<Teuchos::Array<std::string>>("regions", regions)
          .set<std::string>("component", components[0]);
      } else {
        field_ic.sublist("function")
          .sublist(reg_str)
          .set<Teuchos::Array<std::string>>("regions", regions)
          .set<Teuchos::Array<std::string>>("components", components);
      }

    } else {
      double val =
        GetAttributeValueD_(node, data_key.c_str(), TYPE_NUMERICAL, DVAL_MIN, DVAL_MAX, unit);

      field_ic.sublist("function")
        .sublist(reg_str)
        .set<Teuchos::Array<std::string>>("regions", regions)
        .set<std::string>("component", "cell")
        .sublist("function")
        .sublist("function-constant")
        .set<double>("value", val);
    }
  }
}


/* ******************************************************************
* Mesh patition sublist based on materials
****************************************************************** */
Teuchos::ParameterList
InputConverterU::TranslateMaterialsPartition_()
{
  Teuchos::ParameterList out_list;
  Teuchos::ParameterList& tmp_list = out_list.sublist("materials");

  MemoryManager mm;
  bool flag;
  std::vector<std::string> regions;

  DOMNodeList* node_list = doc_->getElementsByTagName(mm.transcode("materials"));
  DOMNodeList* children = node_list->item(0)->getChildNodes();

  for (int i = 0; i < children->getLength(); i++) {
    DOMNode* inode = children->item(i);
    if (DOMNode::ELEMENT_NODE == inode->getNodeType()) {
      DOMNode* node = GetUniqueElementByTagsString_(inode, "assigned_regions", flag);
      if (flag) {
        char* text_content = mm.transcode(node->getTextContent());
        std::vector<std::string> names = CharToStrings_(text_content);

        for (int n = 0; n < names.size(); ++n) { regions.push_back(names[n]); }
      }
    }
  }
  tmp_list.set<Teuchos::Array<std::string>>("region list", regions);

  return out_list;
}


/* ******************************************************************
* Create initialization list for concentration. This routine is called
* when geochemistry list exists for initial conditions.
****************************************************************** */
void
InputConverterU::TranslateStateICsAmanziGeochemistry_(Teuchos::ParameterList& out_list,
                                                      std::string& constraint,
                                                      std::vector<std::string>& regions,
                                                      const std::string& domain)
{
  if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH) {
    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "Compatibility mode: translating ICs for native chemistry: " << domain
               << std::endl;
  }

  bool flag;
  DOMNode* node;
  DOMElement* element;

  node = GetPKChemistryPointer_(flag);
  std::string engine = GetAttributeValueS_(node, "engine");

  node = GetUniqueElementByTagsString_("geochemistry, constraints", flag);
  if (flag && engine == "amanzi") {
    std::string name;
    element = GetUniqueChildByAttribute_(node, "name", constraint, flag, true);
    std::vector<DOMNode*> children = GetSameChildNodes_(element, name, flag);
    if (children.size() != phases_[LIQUID].dissolved.size()) {
      Errors::Message msg;
      msg << "Constraint \"" << constraint << "\" is not backward compatible: "
          << " check the number of components.";
      Exceptions::amanzi_throw(msg);
    }

    Key tcc_key = Keys::getKey(domain, "total_component_concentration");
    Teuchos::ParameterList& ic_list = out_list.sublist(tcc_key).sublist("function").sublist("All");

    ic_list.set<Teuchos::Array<std::string>>("regions", regions)
      .set<std::string>("component", "cell");

    Teuchos::ParameterList& tmp_list = ic_list.sublist("function")
                                         .set<int>("number of dofs", children.size())
                                         .set<std::string>("function type", "composite function");

    int nsolutes = children.size();
    Teuchos::Array<std::string> types(nsolutes);

    for (int i = 0; i < children.size(); ++i) {
      std::string species = GetAttributeValueS_(children[i], "name");
      double val = GetAttributeValueD_(children[i], "value", TYPE_NUMERICAL, DVAL_MIN, DVAL_MAX);
      types[i] = GetAttributeValueS_(children[i], "type");

      // find position of species in the list of component names
      int k(-1);
      for (int n = 0; n < comp_names_all_.size(); ++n) {
        if (comp_names_all_[n] == species) {
          k = n;
          break;
        }
      }

      std::stringstream dof_str;
      dof_str << "dof " << k + 1 << " function";
      tmp_list.sublist(dof_str.str()).sublist("function-constant").set<double>("value", val);
    }

    out_list.sublist(tcc_key).set<Teuchos::Array<std::string>>("names", types);
  }
}


/* ******************************************************************
* Add independent field evaluator
****************************************************************** */
void
InputConverterU::AddIndependentFieldEvaluator_(Teuchos::ParameterList& out_ev,
                                               const std::string& field,
                                               const std::string& region,
                                               const std::string& comp,
                                               double val)
{
  out_ev.sublist(field)
    .sublist("function")
    .sublist("All")
    .set<std::string>("region", region)
    .set<std::string>("component", comp)
    .sublist("function")
    .sublist("function-constant")
    .set<double>("value", val);

  out_ev.sublist(field).set<std::string>("evaluator type", "independent variable")
    .set<bool>("constant in time", true);
}


/* ******************************************************************
* Add secondary field evaluator
****************************************************************** */
void
InputConverterU::AddSecondaryFieldEvaluator_(Teuchos::ParameterList& out_ev,
                                             const Key& field,
                                             const Key& key,
                                             const std::string& type,
                                             const std::string& eos_table_name)
{
  out_ev.sublist(field).set<std::string>("evaluator type", type).set<std::string>(key, field);

  out_ev.sublist(field)
    .sublist("EOS parameters")
    .set<std::string>("eos type", "liquid water " + eos_model_);

  // modifies
  if (eos_lookup_table_ != "") {
    out_ev.sublist(field)
      .sublist("EOS parameters")
      .set<std::string>("eos type", "lookup table")
      .set<std::string>("format", "Amanzi")
      .set<std::string>("table name", eos_lookup_table_)
      .set<std::string>("field name", eos_table_name);
  }

  // extensions
  Key prefix = Keys::split(field, '-').first;
  std::string basename = Keys::split(field, '-').second;

  if (basename == "molar_density_liquid") {
    out_ev.sublist(field)
      .set<std::string>("eos basis", "both")
      .set<std::string>("mass density key", Keys::getKey(prefix, "mass_density_liquid"));
  }

  if (basename == "molar_density_gas") {
    out_ev.sublist(field)
      .set<std::string>("eos basis", "molar")
      .set<std::string>("mass density key", Keys::getKey(prefix, "mass_density_gas"));
  }
}


/* ******************************************************************
* Add constant field
****************************************************************** */
void
InputConverterU::AddConstantFieldInitialization_(Teuchos::ParameterList& out_ic,
                                                 const std::string& field,
                                                 const std::string& region,
                                                 double val)
{
  out_ic.sublist(field)
    .sublist("function")
    .sublist("All")
    .set<std::string>("region", region)
    .set<std::string>("component", "cell")
    .sublist("function")
    .sublist("function-constant")
    .set<double>("value", val);
}

} // namespace AmanziInput
} // namespace Amanzi
