/*
  Input Converter 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Markus Berndt (original version)
           Konstantin Lipnikov (lipnikov@lanl.gov)
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

#include "InputConverterU.hh"
#include "InputConverterU_Defs.hh"

namespace Amanzi {
namespace AmanziInput {

XERCES_CPP_NAMESPACE_USE

/* ******************************************************************
* STATE sublist
****************************************************************** */
Teuchos::ParameterList InputConverterU::TranslateState_()
{
  Teuchos::ParameterList out_list;

  if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH) {
    *vo_->os() << "Translating state" << std::endl;
  }

  // first we write initial conditions for scalars and vectors, not region-specific
  Teuchos::ParameterList& out_ev = out_list.sublist("field evaluators");
  Teuchos::ParameterList& out_ic = out_list.sublist("initial conditions");

  MemoryManager mm;
  Errors::Message msg;
  char* tagname;
  char* text_content;
  
  // --- gravity
  Teuchos::Array<double> gravity(dim_);
  for (int i = 0; i != dim_-1; ++i) gravity[i] = 0.0;
  gravity[dim_-1] = -const_gravity_;
  out_ic.sublist("gravity").set<Teuchos::Array<double> >("value", gravity);

  // --- viscosity
  bool flag;
  DOMNode* node = GetUniqueElementByTagsString_("phases, liquid_phase, viscosity", flag);
  double viscosity = GetTextContentD_(node, "Pa*s");
  out_ic.sublist("fluid_viscosity").set<double>("value", viscosity);

  // --- constant density
  node = GetUniqueElementByTagsString_("phases, liquid_phase, density", flag);
  rho_ = GetTextContentD_(node, "kg/m^3");
  out_ic.sublist("fluid_density").set<double>("value", rho_);

  out_ev.sublist("mass_density_liquid").sublist("function").sublist("All")
      .set<std::string>("region", "All")
      .set<std::string>("component", "cell")
      .sublist("function").sublist("function-constant")
      .set<double>("value", rho_);
  out_ev.sublist("mass_density_liquid")
        .set<std::string>("field evaluator type", "independent variable");

  // --- region specific initial conditions from material properties
  std::map<std::string, int> reg2mat;
  int mat(0);

  // primary continuum
  int nmat(0);
  DOMNodeList* node_list = doc_->getElementsByTagName(mm.transcode("materials"));
  DOMNodeList* children = node_list->item(0)->getChildNodes();

  for (int i = 0; i < children->getLength(); i++) {
    DOMNode* inode = children->item(i);
    if (DOMNode::ELEMENT_NODE == inode->getNodeType()) {
      std::string mat_name = GetAttributeValueS_(inode, "name");
      int mat_id = GetAttributeValueL_(inode, "id", TYPE_NUMERICAL, 0, INT_MAX, false, -1);
      nmat++;

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

      // -- porosity: skip if compressibility model was already provided.
      if (!compressibility_) {
        node = GetUniqueElementByTagsString_(inode, "mechanical_properties, porosity", flag);
        if (flag) {
          TranslateFieldEvaluator_(node, "porosity", "-", reg_str, regions, out_ic, out_ev);
        } else {
          msg << "Porosity element must be specified under mechanical_properties";
          Exceptions::amanzi_throw(msg);
        }
      }

      // -- transport porosity 
      node = GetUniqueElementByTagsString_(inode, "mechanical_properties, transport_porosity", flag);
      if (flag) {
        use_transport_porosity_ = true;
        TranslateFieldEvaluator_(node, "transport_porosity", "-", reg_str, regions, out_ic, out_ev);
      } else if (use_transport_porosity_) {
        msg << "Transport porosity element must be specified for all materials or none.";
        Exceptions::amanzi_throw(msg);
      }

      // -- permeability. We parse matrix and fractures together.
      bool perm_err(false), perm_init_from_file(false), conductivity(false);
      std::string perm_file, perm_format, unit("m^2");

      node = GetUniqueElementByTagsString_(inode, "permeability", flag);
      if (!flag) {
        conductivity = true;
        node = GetUniqueElementByTagsString_(inode, "hydraulic_conductivity", flag);
        unit = "m/s"; 
      }

      // First we get either permeability value or the file name
      int file(0);
      std::string file_name, type, format;
      std::vector<std::string> attr_names;
      double kx, ky, kz;

      kx = GetAttributeValueD_(node, "x", TYPE_NUMERICAL, 0.0, DVAL_MAX, unit, false, -1.0);
      ky = GetAttributeValueD_(node, "y", TYPE_NUMERICAL, 0.0, DVAL_MAX, unit, false, -1.0);
      kz = GetAttributeValueD_(node, "z", TYPE_NUMERICAL, 0.0, DVAL_MAX, unit, false, -1.0);

      type = GetAttributeValueS_(node, "type", TYPE_NONE, false, "");
      if (type == "file") file++;
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
      }
      else if (file == 3) {
        permeability_ic.sublist("exodus file initialization")
            .set<std::string>("file", file_name)
            .set<Teuchos::Array<std::string> >("attributes", attr_names);
        kx = ky = kz = 1.0;
      } else if (file == 0) {
        if (ky < 0) ky = kz;  // x-z system was defined
        Teuchos::ParameterList& aux_list = permeability_ic.sublist("function").sublist(reg_str)
            .set<Teuchos::Array<std::string> >("regions", regions)
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
        ThrowErrorIllformed_("materials", "permeability/hydraulic_conductivity", "file/filename/x/y/z");
      }

      // -- specific_yield
      node = GetUniqueElementByTagsString_(inode, "mechanical_properties, specific_yield", flag);
      if (flag) {
        TranslateFieldIC_(node, "specific_yield", "-", reg_str, regions, out_ic, out_ev);
      }

      // -- specific storage
      node = GetUniqueElementByTagsString_(inode, "mechanical_properties, specific_storage", flag);
      if (flag) {
        TranslateFieldIC_(node, "specific_storage", "m^-1", reg_str, regions, out_ic, out_ev);
      }

      // -- particle density
      node = GetUniqueElementByTagsString_(inode, "mechanical_properties, particle_density", flag);
      if (flag) {
        TranslateFieldEvaluator_(node, "particle_density", "kg*m^-3", reg_str, regions, out_ic, out_ev);
      }

      // -- liquid heat capacity
      node = GetUniqueElementByTagsString_(inode, "thermal_properties, liquid_heat_capacity", flag);
      if (flag) {
        if (nmat > 1) {
          msg << "Heat capacity is supported for problems with one material";
          Exceptions::amanzi_throw(msg);
        }
        double cv = GetAttributeValueD_(node, "cv", TYPE_NUMERICAL, DVAL_MIN, DVAL_MAX, "kg*m^2/s^2/mol/K");
        std::string model = GetAttributeValueS_(node, "model", "linear");

        Teuchos::ParameterList& field_ev = out_ev.sublist("internal_energy_liquid");
        field_ev.set<std::string>("field evaluator type", "iem")
            .set<std::string>("internal energy key", "internal_energy_liquid");
        field_ev.sublist("IEM parameters")
            .set<std::string>("iem type", model)
            .set<double>("heat capacity", cv);
      }

      // -- rock heat capacity
      node = GetUniqueElementByTagsString_(inode, "thermal_properties, rock_heat_capacity", flag);
      if (flag) {
        if (nmat > 1) {
          msg << "Heat capacity is supported for problems with one material";
          Exceptions::amanzi_throw(msg);
        }
        double cv = GetAttributeValueD_(node, "cv", TYPE_NUMERICAL, DVAL_MIN, DVAL_MAX, "m^2/s^2/K");
        std::string model = GetAttributeValueS_(node, "model", "linear");

        Teuchos::ParameterList& field_ev = out_ev.sublist("internal_energy_rock");
        field_ev.set<std::string>("field evaluator type", "iem")
            .set<std::string>("internal energy key", "internal_energy_rock");
        field_ev.sublist("IEM parameters")
            .set<std::string>("iem type", model)
            .set<double>("heat capacity", cv);
      }
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
      if (flag) {
        TranslateFieldIC_(node, "porosity_matrix", "-", reg_str, regions, out_ic, out_ev);
      }
    }
  }

  // optional fracture network
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
      // -- porosity
      node = GetUniqueElementByTagsString_(inode, "mechanical_properties, porosity", flag);
      if (flag) {
        TranslateFieldEvaluator_(node, "fracture-porosity", "-", reg_str, regions, out_ic, out_ev);
      }

      // -- aperture
      node = GetUniqueElementByTagsString_(inode, "fracture_permeability", flag);
      if (flag) {
        TranslateFieldEvaluator_(node, "fracture-aperture", "m", reg_str, regions, out_ic, out_ev, "aperture", "fracture");
        TranslateFieldIC_(node, "fracture-normal_permeability", "m^2*s/kg", reg_str, regions, out_ic, out_ev, "normal");
      } else { 
        msg << "fracture_permeability element must be specified for all materials in fracture network.";
        Exceptions::amanzi_throw(msg);
      }
    }
  }

  // initialization of fields via the initial_conditions list
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

      // create regions string
      std::string reg_str = CreateNameFromVector_(regions);

      // -- uniform pressure
      node = GetUniqueElementByTagsString_(inode, "liquid_phase, liquid_component, uniform_pressure", flag);
      if (flag) {
        double p = GetAttributeValueD_(node, "value", TYPE_NUMERICAL, DVAL_MIN, DVAL_MAX);

        Teuchos::ParameterList& pressure_ic = out_ic.sublist("pressure");
        pressure_ic.sublist("function").sublist(reg_str)
            .set<Teuchos::Array<std::string> >("regions", regions)
            .set<std::string>("component", "cell")
            .sublist("function").sublist("function-constant")
            .set<double>("value", p);
      }

      // -- linear pressure
      node = GetUniqueElementByTagsString_(inode, "liquid_phase, liquid_component, linear_pressure", flag);
      if (flag) {
        double p = GetAttributeValueD_(node, "value", TYPE_NUMERICAL, DVAL_MIN, DVAL_MAX, "Pa");
        std::vector<double> grad = GetAttributeVectorD_(node, "gradient", dim_, "Pa/m");
        std::vector<double> refc = GetAttributeVectorD_(node, "reference_coord", dim_, "m");

        grad.insert(grad.begin(), 0.0);
        refc.insert(refc.begin(), 0.0);

        Teuchos::ParameterList& pressure_ic = out_ic.sublist("pressure");
        pressure_ic.sublist("function").sublist(reg_str)
            .set<Teuchos::Array<std::string> >("regions", regions)
            .set<std::string>("component", "cell")
            .sublist("function").sublist("function-linear")
            .set<double>("y0", p)
            .set<Teuchos::Array<double> >("x0", refc)
            .set<Teuchos::Array<double> >("gradient", grad);
      }

      // -- uniform saturation
      node = GetUniqueElementByTagsString_(inode, "liquid_phase, liquid_component, uniform_saturation", flag);
      if (flag) {
        double s = GetAttributeValueD_(node, "value", TYPE_NUMERICAL, 0.0, 1.0, "-");

        Teuchos::ParameterList& saturation_ic = out_ic.sublist("saturation_liquid");
        saturation_ic.sublist("function").sublist(reg_str)
            .set<Teuchos::Array<std::string> >("regions", regions)
            .set<std::string>("component", "cell")
            .sublist("function").sublist("function-constant")
            .set<double>("value", s);
      }

      // -- linear saturation
      node = GetUniqueElementByTagsString_(inode, "liquid_phase, liquid_component, linear_saturation", flag);
      if (flag) {
        double s = GetAttributeValueD_(node, "value", TYPE_NUMERICAL, 0.0, 1.0, "-");
        std::vector<double> grad = GetAttributeVectorD_(node, "gradient", dim_, "m^-1");
        std::vector<double> refc = GetAttributeVectorD_(node, "reference_coord", dim_, "m");

        grad.insert(grad.begin(), 0.0);
        refc.insert(refc.begin(), 0.0);

        Teuchos::ParameterList& saturation_ic = out_ic.sublist("saturation_liquid");
        saturation_ic.sublist("function").sublist(reg_str)
            .set<Teuchos::Array<std::string> >("regions", regions)
            .set<std::string>("component", "cell")
            .sublist("function").sublist("function-linear")
            .set<double>("y0", s)
            .set<Teuchos::Array<double> >("x0", refc)
            .set<Teuchos::Array<double> >("gradient", grad);
      }

      // -- darcy_flux
      node = GetUniqueElementByTagsString_(inode, "liquid_phase, liquid_component, velocity", flag);
      if (flag) {
        std::vector<double> velocity;
        velocity.push_back(GetAttributeValueD_(node, coords_[0].c_str()));
        velocity.push_back(GetAttributeValueD_(node, coords_[1].c_str()));
        if (dim_ == 3) velocity.push_back(GetAttributeValueD_(node, coords_[2].c_str()));

        Teuchos::ParameterList& darcy_flux_ic = out_ic.sublist("darcy_flux");
        Teuchos::ParameterList& tmp_list =
            darcy_flux_ic.set<bool>("dot with normal", true)
            .sublist("function").sublist(reg_str)
            .set<Teuchos::Array<std::string> >("regions", regions)
            .set<std::string>("component", "face")
            .sublist("function")
            .set<int>("number of dofs", dim_)
            .set<std::string>("function type", "composite function");

        for (int k = 0; k != dim_; ++k) {
          std::stringstream dof_str;
          dof_str << "dof " << k+1 << " function";
          tmp_list.sublist(dof_str.str()).sublist("function-constant")
                                         .set<double>("value", velocity[k]);
        }
      }

      // -- total_component_concentration (liquid phase)
      int ncomp_l = phases_["water"].size();
      int ncomp_g = phases_["air"].size();
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
            int m = GetPosition_(phases_["water"], text);
            GetAttributeValueD_(jnode, "value", TYPE_NUMERICAL, 0.0, DVAL_MAX, "molar");  // just a check
            vals[m] = ConvertUnits_(GetAttributeValueS_(jnode, "value"), unit, solute_molar_mass_[text]);
          }
        }

        Teuchos::ParameterList& tcc_ic = out_ic.sublist("total_component_concentration");
        Teuchos::ParameterList& dof_list = tcc_ic.sublist("function").sublist(reg_str)
            .set<Teuchos::Array<std::string> >("regions", regions)
            .set<std::string>("component", "cell")
            .sublist("function")
            .set<int>("number of dofs", ncomp_all)
            .set<std::string>("function type", "composite function");

        for (int k = 0; k < ncomp_l; k++) {
          std::string name = phases_["water"][k];
          std::stringstream dof_str;
          dof_str << "dof " << k + 1 << " function";
          dof_list.sublist(dof_str.str()).sublist("function-constant").set<double>("value", vals[k]);
        }
      }

      // -- total_component_concentration (gas phase)
      node = GetUniqueElementByTagsString_(inode, "gas_phase, solute_component", flag);
      if (flag) {
        std::vector<double> vals(ncomp_g, 0.0);

        DOMNodeList* children2 = node->getChildNodes();
        for (int j = 0; j < children2->getLength(); ++j) {
          DOMNode* jnode = children2->item(j);
          tagname = mm.transcode(jnode->getNodeName());

          if (strcmp(tagname, "uniform_conc") == 0) {
            std::string text = GetAttributeValueS_(jnode, "name");
            int m = GetPosition_(phases_["air"], text);
            vals[m] = GetAttributeValueD_(jnode, "value", TYPE_NUMERICAL, 0.0, DVAL_MAX);
          }
        }

        Teuchos::ParameterList& tcc_ic = out_ic.sublist("total_component_concentration");
        Teuchos::ParameterList& dof_list = tcc_ic.sublist("function").sublist(reg_str).sublist("function");
        for (int k = 0; k < ncomp_g; k++) {
          std::string name = phases_["air"][k];
          std::stringstream dof_str;
          dof_str << "dof " << ncomp_l + k + 1 << " function";
          dof_list.sublist(dof_str.str()).sublist("function-constant").set<double>("value", vals[k]);
        }
      }

      // -- uniform temperature
      node = GetUniqueElementByTagsString_(inode, "uniform_temperature", flag);
      if (flag) {
        double val = GetAttributeValueD_(node, "value", TYPE_NUMERICAL, 0.0, 1000.0, "K");

        Teuchos::ParameterList& temperature_ic = out_ic.sublist("temperature");
        temperature_ic.sublist("function").sublist(reg_str)
            .set<Teuchos::Array<std::string> >("regions", regions)
            .set<std::string>("component", "cell")
            .sublist("function").sublist("function-constant")
            .set<double>("value", val);
      }

      // -- geochemical condition
      node = GetUniqueElementByTagsString_(inode, "liquid_phase, geochemistry_component, constraint", flag);
      if (flag) {
        std::string name = GetAttributeValueS_(node, "name");

        out_ic.sublist("geochemical conditions").sublist(name)
            .set<Teuchos::Array<std::string> >("regions", regions);

        TranslateStateICsAmanziGeochemistry_(out_ic, name, regions);
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
      node = GetUniqueElementByTagsString_(inode, "liquid_phase, liquid_component, uniform_pressure", flag);
      if (flag) {
        double p = GetAttributeValueD_(node, "value", TYPE_NUMERICAL, DVAL_MIN, DVAL_MAX);

        Teuchos::ParameterList& pressure_ic = out_ic.sublist("fracture-pressure");
        pressure_ic.sublist("function").sublist(reg_str)
            .set<Teuchos::Array<std::string> >("regions", regions)
            .set<std::string>("component", "cell")
            .sublist("function").sublist("function-constant")
            .set<double>("value", p);
      }

      // -- total_component_concentration (liquid phase)
      int ncomp_l = phases_["water"].size();

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
            int m = GetPosition_(phases_["water"], text);
            GetAttributeValueD_(jnode, "value", TYPE_NUMERICAL, 0.0, DVAL_MAX, "molar");  // just a check
            vals[m] = ConvertUnits_(GetAttributeValueS_(jnode, "value"), unit, solute_molar_mass_[text]);
          }
        }

        Teuchos::ParameterList& tcc_ic = out_ic.sublist("fracture-total_component_concentration");
        Teuchos::ParameterList& dof_list = tcc_ic.sublist("function").sublist(reg_str)
            .set<Teuchos::Array<std::string> >("regions", regions)
            .set<std::string>("component", "cell")
            .sublist("function")
            .set<int>("number of dofs", ncomp_l)
            .set<std::string>("function type", "composite function");

        for (int k = 0; k < ncomp_l; k++) {
          std::string name = phases_["water"][k];
          std::stringstream dof_str;
          dof_str << "dof " << k + 1 << " function";
          dof_list.sublist(dof_str.str()).sublist("function-constant").set<double>("value", vals[k]);
        }
      }
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
    out_list.set<Teuchos::Array<std::string> >("blacklist", CharToStrings_(text_content));
  }

  node = GetUniqueElementByTagsString_("output, vis, whitelist", flag);
  if (flag) {
    text_content = mm.transcode(node->getTextContent());
    out_list.set<Teuchos::Array<std::string> >("whitelist", CharToStrings_(text_content));
  }

  // temporary fix for fractures
  if (fractures_) {
    Teuchos::ParameterList empty;
    out_list.sublist("initial conditions").sublist("permeability") = empty;
  }

  return out_list;
}


/* ******************************************************************
* Select proper evaluator based on the list of input parameters.
****************************************************************** */
void InputConverterU::TranslateFieldEvaluator_(
    DOMNode* node, const std::string& field, const std::string& unit,
    const std::string& reg_str, const std::vector<std::string>& regions,
    Teuchos::ParameterList& out_ic, Teuchos::ParameterList& out_ev,
    std::string data_key, std::string domain)
{
  std::string type = GetAttributeValueS_(node, "type", TYPE_NONE, false, "");
  if (type == "file") {  // Amanzi restart file
    std::string filename = GetAttributeValueS_(node, "filename");
    Teuchos::ParameterList& field_ic = out_ic.sublist(field);
    field_ic.set<std::string>("restart file", filename);

    Teuchos::ParameterList& field_ev = out_ev.sublist(field);
    field_ev.set<std::string>("field evaluator type", "constant variable");
  } else if (type == "h5file") {  // regular h5 file
    std::string filename = GetAttributeValueS_(node, "filename");
    bool temporal = GetAttributeValueS_(node, "constant_in_time", TYPE_NUMERICAL, false, "true") == "true";

    Teuchos::ParameterList& field_ev = out_ev.sublist(field);
    field_ev.set<std::string>("field evaluator type", "independent variable from file")
        .set<std::string>("filename", filename)
        .set<std::string>("domain name", domain)
        .set<std::string>("component name", "cell")
        .set<std::string>("mesh entity", "cell")
        .set<std::string>("variable name", field)
        .set<int>("number of dofs", 1)
        .set<bool>("constant in time", temporal);
  } else {
    double val = GetAttributeValueD_(node, data_key.c_str(), TYPE_NUMERICAL, DVAL_MIN, DVAL_MAX, unit);

    Teuchos::ParameterList& field_ev = out_ev.sublist(field);
    field_ev.sublist("function").sublist(reg_str)
        .set<Teuchos::Array<std::string> >("regions", regions)
        .set<std::string>("component", "cell")
        .sublist("function").sublist("function-constant")
        .set<double>("value", val);
    field_ev.set<std::string>("field evaluator type", "independent variable");
  }
}


/* ******************************************************************
* Select proper IC based on the list of input parameters.
****************************************************************** */
void InputConverterU::TranslateFieldIC_(
    DOMNode* node, std::string field, std::string unit,
    const std::string& reg_str, const std::vector<std::string>& regions,
    Teuchos::ParameterList& out_ic, Teuchos::ParameterList& out_ev,
    std::string data_key)
{
  std::string type = GetAttributeValueS_(node, "type", TYPE_NONE, false, "");
  if (type == "file") {
    std::string filename = GetAttributeValueS_(node, "filename");
    Teuchos::ParameterList& field_ic = out_ic.sublist(field);
    field_ic.set<std::string>("restart file", filename);
  } else {
    double val = GetAttributeValueD_(node, data_key.c_str(), TYPE_NUMERICAL, DVAL_MIN, DVAL_MAX, unit);

    Teuchos::ParameterList& field_ic = out_ic.sublist(field);
    field_ic.sublist("function").sublist(reg_str)
        .set<Teuchos::Array<std::string> >("regions", regions)
        .set<std::string>("component", "cell")
        .sublist("function").sublist("function-constant")
        .set<double>("value", val);
  }
}


/* ******************************************************************
* Mesh patition sublist based on materials
****************************************************************** */
Teuchos::ParameterList InputConverterU::TranslateMaterialsPartition_()
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

        for (int n = 0; n < names.size(); ++n) {
          regions.push_back(names[n]);
        } 
      }
    }
  }
  tmp_list.set<Teuchos::Array<std::string> >("region list", regions);
  
  return out_list;
}


/* ******************************************************************
* Create initialization list for concentration. This routine is called
* when geochemistry list exists for initial conditions.
****************************************************************** */
void InputConverterU::TranslateStateICsAmanziGeochemistry_(
    Teuchos::ParameterList& out_list, std::string& constraint,
    std::vector<std::string>& regions)
{
  if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH) {
    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "Compatibility mode: translating ICs for native chemistry" << std::endl;
  }

  bool flag;
  DOMNode* node;
  DOMElement* element;

  node = GetUniqueElementByTagsString_("process_kernels, chemistry", flag);
  std::string engine = GetAttributeValueS_(node, "engine");

  node = GetUniqueElementByTagsString_("geochemistry, constraints", flag);
  if (flag && engine == "amanzi") {
    std::string name;
    element = GetUniqueChildByAttribute_(node, "name", constraint, flag, true);
    std::vector<DOMNode*> children = GetSameChildNodes_(element, name, flag);
    if (children.size() != phases_["water"].size()) {
      Errors::Message msg;
      msg << "Constraint \"" << constraint << "\" is not backward compatible: "
          << " check the number of components.";
      Exceptions::amanzi_throw(msg);
    }

    Teuchos::ParameterList& ic_list = out_list.sublist("total_component_concentration")
        .sublist("function").sublist("All");

    ic_list.set<Teuchos::Array<std::string> >("regions", regions)
        .set<std::string>("component", "cell");

    Teuchos::ParameterList& tmp_list = ic_list.sublist("function")
        .set<int>("number of dofs", children.size())
        .set<std::string>("function type", "composite function");

    for (int i = 0; i < children.size(); ++i) {
      std::string species = GetAttributeValueS_(children[i], "name");
      double val = GetAttributeValueD_(children[i], "value", TYPE_NUMERICAL, DVAL_MIN, DVAL_MAX);

      // find position of species in the list of component names
      int k(-1);
      for (int n = 0; n < comp_names_all_.size(); ++n) {
        if (comp_names_all_[n] == species) {
          k = n;
          break;
        }
      }

      std::stringstream dof_str;
      dof_str << "dof " << k+1 << " function";
      tmp_list.sublist(dof_str.str()).sublist("function-constant")
                                     .set<double>("value", val);
    }
  }
}

}  // namespace AmanziInput
}  // namespace Amanzi
