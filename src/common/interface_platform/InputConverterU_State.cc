/*
  This is the input component of the Amanzi code. 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Markus Berndt (original version)
           Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <sstream>
#include <string>

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
  gravity[dim_-1] = -GRAVITY_MAGNITUDE;
  out_ic.sublist("gravity").set<Teuchos::Array<double> >("value", gravity);

  // --- viscosity
  bool flag;
  DOMNode* node = GetUniqueElementByTagsString_("phases, liquid_phase, viscosity", flag);
  text_content = mm.transcode(node->getTextContent());
  double viscosity = std::strtod(text_content, NULL);
  out_ic.sublist("fluid_viscosity").set<double>("value", viscosity);

  // --- constant density
  node = GetUniqueElementByTagsString_("phases, liquid_phase, density", flag);
  text_content = mm.transcode(node->getTextContent());
  rho_ = std::strtod(text_content, NULL);
  out_ic.sublist("fluid_density").set<double>("value", rho_);

  out_ic.sublist("mass_density_liquid").sublist("function").sublist("All")
      .set<std::string>("region", "All")
      .set<std::string>("component", "cell")
      .sublist("function").sublist("function-constant")
      .set<double>("value", rho_);

  // --- region specific initial conditions from material properties
  std::map<std::string, int> reg2mat;
  int mat(0);

  DOMNodeList* node_list = doc_->getElementsByTagName(mm.transcode("materials"));
  DOMNodeList* children = node_list->item(0)->getChildNodes();

  for (int i = 0; i < children->getLength(); i++) {
    DOMNode* inode = children->item(i);
    if (DOMNode::ELEMENT_NODE == inode->getNodeType()) {
      std::string mat_name = GetAttributeValueS_(static_cast<DOMElement*>(inode), "name");

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

      // create regions string
      std::string reg_str;
      for (std::vector<std::string>::const_iterator it = regions.begin(); it != regions.end(); ++it) {
        reg_str = reg_str + *it;
      }

      // -- porosity: skip if compressibility model was already provided.
      if (!compressibility_) {
        double porosity;
        node = GetUniqueElementByTagsString_(inode, "mechanical_properties, porosity", flag);
        if (flag) {
          porosity = GetAttributeValueD_(static_cast<DOMElement*>(node), "value");
        } else {
          msg << "Porosity element must be specified under mechanical_properties";
          Exceptions::amanzi_throw(msg);
        }
        Teuchos::ParameterList& porosity_ev = out_ev.sublist("porosity");
        porosity_ev.sublist("function").sublist(reg_str)
            .set<Teuchos::Array<std::string> >("regions", regions)
            .set<std::string>("component", "cell")
            .sublist("function").sublist("function-constant")
            .set<double>("value", porosity);
        porosity_ev.set<std::string>("field evaluator type", "independent variable");
      }

      // -- permeability.
      double perm_x, perm_y, perm_z;
      bool perm_init_from_file(false), conductivity(false);
      std::string perm_file, perm_attribute, perm_format;

      node = GetUniqueElementByTagsString_(inode, "permeability", flag);
      if (!flag) {
        conductivity = true;
        node = GetUniqueElementByTagsString_(inode, "hydraulic_conductivity", flag);
      }

      // first we get either permeability value or the file name
      int file(0);
      char* file_name;
      char* attr_name;
      double kx(-1.0), ky(-1.0), kz(-1.0);

      DOMNamedNodeMap* attr_tmp = node->getAttributes();
      for (int k=0; k < attr_tmp->getLength(); k++) {
        DOMNode* knode = attr_tmp->item(k);

        if (DOMNode::ATTRIBUTE_NODE == knode->getNodeType()) {
          tagname = mm.transcode(knode->getNodeName());
          text_content = mm.transcode(knode->getTextContent());

          if (strcmp(tagname, "x") == 0) {
            kx = std::strtod(text_content, NULL);
          } else if (strcmp(tagname, "y") == 0) {
            ky = std::strtod(text_content, NULL);
          } else if (strcmp(tagname, "z") == 0) {
            kz = std::strtod(text_content, NULL);
          } else if (strcmp(tagname, "type") == 0) {
            file++;
          } else if (strcmp(tagname, "filename") == 0) {
            file++;
            file_name = new char[std::strlen(text_content) + 1];
            std::strcpy(file_name, text_content);
          } else if (strcmp(tagname, "attribute") == 0) {
            file++;
            attr_name = new char[std::strlen(text_content) + 1];
            std::strcpy(attr_name, text_content);
          }
        }
      }

      if (conductivity) {
        kx *= viscosity / (rho_ * GRAVITY_MAGNITUDE);
        ky *= viscosity / (rho_ * GRAVITY_MAGNITUDE);
        kz *= viscosity / (rho_ * GRAVITY_MAGNITUDE);
      }

      // Second, we copy collected data to XML file.
      // For now permeability is not dumped into checkpoint files.
      Teuchos::ParameterList& permeability_ic = out_ic.sublist("permeability");
      permeability_ic.set<bool>("write checkpoint", false);

      if (file == 3) {
        permeability_ic.sublist("exodus file initialization")
            .set<std::string>("file", file_name)
            .set<std::string>("attribute", attr_name);
        delete file_name;
        delete attr_name;
      } else if (file == 0) {
        if (ky < 0) ky = kz;  // x-z system was defined
        Teuchos::ParameterList& aux_list = permeability_ic.sublist("function").sublist(reg_str)
            .set<Teuchos::Array<std::string> >("regions", regions)
            .set<std::string>("component", "cell")
            .sublist("function");
        aux_list.set<int>("Number of DoFs", dim_)
            .set<std::string>("Function type", "composite function");
        aux_list.sublist("DoF 1 Function").sublist("function-constant").set<double>("value", kx);
        aux_list.sublist("DoF 2 Function").sublist("function-constant").set<double>("value", ky);
        if (dim_ == 3) {
          aux_list.sublist("DoF 3 Function").sublist("function-constant").set<double>("value", kz);
        } else {
          kz = 0.0;
        }
      } else {
        ThrowErrorIllformed_("materials", "permeability/hydraulic conductivity", "file/filename/attribute");
      }
      if (kx < 0.0 || ky < 0.0 || kz < 0.0) {
        ThrowErrorIllformed_("materials", "permeability/hydraulic conductivity", "file/filename/attribute");
      }

      // -- specific_yield
      node = GetUniqueElementByTagsString_(inode, "mechanical_properties, specific_yield", flag);
      if (flag) {
        double specific_yield = GetAttributeValueD_(static_cast<DOMElement*>(node), "value");

        Teuchos::ParameterList& spec_yield_ic = out_ic.sublist("specific_yield");
        spec_yield_ic.sublist("function").sublist(reg_str)
            .set<Teuchos::Array<std::string> >("regions",regions)
            .set<std::string>("component", "cell")
            .sublist("function").sublist("function-constant")
            .set<double>("value", specific_yield);
      }

      // -- specific storage
      node = GetUniqueElementByTagsString_(inode, "mechanical_properties, specific_storage", flag);
      if (flag) {
        double specific_storage = GetAttributeValueD_(static_cast<DOMElement*>(node), "value");

        Teuchos::ParameterList& spec_yield_ic = out_ic.sublist("specific_storage");
        spec_yield_ic.sublist("function").sublist(reg_str)
            .set<Teuchos::Array<std::string> >("regions",regions)
            .set<std::string>("component", "cell")
            .sublist("function").sublist("function-constant")
            .set<double>("value", specific_storage);
      }

      // -- particle density
      node = GetUniqueElementByTagsString_(inode, "mechanical_properties, particle_density", flag);
      if (flag) {
        double particle_density = GetAttributeValueD_(static_cast<DOMElement*>(node), "value");

        Teuchos::ParameterList& part_dens_ev = out_ev.sublist("particle_density");
        part_dens_ev.sublist("function").sublist(reg_str)
            .set<Teuchos::Array<std::string> >("regions", regions)
            .set<std::string>("component", "cell")
            .sublist("function").sublist("function-constant")
            .set<double>("value", particle_density);
        part_dens_ev.set<std::string>("field evaluator type", "independent variable");
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
      std::string reg_str;
      for (std::vector<std::string>::const_iterator it = regions.begin(); it != regions.end(); ++it) {
        reg_str = reg_str + *it;
      }

      // -- uniform pressure
      node = GetUniqueElementByTagsString_(inode, "liquid_phase, liquid_component, uniform_pressure", flag);
      if (flag) {
        double p = GetAttributeValueD_(static_cast<DOMElement*>(node), "value");

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
        double p = GetAttributeValueD_(static_cast<DOMElement*>(node), "value");
        std::vector<double> grad = GetAttributeVector_(static_cast<DOMElement*>(node), "gradient");
        std::vector<double> refc = GetAttributeVector_(static_cast<DOMElement*>(node), "reference_coord");

        Teuchos::Array<double> grad_with_time(grad.size() + 1);
        Teuchos::Array<double> refc_with_time(grad.size() + 1);

        grad_with_time[0] = 0.0;
        refc_with_time[0] = 0.0;

        for (int j = 0; j != grad.size(); ++j) {
          grad_with_time[j + 1] = grad[j];
          refc_with_time[j + 1] = refc[j];
        }

        Teuchos::ParameterList& pressure_ic = out_ic.sublist("pressure");
        pressure_ic.sublist("function").sublist(reg_str)
            .set<Teuchos::Array<std::string> >("regions", regions)
            .set<std::string>("component", "cell")
            .sublist("function").sublist("function-linear")
            .set<double>("y0", p)
            .set<Teuchos::Array<double> >("x0", refc_with_time)
            .set<Teuchos::Array<double> >("gradient", grad_with_time);
      }

      // -- uniform saturation
      node = GetUniqueElementByTagsString_(inode, "liquid_phase, liquid_component, uniform_saturation", flag);
      if (flag) {
        double s = GetAttributeValueD_(static_cast<DOMElement*>(node), "value");

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
        double s = GetAttributeValueD_(static_cast<DOMElement*>(node), "value");
        std::vector<double> grad = GetAttributeVector_(static_cast<DOMElement*>(node), "gradient");
        std::vector<double> refc = GetAttributeVector_(static_cast<DOMElement*>(node), "reference_coord");

        Teuchos::Array<double> grad_with_time(grad.size() + 1);
        Teuchos::Array<double> refc_with_time(grad.size() + 1);

        grad_with_time[0] = 0.0;
        refc_with_time[0] = 0.0;

        for (int j = 0; j != grad.size(); ++j) {
          grad_with_time[j + 1] = grad[j];
          refc_with_time[j + 1] = refc[j];
        }

        Teuchos::ParameterList& saturation_ic = out_ic.sublist("saturation_liquid");
        saturation_ic.sublist("function").sublist(reg_str)
            .set<Teuchos::Array<std::string> >("regions", regions)
            .set<std::string>("component", "cell")
            .sublist("function").sublist("function-linear")
            .set<double>("y0", s)
            .set<Teuchos::Array<double> >("x0", refc_with_time)
            .set<Teuchos::Array<double> >("gradient", grad_with_time);
      }

      // -- darcy_flux
      node = GetUniqueElementByTagsString_(inode, "liquid_phase, liquid_component, velocity", flag);
      if (flag) {
        std::vector<double> velocity;
        velocity.push_back(GetAttributeValueD_(static_cast<DOMElement*>(node), coords_[0].c_str()));
        velocity.push_back(GetAttributeValueD_(static_cast<DOMElement*>(node), coords_[1].c_str()));
        if (dim_ == 3) velocity.push_back(GetAttributeValueD_(static_cast<DOMElement*>(node), coords_[2].c_str()));

        Teuchos::ParameterList& darcy_flux_ic = out_ic.sublist("darcy_flux");
        Teuchos::ParameterList& tmp_list =
            darcy_flux_ic.set<bool>("dot with normal", true)
            .sublist("function").sublist(reg_str)
            .set<Teuchos::Array<std::string> >("regions", regions)
            .set<std::string>("component", "face")
            .sublist("function")
            .set<int>("Number of DoFs", dim_)
            .set<std::string>("Function type", "composite function");

        for (int k = 0; k != dim_; ++k) {
          std::stringstream dof_str;
          dof_str << "DoF " << k+1 << " Function";
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

        DOMNodeList* children = node->getChildNodes();
        for (int j = 0; j < children->getLength(); ++j) {
          DOMNode* jnode = children->item(j);
          tagname = mm.transcode(jnode->getNodeName());

          if (strcmp(tagname, "uniform_conc") == 0) {
            std::string text = GetAttributeValueS_(static_cast<DOMElement*>(jnode), "name");
            int m = GetPosition_(phases_["water"], text);
            vals[m] = GetAttributeValueD_(static_cast<DOMElement*>(jnode), "value");
          }
        }

        Teuchos::ParameterList& tcc_ic = out_ic.sublist("total_component_concentration");
        Teuchos::ParameterList& dof_list = tcc_ic.sublist("function").sublist(reg_str)
            .set<Teuchos::Array<std::string> >("regions", regions)
            .set<std::string>("component", "cell")
            .sublist("function")
            .set<int>("Number of DoFs", ncomp_all)
            .set<std::string>("Function type", "composite function");

        for (int k = 0; k < ncomp_l; k++) {
          std::string name = phases_["water"][k];
          std::stringstream dof_str;
          dof_str << "DoF " << k + 1 << " Function";
          dof_list.sublist(dof_str.str()).sublist("function-constant").set<double>("value", vals[k]);
        }
      }

      // -- total_component_concentration (gas phase)
      node = GetUniqueElementByTagsString_(inode, "gas_phase, solute_component", flag);
      if (flag) {
        std::vector<double> vals(ncomp_g, 0.0);

        DOMNodeList* children = node->getChildNodes();
        for (int j = 0; j < children->getLength(); ++j) {
          DOMNode* jnode = children->item(j);
          tagname = mm.transcode(jnode->getNodeName());

          if (strcmp(tagname, "uniform_conc") == 0) {
            std::string text = GetAttributeValueS_(static_cast<DOMElement*>(jnode), "name");
            int m = GetPosition_(phases_["air"], text);
            vals[m] = GetAttributeValueD_(static_cast<DOMElement*>(jnode), "value");
          }
        }

        Teuchos::ParameterList& tcc_ic = out_ic.sublist("total_component_concentration");
        Teuchos::ParameterList& dof_list = tcc_ic.sublist("function").sublist(reg_str).sublist("function");
        for (int k = 0; k < ncomp_g; k++) {
          std::string name = phases_["air"][k];
          std::stringstream dof_str;
          dof_str << "DoF " << ncomp_l + k + 1 << " Function";
          dof_list.sublist(dof_str.str()).sublist("function-constant").set<double>("value", vals[k]);
        }
      }

      // -- uniform temperature
      node = GetUniqueElementByTagsString_(inode, "uniform_temperature", flag);
      if (flag) {
        double val = GetAttributeValueD_(static_cast<DOMElement*>(node), "value");

        Teuchos::ParameterList& temperature_ic = out_ic.sublist("temperature");
        temperature_ic.sublist("function").sublist(reg_str)
            .set<Teuchos::Array<std::string> >("regions", regions)
            .set<std::string>("component", "cell")
            .sublist("function").sublist("function-constant")
            .set<double>("value", val);
      }

      // -- geochemical condition
      node = GetUniqueElementByTagsString_(inode, "liquid_phase, geochemistry, constraint", flag);
      if (flag) {
        std::string name = GetAttributeValueS_(static_cast<DOMElement*>(node), "name");

        out_ic.sublist("geochemical conditions").sublist(name)
            .set<Teuchos::Array<std::string> >("regions", regions);
      }
    }
  }

  // atmospheric pressure
  out_ic.sublist("atmospheric_pressure").set<double>("value", ATMOSPHERIC_PRESSURE);

  // add mesh partitions to the state list
  out_list.sublist("mesh partitions") = TranslateMaterialsPartition_();

  return out_list;
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
      DOMNamedNodeMap* attr_map = inode->getAttributes();

      DOMNode* node = GetUniqueElementByTagsString_(inode, "assigned_regions", flag);
      if (flag) {
        char* text_content = mm.transcode(node->getTextContent());
        std::vector<std::string> names = CharToStrings_(text_content);

        for (int i = 0; i < names.size(); i++) {
          regions.push_back(names[i]);
        } 
      }
    }
  }
  tmp_list.set<Teuchos::Array<std::string> >("region list", regions);
  
  return out_list;
}

}  // namespace AmanziInput
}  // namespace Amanzi
