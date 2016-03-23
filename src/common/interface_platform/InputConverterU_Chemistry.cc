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
#include <sys/stat.h>

//TPLs
#include "Teuchos_ParameterList.hpp"
#include "xercesc/dom/DOM.hpp"

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
Teuchos::ParameterList InputConverterU::TranslateChemistry_()
{
  Teuchos::ParameterList out_list;

  if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH)
    *vo_->os() << "Translating chemistry" << std::endl;

  MemoryManager mm;
  char* text;
  DOMNodeList *node_list, *children;
  DOMNode* node;
  DOMElement* element;

  // chemical engine
  bool flag;
  node = GetUniqueElementByTagsString_("process_kernels, chemistry", flag);
  std::string engine = GetAttributeValueS_(static_cast<DOMElement*>(node), "engine");

  // process engine
  bool native(false);
  if (engine ==  "amanzi") {
    native = true;
    out_list.set<std::string>("chemistry model", "Amanzi");

    std::string bgdfilename, format("simple");
    node = GetUniqueElementByTagsString_("process_kernels, chemistry", flag);
    if (flag) {
      element = static_cast<DOMElement*>(node);
      bgdfilename = GetAttributeValueS_(element, "input_filename");
      format = GetAttributeValueS_(element, "format", TYPE_NONE, false, format);
    } else {
      bgdfilename = CreateBGDFile(xmlfilename_);
    }

    Teuchos::ParameterList& bgd_list = out_list.sublist("Thermodynamic Database");
    bgd_list.set<std::string>("file", bgdfilename);
    bgd_list.set<std::string>("format", format);

  } else {
    bool valid_engine(true);
    std::string file_location;

    if (engine == "pflotran") {
      out_list.set<std::string>("Engine", "PFloTran");
      file_location = "process_kernels, chemistry";
    } else if (engine == "crunchflow") {
      out_list.set<std::string>("Engine", "CrunchFlow");
      file_location = "process_kernels, chemistry";
    } else {
      valid_engine = false;
    }

    // Pass along chemistry engine info.
    if (valid_engine) {
      out_list.set<std::string>("chemistry model", "Alquimia");

      // Find the name of the engine-specific input file.
      node = GetUniqueElementByTagsString_(file_location, flag);
      if (flag) {
        element = static_cast<DOMElement*>(node);
        std::string inpfilename = GetAttributeValueS_(element, "input_filename");
        out_list.set<std::string>("Engine Input File", inpfilename);
      } else {
        Errors::Message msg;
        msg << "Unique tag string \"" << file_location << "\" must exists.\n";
        Exceptions::amanzi_throw(msg);
      }
    }
  }
  
  // minerals
  std::vector<std::string> minerals;
  node = GetUniqueElementByTagsString_("phases, solid_phase, minerals", flag);
  if (flag) {
    children = static_cast<DOMElement*>(node)->getElementsByTagName(mm.transcode("mineral"));

    for (int i = 0; i < children->getLength(); ++i) {
      DOMNode* inode = children->item(i);
      std::string name = TrimString_(mm.transcode(inode->getTextContent()));
      minerals.push_back(name);
    }
  }

  // region specific initial conditions
  std::vector<std::string> sorption_sites;
  Teuchos::ParameterList& ic_list = out_list.sublist("initial conditions");

  node_list = doc_->getElementsByTagName(mm.transcode("materials"));
  element = static_cast<DOMElement*>(node_list->item(0));
  children = element->getElementsByTagName(mm.transcode("material"));

  for (int i = 0; i < children->getLength(); ++i) {
    DOMNode* inode = children->item(i);

    node = GetUniqueElementByTagsString_(inode, "assigned_regions", flag);
    text = mm.transcode(node->getTextContent());
    std::vector<std::string> regions = CharToStrings_(text);

    // mineral volume fraction and specific surface area.
    if (minerals.size() > 0) {
      out_list.set<Teuchos::Array<std::string> >("minerals", minerals);

      // if (pk_model_["chemistry"] == "amanzi") {
      Teuchos::ParameterList& volfrac = ic_list.sublist("mineral_volume_fractions");
      Teuchos::ParameterList& surfarea = ic_list.sublist("mineral_specific_surface_area");

      for (std::vector<std::string>::const_iterator it = regions.begin(); it != regions.end(); it++) {
        Teuchos::ParameterList& aux1_list = volfrac.sublist("function").sublist(*it)
            .set<std::string>("region", *it)
            .set<std::string>("component", "cell")
            .sublist("function");
        aux1_list.set<int>("number of dofs", minerals.size())
            .set("function type", "composite function");

        Teuchos::ParameterList& aux2_list = surfarea.sublist("function").sublist(*it)
            .set<std::string>("region", *it)
            .set<std::string>("component", "cell")
            .sublist("function");
        aux2_list.set<int>("number of dofs", minerals.size())
            .set("function type", "composite function");

        for (int j = 0; j < minerals.size(); ++j) {
          std::stringstream ss;
          ss << "DoF " << j + 1 << " Function";
 
          node = GetUniqueElementByTagsString_(inode, "minerals", flag);
          double mvf(0.0), msa(0.0);
          if (flag) {
            element = GetUniqueChildByAttribute_(node, "name", minerals[j], flag, true);
            mvf = GetAttributeValueD_(element, "volume_fraction", TYPE_NUMERICAL, false, 0.0);
            msa = GetAttributeValueD_(element, "specific_surface_area", TYPE_NUMERICAL, false, 0.0);
          }
          aux1_list.sublist(ss.str()).sublist("function-constant").set<double>("value", mvf);
          aux2_list.sublist(ss.str()).sublist("function-constant").set<double>("value", msa);
        }
      }
    }

    // ion exchange
    node = GetUniqueElementByTagsString_(inode, "ion_exchange, cec", flag);
    if (flag) {
      Teuchos::ParameterList& sites = ic_list.sublist("ion_exchange_sites");
      double cec = std::strtod(mm.transcode(node->getTextContent()), NULL);

      for (std::vector<std::string>::const_iterator it = regions.begin(); it != regions.end(); ++it) {
        sites.sublist("function").sublist(*it)
            .set<std::string>("region", *it)
            .set<std::string>("component", "cell")
            .sublist("function").sublist("function-constant")
            .set<double>("value", cec);
      }
    }

    // sorption 
    int nsolutes = phases_["water"].size();
    node = GetUniqueElementByTagsString_(inode, "sorption_isotherms", flag);
    if (flag && nsolutes > 0) {
      Teuchos::ParameterList& kd = ic_list.sublist("isotherm_kd");
      Teuchos::ParameterList& langmuir_b = ic_list.sublist("isotherm_langmuir_b");
      Teuchos::ParameterList& freundlich_n = ic_list.sublist("isotherm_freundlich_n");

      for (std::vector<std::string>::const_iterator it = regions.begin(); it != regions.end(); ++it) {
        Teuchos::ParameterList& aux1_list = kd.sublist("function").sublist(*it)
            .set<std::string>("region", *it)
            .set<std::string>("component", "cell")
            .sublist("function");
        aux1_list.set<int>("number of dofs", nsolutes).set("function type", "composite function");

        Teuchos::ParameterList& aux2_list = langmuir_b.sublist("function").sublist(*it)
            .set<std::string>("region", *it)
            .set<std::string>("component", "cell")
            .sublist("function");
        aux2_list.set<int>("number of dofs", nsolutes).set("function type", "composite function");

        Teuchos::ParameterList& aux3_list = langmuir_b.sublist("function").sublist(*it)
            .set<std::string>("region", *it)
            .set<std::string>("component", "cell")
            .sublist("function");
        aux3_list.set<int>("number of dofs", nsolutes).set("function type", "composite function");

        for (int j = 0; j < nsolutes; ++j) {
          std::string solute_name = phases_["water"][j];
          element = GetUniqueChildByAttribute_(node, "name", solute_name, flag, false);
          if (flag) {
            DOMNode* jnode = GetUniqueElementByTagsString_(element, "kd_model", flag);
            element = static_cast<DOMElement*>(jnode);
          }

          std::stringstream ss;
          ss << "DoF " << j + 1 << " Function";

          double val = (!flag) ? 0.0 : GetAttributeValueD_(element, "kd", TYPE_NUMERICAL, false, 0.0);
          aux1_list.sublist(ss.str()).sublist("function-constant").set<double>("value", val);

          val = (!flag) ? 0.0 : GetAttributeValueD_(element, "b", TYPE_NUMERICAL, false, 0.0);
          aux2_list.sublist(ss.str()).sublist("function-constant").set<double>("value", val);

          val = (!flag) ? 0.0 : GetAttributeValueD_(element, "n", TYPE_NUMERICAL, false, 0.0);
          aux3_list.sublist(ss.str()).sublist("function-constant").set<double>("value", val);
        }
      }
    }

    // surface complexation
    node = GetUniqueElementByTagsString_(inode, "surface_complexation", flag);
    if (flag) {
      // sorption_sites.push_back("siteA");  // no translation rules so far
      Teuchos::ParameterList& complexation = ic_list.sublist("surface_complexation");

      double val = GetAttributeValueD_(static_cast<DOMElement*>(node), "density");

      for (std::vector<std::string>::const_iterator it = regions.begin(); it != regions.end(); ++it) {
        complexation.sublist("function").sublist(*it)
            .set<std::string>("region", *it)
            .set<std::string>("component", "cell")
            .sublist("function").sublist("function-constant")
            .set<double>("value", val);
      }
    }

    // free ion
    {
      Teuchos::ParameterList& free_ion = ic_list.sublist("free_ion_species");

      for (std::vector<std::string>::const_iterator it = regions.begin(); it != regions.end(); ++it) {
        Teuchos::ParameterList& aux1_list = free_ion.sublist("function").sublist(*it)
            .set<std::string>("region", *it)
            .set<std::string>("component", "cell")
            .sublist("function");
        aux1_list.set<int>("number of dofs", nsolutes).set("function type", "composite function");

        for (int j = 0; j < nsolutes; ++j) {
          std::string solute_name = phases_["water"][j];

          std::stringstream ss;
          ss << "DoF " << j + 1 << " Function";

          aux1_list.sublist(ss.str()).sublist("function-constant").set<double>("value", 1e-9);
        }
      }
    }
  }

  // general parameters
  int max_itrs(100), cut_threshold(8), increase_threshold(4);
  double tol(1e-12), dt_max(1e+10), dt_min(1e+10), dt_init(1e+7), dt_cut(2.0), dt_increase(1.2);
  std::string activity_model("unit"), dt_method("fixed");
  std::vector<std::string> aux_data;

  node = GetUniqueElementByTagsString_("unstructured_controls, unstr_chemistry_controls", flag);
  if (flag) {
    children = node->getChildNodes();
    for (int i = 0; i < children->getLength(); ++i) {
      DOMNode* inode = children->item(i);
      if (inode->getNodeType() != DOMNode::ELEMENT_NODE) continue;

      text = mm.transcode(inode->getNodeName());
      if (strcmp(text, "activity_model") == 0) {
        activity_model = GetTextContentS_(inode, "unit, debye-huckel");
      } else if (strcmp(text, "maximum_newton_iterations") == 0) {
        max_itrs = strtol(mm.transcode(inode->getTextContent()), NULL, 10);
      } else if (strcmp(text, "tolerance") == 0) {
        tol = strtod(mm.transcode(inode->getTextContent()), NULL);
      } else if (strcmp(text, "min_time_step") == 0) {
        dt_min = strtod(mm.transcode(inode->getTextContent()), NULL);
      } else if (strcmp(text, "max_time_step") == 0) {
        dt_max = strtod(mm.transcode(inode->getTextContent()), NULL);
      } else if (strcmp(text, "initial_time_step") == 0) {
        dt_init = strtod(mm.transcode(inode->getTextContent()), NULL);
      } else if (strcmp(text, "time_step_control_method") == 0) {
        dt_method = mm.transcode(inode->getTextContent());
      } else if (strcmp(text, "time_step_cut_threshold") == 0) {
        cut_threshold = strtol(mm.transcode(inode->getTextContent()), NULL, 10);
      } else if (strcmp(text, "time_step_cut_factor") == 0) {
        dt_cut = strtod(mm.transcode(inode->getTextContent()), NULL);
      } else if (strcmp(text, "time_step_increase_threshold") == 0) {
        increase_threshold = strtol(mm.transcode(inode->getTextContent()), NULL, 10);
      } else if (strcmp(text, "time_step_increase_factor") == 0) {
        dt_increase = strtod(mm.transcode(inode->getTextContent()), NULL);
      } else if (strcmp(text, "auxiliary_data") == 0) {
        aux_data = CharToStrings_(mm.transcode(inode->getTextContent()));
      }
    }
  }
  out_list.set<std::string>("activity model", activity_model);
  out_list.set<int>("maximum Newton iterations", max_itrs);
  out_list.set<double>("tolerance", tol);
  out_list.set<double>("max time step (s)", dt_max);
  out_list.set<double>("min time step (s)", dt_min);
  out_list.set<double>("initial time step (s)", dt_init);
  out_list.set<int>("time step cut threshold", cut_threshold);
  out_list.set<double>("time step cut factor", dt_cut);
  out_list.set<int>("time step increase threshold", increase_threshold);
  out_list.set<double>("time step increase factor", dt_increase);
  out_list.set<std::string>("time step control method", dt_method);
  if (aux_data.size() > 0)
      out_list.set<Teuchos::Array<std::string> >("auxiliary data", aux_data);
  if (sorption_sites.size() > 0)
      out_list.set<Teuchos::Array<std::string> >("sorption sites", sorption_sites);

  // miscalleneous
  out_list.set<int>("number of component concentrations", comp_names_all_.size());

  out_list.sublist("VerboseObject") = verb_list_.sublist("VerboseObject");
  return out_list;
}


/* ******************************************************************
* Adds kd values to a bgd file (using the first material).
* Returns the name of this file.  
****************************************************************** */
std::string InputConverterU::CreateBGDFile(std::string& filename)
{
  DOMNode* node;
  DOMElement* element;

  std::ofstream bgd_file;
  std::stringstream species;
  std::stringstream isotherms;

  bool flag;
  node = GetUniqueElementByTagsString_("materials, material, sorption_isotherms", flag);
  if (flag) {
    std::string name;
    std::vector<DOMNode*> children = GetSameChildNodes_(node, name, flag, false);
    for (int i = 0; i < children.size(); ++i) {
      DOMNode* inode = children[i];
      element = static_cast<DOMElement*>(inode);
      name = GetAttributeValueS_(element, "name");
      species << name << " ;   0.00 ;   0.00 ;   1.00 \n";

      DOMNode* knode = GetUniqueElementByTagsString_(inode, "kd_model", flag);
      element = static_cast<DOMElement*>(knode);
      if (flag) {
        double kd = GetAttributeValueD_(element, "kd");
        std::string model = GetAttributeValueS_(element, "model");
        if (model == "langmuir") {
          double b = GetAttributeValueD_(element, "b");
          isotherms << name << " ; langmuir ; " << kd << " " << b << std::endl;
        } else if (model == "freundlich") {
          double n = GetAttributeValueD_(element, "n");
          isotherms << name << " ; freundlich ; " << kd << " " << n << std::endl;
        } else {
          isotherms << name << " ; linear ; " << kd << std::endl;
        }
      }
    }
  }
  
  // build bgd filename
  std::string bgdfilename(filename);
  std::string new_extension(".bgd");
  size_t pos = bgdfilename.find(".xml");
  bgdfilename.replace(pos, (size_t)4, new_extension, (size_t)0, (size_t)4);

  // open output bgd file
  if (rank_ == 0) {
    struct stat buffer;
    int status = stat(bgdfilename.c_str(), &buffer);
    if (!status) {
      if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH) {
        Teuchos::OSTab tab = vo_->getOSTab();
        *vo_->os() << "File \"" << bgdfilename.c_str() 
                   << "\" exists, skipping its creation." << std::endl;
      }
    } else {
      bgd_file.open(bgdfilename.c_str());

      // <Primary Species
      bgd_file << "<Primary Species\n";
      bgd_file << species.str();

      //<Isotherms
      bgd_file << "<Isotherms\n";
      bgd_file << "# Note, these values will be overwritten by the xml file\n";
      bgd_file << isotherms.str();

      bgd_file.close();
    }
  }

  return bgdfilename;
}

}  // namespace AmanziInput
}  // namespace Amanzi


