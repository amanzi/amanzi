/*
  Input Converter

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
      bgdfilename = GetAttributeValueS_(element, "input_filename", TYPE_NONE, false, "");
      if (bgdfilename == "") {
        bgdfilename = CreateBGDFile(xmlfilename_);
      }
      format = GetAttributeValueS_(element, "format", TYPE_NONE, false, format);
    }

    Teuchos::ParameterList& bgd_list = out_list.sublist("thermodynamic database");
    bgd_list.set<std::string>("file", bgdfilename);
    bgd_list.set<std::string>("format", format);

  } else {
    bool valid_engine(true);
    std::string file_location;

    if (engine == "pflotran") {
      out_list.set<std::string>("engine", "PFloTran");
      file_location = "process_kernels, chemistry";
    } else if (engine == "crunchflow") {
      out_list.set<std::string>("engine", "CrunchFlow");
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
	if (element->hasAttribute(mm.transcode("input_filename"))) {
            std::string inpfilename = GetAttributeValueS_(element, "input_filename");
            out_list.set<std::string>("engine input file", inpfilename);
	} else {
	    std::string inpfilename = CreateINFile(xmlfilename_);
            out_list.set<std::string>("engine input file", inpfilename);
	}
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
          ss << "dof " << j + 1 << " function";
 
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
    node = GetUniqueElementByTagsString_(inode, "ion_exchange, cations", flag);
    if (flag) {
      Teuchos::ParameterList& sites = ic_list.sublist("ion_exchange_sites");
      double cec = GetAttributeValueD_(static_cast<DOMElement*>(node), "cec");

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
          ss << "dof " << j + 1 << " function";

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
      std::string name;
      std::vector<DOMNode*> sites = GetSameChildNodes_(node, name, flag);

      if (flag) {
        std::stringstream site_str;
        site_str << "MESH BLOCK " << i+1;
       
        sorption_sites.clear();

        Teuchos::ParameterList& complexation = ic_list.sublist("sorption_sites");
        Teuchos::ParameterList& tmp_list = complexation.sublist("function")
            .sublist(site_str.str())
            .set<Teuchos::Array<std::string> >("regions", regions)
            .set<std::string>("component", "cell")
            .sublist("function")
            .set<int>("number of dofs", sites.size())
            .set<std::string>("function type", "composite function");

        for (int k = 0; k < sites.size(); ++k) {
          element = static_cast<DOMElement*>(sites[k]);
          double val = GetAttributeValueD_(element, "density", TYPE_NUMERICAL);
          sorption_sites.push_back(GetAttributeValueS_(element, "name", TYPE_NONE));

          std::stringstream dof_str;
          dof_str << "dof " << k+1 << " function";
          tmp_list.sublist(dof_str.str()).sublist("function-constant")
                                         .set<double>("value", val);
        }
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
          ss << "dof " << j + 1 << " function";

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

  out_list.sublist("verbose object") = verb_list_.sublist("verbose object");
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


/* ******************************************************************
* Adds kd values to a bgd file (using the first material).
* Returns the name of this file.  
****************************************************************** */
std::string InputConverterU::CreateINFile(std::string& filename)
{
  MemoryManager mm;
  DOMNode* node;
  DOMElement* element;

  std::ofstream in_file;
  std::stringstream primaries;
  std::stringstream secondaries;
  std::stringstream minerals;
  std::stringstream mineral_list;
  std::stringstream gases;
  std::stringstream mineral_kinetics;
  std::stringstream isotherms;
  std::stringstream cations;
  std::stringstream complexes;
  std::stringstream constraints;
  std::stringstream reactionrates;
  std::stringstream decayrates;

  // database filename
  bool flag;
  node = GetUniqueElementByTagsString_("process_kernels, chemistry", flag);
  element = static_cast<DOMElement*>(node);
  std::string datfilename = GetAttributeValueS_(element, "database", TYPE_NONE, true, "");

  // get species names
  int nsolutes = phases_["water"].size();
  for (int i = 0; i < nsolutes; ++i) {
    std::string name = phases_["water"][i];
    primaries << "    " << name << "\n";
  }

  // check for forward/backward rates on primaries (for solutes/non-reactive primaries only)
  node = GetUniqueElementByTagsString_("liquid_phase, dissolved_components, primaries", flag);
  if (flag) {
    std::string primary;
    std::vector<DOMNode*> children = GetSameChildNodes_(node, primary, flag, false);

    for (int i = 0; i < children.size(); ++i) {
      DOMNode* inode = children[i];
      element = static_cast<DOMElement*>(inode);
      if (element->hasAttribute(mm.transcode("forward_rate"))) {
        double frate = GetAttributeValueD_(element, "forward_rate");
        double brate = GetAttributeValueD_(element, "backward_rate");
        std::string name = TrimString_(mm.transcode(inode->getTextContent()));
        reactionrates << "    REACTION " << name << " <->\n";
        reactionrates << "    FORWARD_RATE " << frate << "\n";
        reactionrates << "    BACKWARD_RATE " << brate << "\n";
      }
      if (element->hasAttribute(mm.transcode("first_order_decay_constant"))) {
        double decay = GetAttributeValueD_(element, "first_order_decay_constant");
        std::string name = TrimString_(mm.transcode(inode->getTextContent()));
        decayrates << "    REACTION " << name << " <-> \n";
        decayrates << "    RATE_CONSTANT " << decay << "\n";
      }
    }
  }


  node = GetUniqueElementByTagsString_("liquid_phase, dissolved_components, secondaries", flag);
  if (flag) {
    //std::vector<DOMNode*> children = static_cast<DOMElement*>(node)->getElementsByTagName(mm.transcode("secondary"));
    std::string secondary;
    std::vector<DOMNode*> children = GetSameChildNodes_(node, secondary, flag, false);

    for (int i = 0; i < children.size(); ++i) {
      DOMNode* inode = children[i];
      std::string name = TrimString_(mm.transcode(inode->getTextContent()));
      secondaries << "    " << name << "\n";
    }
  }


  node = GetUniqueElementByTagsString_("phases, solid_phase, minerals", flag);
  if (flag) {
    //std::vector<DOMNode*> children = static_cast<DOMElement*>(node)->getElementsByTagName(mm.transcode("mineral"));
    std::string mineral;
    std::vector<DOMNode*> children = GetSameChildNodes_(node, mineral, flag, false);

    for (int i = 0; i < children.size(); ++i) {
      DOMNode* inode = children[i];
      std::string name = TrimString_(mm.transcode(inode->getTextContent()));
      minerals << "    " << name << "\n";
      mineral_list << name << ", ";

      element = static_cast<DOMElement*>(inode);
      double rate = GetAttributeValueD_(element, "rate_constant", TYPE_NUMERICAL, false, 0.0);
      if (rate > 0.0) {
        mineral_kinetics << "    " << name << "\n";
        mineral_kinetics << "      RATE_CONSTANT " << rate << "\n";
        //mineral_kinetics << "      RATE_CONSTANT " << rate << " mol/cm^2-sec\n";
        mineral_kinetics << "    /\n";
      }
    }
  }

  node = GetUniqueElementByTagsString_("phases, gas_phase, gases", flag);
  if (flag) {
    //std::vector<DOMNode*> children = static_cast<DOMElement*>(node)->getElementsByTagName(mm.transcode("gas"));
    std::string gas;
    std::vector<DOMNode*> children = GetSameChildNodes_(node, gas, flag, false);

    for (int i = 0; i < children.size(); ++i) {
      DOMNode* inode = children[i];
      std::string name = TrimString_(mm.transcode(inode->getTextContent()));
      gases << "    " << name << "\n";
    }
  }

  // soprtion isotherms
  node = GetUniqueElementByTagsString_("materials, material, sorption_isotherms", flag);
  if (flag) {
    std::string name;
    std::vector<DOMNode*> children = GetSameChildNodes_(node, name, flag, false);
    for (int i = 0; i < children.size(); ++i) {
      DOMNode* inode = children[i];
      element = static_cast<DOMElement*>(inode);
      name = GetAttributeValueS_(element, "name");

      DOMNode* knode = GetUniqueElementByTagsString_(inode, "kd_model", flag);
      element = static_cast<DOMElement*>(knode);
      if (flag) {
        double kd = GetAttributeValueD_(element, "kd");
        std::string model = GetAttributeValueS_(element, "model");
        if (model == "linear") {
	  isotherms << "      " << name << "\n";
	  isotherms << "        TYPE LINEAR\n";
	  isotherms << "        DISTRIBUTION_COEFFICIENT " << kd << "\n";
	} else if (model == "langmuir") {
	  isotherms << "      " << name << "\n";
	  isotherms << "        TYPE LANGMUIR\n";
	  isotherms << "        DISTRIBUTION_COEFFICIENT " << kd << "\n";
          double b = GetAttributeValueD_(element, "b");
          isotherms << "        LANGMUIR_B " << b << "\n";
        } else if (model == "freundlich") {
	  isotherms << "      " << name << "\n";
	  isotherms << "        TYPE FREUNDLICH\n";
	  isotherms << "        DISTRIBUTION_COEFFICIENT " << kd << "\n";
          double n = GetAttributeValueD_(element, "n");
          isotherms << "        FREUNDLICH_N " << n << "\n";
        }
	isotherms << "      /\n";
      }
    }
  }
  
  // soprtion ion exchange
  node = GetUniqueElementByTagsString_("materials, material, ion_exchange, cations", flag);
  if (flag) {
    element = static_cast<DOMElement*>(node);
    double cec = GetAttributeValueD_(element, "cec");
    cations << "      CEC " << cec << "\n";
    cations << "      CATIONS\n";

    std::string name;
    std::vector<DOMNode*> children = GetSameChildNodes_(node, name, flag, false);
    for (int i = 0; i < children.size(); ++i) {
      DOMNode* inode = children[i];
      element = static_cast<DOMElement*>(inode);
      name = GetAttributeValueS_(element, "name");
      double value = GetAttributeValueD_(element, "value");
      cations << "        " << name << " " << value << "\n";
    }
  }

  // soprtion surface complexation
  //node = GetUniqueElementByTagsString_("materials, material, surface_complexation, mineral", flag);
  //if (flag) {
  //  std::string name = GetTextContentS_(node, mineral_list.str().c_str());
  //  complexes << "      MINERAL " << name << "\n";
  //}
  node = GetUniqueElementByTagsString_("materials, material, surface_complexation", flag);
  if (flag) {
    std::string name;
    std::vector<DOMNode*> children = GetSameChildNodes_(node, name, flag, false);
    
    for (int i = 0; i < children.size(); ++i) {
      DOMNode* inode = children[i];
      element = static_cast<DOMElement*>(inode);
      name = GetAttributeValueS_(element, "name", TYPE_NONE);
      double density = GetAttributeValueD_(element, "density");
      complexes << "    SURFACE_COMPLEXATION_RXN\n";
      complexes << "      EQUILIBRIUM\n";
      complexes << "      SITE " << name << " " << density << "\n";
      
      //DOMNode* cnode = GetUniqueElementByTagsString_("materials, material, surface_complexation, complexes", flag);
      std::vector<DOMNode*> kids = GetChildren_(inode, "complexes", flag);
      if (flag) {
        complexes << "      COMPLEXES\n";
        for (int j = 0; j < kids.size(); ++j) {
          DOMNode* jnode = kids[j];
          std::vector<std::string> complexe_names = CharToStrings_(mm.transcode(jnode->getTextContent()));
          for (std::vector<std::string>::const_iterator it = complexe_names.begin(); it != complexe_names.end(); it++) {
            complexes << "        " << *it << "\n";
          }
        }
        complexes << "      /\n";
        complexes << "    /\n";
      }
    }
  }

  // constraints
  node = GetUniqueElementByTagsString_("geochemistry, constraints", flag);
  if (flag) {

    // loop over constraints
    std::string name;
    std::vector<DOMNode*> constraint_list = GetSameChildNodes_(node, name, flag, false);
    for (int i = 0; i < constraint_list.size(); ++i) {
      std::stringstream primary;
      std::stringstream mineral;
      std::stringstream gas;

      DOMNode* inode = constraint_list[i];
      element = static_cast<DOMElement*>(inode);
      std::string const_name = GetAttributeValueS_(element, "name");
       
      // loop over individual items in current constraint
      std::string constraint;
      bool found;
      std::vector<DOMNode*> children = GetChildren_(inode, "primary", found);
      if (found) {
        for (int j = 0; j < children.size(); ++j) {
          DOMNode* jnode = children[j];
          element = static_cast<DOMElement*>(jnode);

	  std::string constname = GetAttributeValueS_(element, "name");
	  std::string ctype = GetAttributeValueS_(element, "type");
	  double constvalue = GetAttributeValueD_(element, "value");
	  std::string typeletter;
	  if (ctype == "free_ion") {
	      typeletter = "F";
	  } else if (ctype == "pH" || ctype == "ph") {
	      typeletter = "P";
	  } else if (ctype == "total") {
	      typeletter = "T";
	  } else if (ctype == "total+sorbed") {
	      typeletter = "S";
	  } else if (ctype == "charge") {
	      typeletter = "Z";
	  } else if (ctype == "mineral") {
	      std::string mname = GetAttributeValueS_(element, "mineral");
	      typeletter = "M " + mname;
	  } else if (ctype == "gas") {
	      std::string gname = GetAttributeValueS_(element, "gas");
	      typeletter = "G " + gname;
	  }
	  primary << "    " << constname << " " << constvalue << " " << typeletter << "\n";
	}
      }

      children = GetChildren_(inode, "mineral", found);
      if (found) {
        for (int j = 0; j < children.size(); ++j) {
          DOMNode* jnode = children[j];
          element = static_cast<DOMElement*>(jnode);
	  std::string constname = GetAttributeValueS_(element, "name");
	  double constvf = GetAttributeValueD_(element, "volume_fraction");
	  double constsa = GetAttributeValueD_(element, "specific_surface_area");
          mineral << "    " << constname << " " << constvf << " " << constsa << "\n";
          //mineral << "    " << constname << " " << constvf << " " << constsa << " cm^2/cm^3\n";
	}
      }

      children = GetChildren_(inode, "gas", found);
      if (found) {
        for (int j = 0; j < children.size(); ++j) {
          DOMNode* jnode = children[j];
          element = static_cast<DOMElement*>(jnode);
	  std::string constname = GetAttributeValueS_(element, "name");
	  gas << "    " << constname << "\n";
	}
      }

      constraints << "CONSTRAINT " << const_name << "\n";
      constraints << "  CONCENTRATIONS\n";
      constraints << primary.str();
      constraints << "  /\n";
      if (!mineral.str().empty()) {
        constraints << "  MINERALS\n";
        constraints << mineral.str();
        constraints << "  /\n";
      }
      constraints << "END\n";
    }
  }


  // build in filename
  std::string infilename(filename);
  std::string new_extension(".in");
  size_t pos = infilename.find(".xml");
  infilename.replace(pos, (size_t)4, new_extension, (size_t)0, (size_t)4);

  // open output in file
  if (rank_ == 0) {
    struct stat buffer;
    int status = stat(infilename.c_str(), &buffer);
    if (!status) {
      if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH) {
        Teuchos::OSTab tab = vo_->getOSTab();
        *vo_->os() << "File \"" << infilename.c_str() 
                   << "\" exists, skipping its creation." << std::endl;
      }
    } else {
      in_file.open(infilename.c_str());

      in_file << "CHEMISTRY\n";
      // Chemistry - Species
      in_file << "  PRIMARY_SPECIES\n";
      in_file << primaries.str();
      in_file << "  /\n";
      if (!reactionrates.str().empty()) {
        in_file << "  GENERAL_REACTION\n";
        in_file << reactionrates.str();
        in_file << "  /\n";
      }
      if (!decayrates.str().empty()) {
        in_file << "  RADIOACTIVE_DECAY_REACTION\n";
        in_file << decayrates.str();
        in_file << "  /\n";
      }
      if (!secondaries.str().empty()) {
        in_file << "  SECONDARY_SPECIES\n";
        in_file << secondaries.str();
        in_file << "  /\n";
      }
      //in_file << "  REDOX_SPECIES\n";
      //in_file << species.str();
      //in_file << "  /\n";
      if (!gases.str().empty()) {
        in_file << "  GAS_SPECIES\n";
        in_file << gases.str();
        in_file << "  /\n";
      }
      if (!minerals.str().empty()) {
        in_file << "  MINERALS\n";
        in_file << minerals.str();
        in_file << "  /\n";
      }

      // Chemistry - Mineral Kinetics
      if (!mineral_kinetics.str().empty()) {
        in_file << "  MINERAL_KINETICS\n";
        in_file << mineral_kinetics.str();
        in_file << "  /\n";
      }

      // Chemistry - Sorption
      if (!isotherms.str().empty() || !complexes.str().empty() || !cations.str().empty()) {
        in_file << "  SORPTION\n";
        if (!isotherms.str().empty()) {
          in_file << "    ISOTHERM_REACTIONS\n";
          in_file << isotherms.str();
          in_file << "    /\n";
        }
        if (!complexes.str().empty()) {
          in_file << complexes.str();
        }
        if (!cations.str().empty()) {
          in_file << "    ION_EXCHANGE_RXN\n";
          in_file << cations.str();
          in_file << "      /\n";
          in_file << "    /\n";
        }
        in_file << "  /\n";
      }

      // Chemistry - Controls
      in_file << "  DATABASE " << datfilename.c_str() << "\n";
      in_file << "  LOG_FORMULATION\n";
      in_file << "  ACTIVITY_COEFFICIENTS TIMESTEP\n";
      in_file << "END\n";

      // Constraints
      in_file << constraints.str();
      //in_file << "END\n";

      in_file.close();
    }
  }

  return infilename;
}
}  // namespace AmanziInput
}  // namespace Amanzi


