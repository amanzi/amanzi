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
	    std::string inpfilename = CreateINFile_(xmlfilename_, rank_);
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

        Teuchos::ParameterList& aux3_list = freundlich_n.sublist("function").sublist(*it)
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
  }

  // general parameters
  int max_itrs(100), cut_threshold(8), increase_threshold(4);
  double tol(1e-12), dt_max(1e+10), dt_min(1e+10), dt_init(1e+7), dt_cut(2.0), dt_increase(1.2);

  double ion_value;
  bool ion_guess(false);

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
      } else if (strcmp(text, "free_ion_guess") == 0) {
        ion_guess = true;
        ion_value = strtod(mm.transcode(inode->getTextContent()), NULL);
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

  // free ion has optional initialization. Default initialization is tight
  // to the valus of the initial value of the total component concentration.
  if (ion_guess) {
    int nsolutes = phases_["water"].size();
    Teuchos::ParameterList& free_ion = ic_list.sublist("free_ion_species");

    Teuchos::ParameterList& aux1_list = free_ion.sublist("function").sublist("All")
        .set<std::string>("region", "All")
        .set<std::string>("component", "cell")
        .sublist("function");
    aux1_list.set<int>("number of dofs", nsolutes).set("function type", "composite function");

    for (int j = 0; j < nsolutes; ++j) {
      std::string solute_name = phases_["water"][j];

      std::stringstream ss;
      ss << "dof " << j + 1 << " function";

      aux1_list.sublist(ss.str()).sublist("function-constant").set<double>("value", ion_value);
    }
  }

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
  MemoryManager mm;

  // get and write list of primaries/solutes
  bool flag;
  node = GetUniqueElementByTagsString_("phases, liquid_phase, dissolved_components, primaries", flag);
  if (flag) {
    std::string name;
    bool flag2;
    std::vector<DOMNode*> children = GetSameChildNodes_(node, name, flag, false);
    for (int i = 0; i < children.size(); ++i) {
      DOMNode* inode = children[i];
      name = mm.transcode(inode->getTextContent());
      species << name << " ;   0.00 ;   0.00 ;   1.00 \n";
    }
  } else {
    node = GetUniqueElementByTagsString_("phases, liquid_phase, dissolved_components, solutes", flag);

    if (flag) {
      std::string name;
      bool flag2;
      std::vector<DOMNode*> children = GetSameChildNodes_(node, name, flag, false);
      for (int i = 0; i < children.size(); ++i) {
        DOMNode* inode = children[i];
        name = mm.transcode(inode->getTextContent());
        species << name << " ;   0.00 ;   0.00 ;   1.00 \n";
      }
    }
  }
  
  // loop over materials to get kd information
  node = GetUniqueElementByTagsString_("materials", flag);
  if (flag ) {
    
    // get kd information from all materials
    Teuchos::ParameterList IsothermsPL;
    std::string name;
    std::vector<DOMNode*> mat_list = GetSameChildNodes_(node, name, flag, false);
    for (int i = 0; i < mat_list.size(); ++i) {
      bool flag2;
      DOMElement* child_elem;
      
      DOMNode* inode = mat_list[i];
      element = static_cast<DOMElement*>(inode);
      name = GetAttributeValueS_(element, "name");
      
      // look for sorption_isotherms
      child_elem = GetChildByName_(inode, "sorption_isotherms", flag2, false);
      if (flag2) {
        
        // loop over sublist of primaries to get Kd information
        std::vector<DOMNode*> primary_list = GetSameChildNodes_(static_cast<DOMNode*>(child_elem), name, flag2, false);
        if (flag2) {
          
          for (int j = 0; j < primary_list.size(); ++j) {
            DOMNode* jnode = primary_list[j];
            std::string primary_name = GetAttributeValueS_(static_cast<DOMElement*>(jnode), "name");
            DOMNode* kd_node = GetUniqueElementByTagsString_(jnode, "kd_model", flag2);
            DOMElement* kd_elem = static_cast<DOMElement*>(kd_node);
            
            if (flag2) {
              Teuchos::ParameterList kd_list;
              double kd = GetAttributeValueD_(kd_elem, "kd");
              std::string model = GetAttributeValueS_(kd_elem, "model");
              kd_list.set<std::string>("model",model);
              kd_list.set<double>("kd",kd);
              
              if (model == "langmuir") {
                double b = GetAttributeValueD_(kd_elem, "b");
                kd_list.set<double>("b",b);
              } else if (model == "freundlich") {
                double n = GetAttributeValueD_(kd_elem, "n");
                kd_list.set<double>("n",n);
              }
              
              IsothermsPL.sublist(primary_name) = kd_list;
            }
          }
        }
      }
    }
    
    // create text for kds
    for (Teuchos::ParameterList::ConstIterator iter = IsothermsPL.begin(); iter != IsothermsPL.end(); ++iter) {
      
      std::string primary = IsothermsPL.name(iter);
      Teuchos::ParameterList& curprimary = IsothermsPL.sublist(primary);
      std::string model = curprimary.get<std::string>("model");
      double kd = curprimary.get<double>("kd");
      
      if (model == "langmuir") {
        double b = curprimary.get<double>("b");
        isotherms << primary << " ; langmuir ; " << kd << " " << b << std::endl;
      } else if (model == "freundlich") {
        double n = curprimary.get<double>("n");
        isotherms << primary << " ; freundlich ; " << kd << " " << n << std::endl;
      } else {
        isotherms << primary << " ; linear ; " << kd << std::endl;
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


