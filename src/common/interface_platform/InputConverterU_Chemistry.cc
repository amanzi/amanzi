/*
  Input Converter

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
           Sergi Molins <smolins@lbl.gov>
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
#include "Key.hh"

#include "InputConverterU.hh"
#include "InputConverterU_Defs.hh"

namespace Amanzi {
namespace AmanziInput {

XERCES_CPP_NAMESPACE_USE

/* ******************************************************************
* Create chemistry list.
****************************************************************** */
Teuchos::ParameterList InputConverterU::TranslateChemistry_(const std::string& domain)
{
  Teuchos::ParameterList out_list;

  if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH)
    *vo_->os() << "Translating chemistry, domain=" << domain << std::endl;

  MemoryManager mm;
  char* text;
  DOMNodeList *node_list, *children;
  DOMNode* node;
  DOMElement* element;

  // create header
  out_list.set<std::string>("domain name", (domain == "matrix") ? "domain" : domain);

  // chemical engine
  bool flag;
  node = GetPKChemistryPointer_(flag);
  std::string engine = GetAttributeValueS_(node, "engine");
  std::string bgdfilename = GetAttributeValueS_(node, "input_filename", TYPE_NONE, false, "");

  // process engine
  if (engine ==  "amanzi") {
    out_list.set<std::string>("chemistry model", "Amanzi");

    node = GetUniqueElementByTagsString_("process_kernels, chemistry", flag);

    if (bgdfilename != "") {
      auto pair = Keys::split(bgdfilename, '.');
      if (pair.second != "xml") {
        Errors::Message msg;
        msg << "Incorect suffix for optional XML file with a thermodynamic database.\n";
        Exceptions::amanzi_throw(msg);
      }
      Teuchos::ParameterList& bgd_list = out_list.sublist("thermodynamic database");
      bgd_list.set<std::string>("file", bgdfilename);
      if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH)
        *vo_->os() << " using file:" << bgdfilename << std::endl;
    }
  } else {
    bool valid_engine(true);

    if (engine == "pflotran") {
      out_list.set<std::string>("engine", "PFloTran");
    } else if (engine == "crunchflow") {
      out_list.set<std::string>("engine", "CrunchFlow");
    } else {
      valid_engine = false;
    }

    // Pass along chemistry engine info.
    if (valid_engine) {
      out_list.set<std::string>("chemistry model", "Alquimia");

      // Find the name of the engine-specific input file.
      node = GetPKChemistryPointer_(flag);
      if (flag) {
        std::string inpfilename;
        element = static_cast<DOMElement*>(node);
	if (element->hasAttribute(mm.transcode("input_filename"))) {
          inpfilename = GetAttributeValueS_(element, "input_filename");
          out_list.set<std::string>("engine input file", inpfilename);
	} else {
	  inpfilename = CreateINFile_(xmlfilename_, rank_);
          out_list.set<std::string>("engine input file", inpfilename);
	} 

        if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH)
          *vo_->os() << " using file:" << inpfilename << std::endl;
      } else {
        Errors::Message msg;
        msg << "Chemical process kernel has not been found.\n";
        Exceptions::amanzi_throw(msg);
      }
    }
  }
  
  // minerals
  std::vector<std::string> minerals;
  std::vector<double> min_rate_cnst;
  node = GetUniqueElementByTagsString_("phases, solid_phase, minerals", flag);
  if (flag) {
    children = static_cast<DOMElement*>(node)->getElementsByTagName(mm.transcode("mineral"));

    for (int i = 0; i < children->getLength(); ++i) {
      DOMNode* inode = children->item(i);

      double mrc(0.0);
      mrc = GetAttributeValueD_(inode, "rate_constant", TYPE_NUMERICAL, DVAL_MIN, DVAL_MAX, "", false, 0.0);
      min_rate_cnst.push_back(mrc);

      std::string name = TrimString_(mm.transcode(inode->getTextContent()));
      minerals.push_back(name);
    }
  }

  // first-order decay rate
  std::vector<std::string> aqueous_reactions;
  std::vector<double> kin_rate_cnst;
  node = GetUniqueElementByTagsString_("phases, liquid_phase, dissolved_components, primaries", flag);
  if (flag) {
    children = static_cast<DOMElement*>(node)->getElementsByTagName(mm.transcode("primary"));

    for (int i = 0; i < children->getLength(); ++i) {
      DOMNode* inode = children->item(i);

      double krate(-99.9);
      krate = GetAttributeValueD_(inode, "first_order_decay_constant", TYPE_NUMERICAL, DVAL_MIN, DVAL_MAX, "", false, -99.9);

      if (krate != -99.9) {
	kin_rate_cnst.push_back(krate);
	std::string name = TrimString_(mm.transcode(inode->getTextContent()));
	aqueous_reactions.push_back(name);
      }
    }
  }

  // region specific initial conditions
  std::vector<std::string> sorption_sites;
  Teuchos::ParameterList& ic_list = out_list.sublist("initial conditions");

  if (domain == "fracture") {
    node = GetUniqueElementByTagsString_("fracture_network, materials", flag);
    element = static_cast<DOMElement*>(node);
  } else {
    node_list = doc_->getElementsByTagName(mm.transcode("materials"));
    element = static_cast<DOMElement*>(node_list->item(0));
  }

  children = element->getElementsByTagName(mm.transcode("material"));

  for (int i = 0; i < children->getLength(); ++i) {
    DOMNode* inode = children->item(i);

    node = GetUniqueElementByTagsString_(inode, "assigned_regions", flag);
    text = mm.transcode(node->getTextContent());
    std::vector<std::string> regions = CharToStrings_(text);

    // aqueous kinetic reactions (only first order, really)
    if (aqueous_reactions.size() > 0) {
      out_list.set<Teuchos::Array<std::string> >("aqueous_reactions", aqueous_reactions);

      Teuchos::ParameterList& rate = ic_list.sublist("first_order_decay_constant");

      for (std::vector<std::string>::const_iterator it = regions.begin(); it != regions.end(); it++) {
        Teuchos::ParameterList& aux3_list = rate.sublist("function").sublist(*it)
            .set<std::string>("region", *it)
            .set<std::string>("component", "cell")
            .sublist("function");
        aux3_list.set<int>("number of dofs", aqueous_reactions.size())
            .set("function type", "composite function");

        for (int j = 0; j < aqueous_reactions.size(); ++j) {
          std::stringstream ss;
          ss << "dof " << j + 1 << " function";
 
          node = GetUniqueElementByTagsString_(inode, "aqueous_reactions", flag);
          double arc(0.0);
          if (flag) {
            element = GetUniqueChildByAttribute_(node, "name", aqueous_reactions[j], flag, true);
	    arc = GetAttributeValueD_(element, "first_order_rate_constant", TYPE_NUMERICAL, DVAL_MIN, DVAL_MAX, "", false, 0.0);
          }
          aux3_list.sublist(ss.str()).sublist("function-constant").set<double>("value", kin_rate_cnst[j]);
        }
      }
    }
    
    // mineral volume fraction and specific surface area.
    if (minerals.size() > 0) {
      out_list.set<Teuchos::Array<std::string> >("minerals", minerals);

      // if (pk_model_["chemistry"] == "amanzi") {
      Teuchos::ParameterList& volfrac = ic_list.sublist("mineral_volume_fractions");
      Teuchos::ParameterList& surfarea = ic_list.sublist("mineral_specific_surface_area");
      Teuchos::ParameterList& rate = ic_list.sublist("mineral_rate_constant");

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

        Teuchos::ParameterList& aux3_list = rate.sublist("function").sublist(*it)
            .set<std::string>("region", *it)
            .set<std::string>("component", "cell")
            .sublist("function");
        aux3_list.set<int>("number of dofs", minerals.size())
            .set("function type", "composite function");

        for (int j = 0; j < minerals.size(); ++j) {
          std::stringstream ss;
          ss << "dof " << j + 1 << " function";
 
          node = GetUniqueElementByTagsString_(inode, "minerals", flag);
          double mvf(0.0), msa(0.0);
          if (flag) {
            element = GetUniqueChildByAttribute_(node, "name", minerals[j], flag, true);
            mvf = GetAttributeValueD_(element, "volume_fraction", TYPE_NUMERICAL, 0.0, 1.0, "", false, 0.0);
            msa = GetAttributeValueD_(element, "specific_surface_area", TYPE_NUMERICAL, 0.0, DVAL_MAX, "", false, 0.0);
	    // mrc = GetAttributeValueD_(element, "rate_constant", TYPE_NUMERICAL, 0.0, DVAL_MAX, "", false, 0.0);
          }
          aux1_list.sublist(ss.str()).sublist("function-constant").set<double>("value", mvf);
          aux2_list.sublist(ss.str()).sublist("function-constant").set<double>("value", msa);
          aux3_list.sublist(ss.str()).sublist("function-constant").set<double>("value", min_rate_cnst[j]);
        }
      }
    }

    // ion exchange
    node = GetUniqueElementByTagsString_(inode, "ion_exchange, cations", flag);
    if (flag) {
      Teuchos::ParameterList& sites = ic_list.sublist("ion_exchange_sites");
      double cec = GetAttributeValueD_(node, "cec", TYPE_NUMERICAL, DVAL_MIN, DVAL_MAX, "mol/m^3");

      for (auto it = regions.begin(); it != regions.end(); ++it) {
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

          double val = (!flag) ? 0.0 : GetAttributeValueD_(element, "kd", TYPE_NUMERICAL, 0.0, DVAL_MAX, "", false, 0.0);
          aux1_list.sublist(ss.str()).sublist("function-constant").set<double>("value", val);

          val = (!flag) ? 0.0 : GetAttributeValueD_(element, "b", TYPE_NUMERICAL, 0.0, DVAL_MAX, "", false, 0.0);
          aux2_list.sublist(ss.str()).sublist("function-constant").set<double>("value", val);

          val = (!flag) ? 0.0 : GetAttributeValueD_(element, "n", TYPE_NUMERICAL, 0.0, DVAL_MAX, "", false, 0.0);
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
          double val = GetAttributeValueD_(element, "density", TYPE_NUMERICAL, 0.0, DVAL_MAX);
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

  std::string activity_model("unit"), dt_method("fixed"), pitzer_database;
  std::vector<std::string> aux_data;

  node = GetUniqueElementByTagsString_("unstructured_controls, unstr_chemistry_controls", flag);
  if (flag) {
    children = node->getChildNodes();
    for (int i = 0; i < children->getLength(); ++i) {
      DOMNode* inode = children->item(i);
      if (inode->getNodeType() != DOMNode::ELEMENT_NODE) continue;

      text = mm.transcode(inode->getNodeName());
      if (strcmp(text, "activity_model") == 0) {
        activity_model = GetTextContentS_(inode, "unit, debye-huckel, pitzer-hwm");
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
      } else if (strcmp(text, "pitzer_database") == 0) {
        pitzer_database = mm.transcode(inode->getTextContent());
      }
    }
  }
  out_list.set<std::string>("activity model", activity_model);
  if (pitzer_database.size() > 0)
    out_list.set<std::string>("Pitzer database file", pitzer_database);
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
  out_list.set<double>("initial conditions time", ic_time_);
  out_list.set<int>("number of component concentrations", comp_names_all_.size());

  out_list.sublist("verbose object") = verb_list_.sublist("verbose object");
  return out_list;
}


/* ******************************************************************
* Helper utility for two stuctures of process_kernel lists
****************************************************************** */
DOMNode* InputConverterU::GetPKChemistryPointer_(bool& flag)
{
  MemoryManager mm;
  DOMNode* node = NULL;
  DOMNodeList *children;

  node = GetUniqueElementByTagsString_("process_kernels, chemistry", flag);
  if (!flag) {
    node = GetUniqueElementByTagsString_("process_kernels", flag);
    children = static_cast<DOMElement*>(node)->getElementsByTagName(mm.transcode("pk"));
    for (int i = 0; i < children->getLength(); ++i) {
      node = GetUniqueElementByTagsString_(children->item(i), "chemistry", flag);
      if (flag) break;
    }
  }

  return node;
}

}  // namespace AmanziInput
}  // namespace Amanzi


