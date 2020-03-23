/*
  Input Converter

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Erin Barker (original version)
           Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <sstream>
#include <fstream>
#include <string>
#include <algorithm>

#define  BOOST_FILESYTEM_NO_DEPRECATED

// TPLs
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "xercesc/dom/DOM.hpp"
#include "xercesc/parsers/DOMLSParserImpl.hpp"

// Amanzi's
#include "ErrorHandler.hpp"
#include "errors.hh"
#include "exceptions.hh"
#include "dbc.hh"

#include "InputConverterU.hh"

namespace Amanzi {
namespace AmanziInput {

XERCES_CPP_NAMESPACE_USE

/* ******************************************************************
* Empty
****************************************************************** */
Teuchos::ParameterList InputConverterU::TranslateOutput_()
{
  Teuchos::ParameterList out_list;

  if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH)
    *vo_->os() << "Translating output" << std::endl;
  Teuchos::OSTab tab = vo_->getOSTab();

  MemoryManager mm;

  char *tagname, *text;
  std::string unit;

  DOMNamedNodeMap* attr_map;
  DOMNodeList *node_list;
  DOMNode* node;

  // get definitions node - this node MAY exist ONCE
  // this contains any time macros and cycle macros
  Teuchos::ParameterList tmPL, cmPL;
  node_list = doc_->getElementsByTagName(mm.transcode("macros"));
  DOMNode* inode = node_list->item(0)->getFirstChild();

  while (inode != NULL) {
    if (DOMNode::ELEMENT_NODE == inode->getNodeType()) {
      tagname = mm.transcode(inode->getNodeName());

      // process time macros
      if (strcmp(tagname, "time_macro") == 0) {
        Teuchos::ParameterList tm_parameter;
        std::string name = GetAttributeValueS_(inode, "name");

        // deal differently if "times" or "start-period-stop"
        bool flag;
        std::string child_name;
        std::vector<DOMNode*> multi_list = GetSameChildNodes_(inode, child_name, flag, false);

        if (child_name == "time") {
          std::vector<double> times;
          for (int j = 0; j < multi_list.size(); j++) {
            DOMNode* jnode = multi_list[j];
            if (DOMNode::ELEMENT_NODE == jnode->getNodeType()) {
              text = mm.transcode(jnode->getTextContent());
              times.push_back(ConvertUnits_(TrimString_(text), unit));
            }
          }
          tm_parameter.set<Teuchos::Array<double> >("values", times);
        } else {
          DOMElement* element = static_cast<DOMElement*>(inode);
          DOMNodeList* list = element->getElementsByTagName(mm.transcode("start"));
          node = list->item(0);

          text = mm.transcode(node->getTextContent());
          Teuchos::Array<double> sps;
          sps.append(TimeCharToValue_(text));

          list = element->getElementsByTagName(mm.transcode("timestep_interval"));
          if (list->getLength() > 0) {
            text = mm.transcode(list->item(0)->getTextContent());
            sps.append(TimeCharToValue_(text));

            list = element->getElementsByTagName(mm.transcode("stop"));
            if (list->getLength() > 0) {
              text = mm.transcode(list->item(0)->getTextContent());
              sps.append(TimeCharToValue_(text));
            } else {
              sps.append(-1.0);
            }
            tm_parameter.set<Teuchos::Array<double> >("sps", sps);
          } else {
            tm_parameter.set<Teuchos::Array<double> >("values", sps);
          }
        }
        tmPL.sublist(name) = tm_parameter;
      }

      // process cycle macros
      else if (strcmp(tagname,"cycle_macro") == 0) {
        Teuchos::ParameterList cm_parameter;
        std::string name = GetAttributeValueS_(inode, "name");      

        DOMElement* element = static_cast<DOMElement*>(inode);
        DOMNodeList* list = element->getElementsByTagName(mm.transcode("start"));
        node = list->item(0);

        text = mm.transcode(node->getTextContent());
        Teuchos::Array<int> sps;
        sps.append(std::strtol(text, NULL, 10));

        list = element->getElementsByTagName(mm.transcode("timestep_interval"));
        if (list->getLength() > 0) {
          text = mm.transcode(list->item(0)->getTextContent());
          sps.append(std::strtol(text, NULL, 10));

          list = element->getElementsByTagName(mm.transcode("stop"));
          if (list->getLength() > 0) {
            text = mm.transcode(list->item(0)->getTextContent());
            sps.append(std::strtol(text, NULL, 10));
          } else {
            sps.append(-1);
          }
          cm_parameter.set<Teuchos::Array<int> >("sps", sps);
        } else {
          cm_parameter.set<Teuchos::Array<int> >("values", sps);
        }
        cmPL.sublist(name) = cm_parameter;
      }
    }
    inode = inode->getNextSibling();
  }
  if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH)
    *vo_->os() << "Found " << cmPL.numParams() << " cycle macros and "
               << tmPL.numParams() << " time macros." << std::endl;

  // get output->vis node - this node must exist ONCE
  bool flag;
  node = GetUniqueElementByTagsString_("output, vis", flag);

  if (flag && node->getNodeType() == DOMNode::ELEMENT_NODE) {
    if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH) {
      *vo_->os() << "Translating output: visualization" << std::endl;
    }

    DOMNodeList* children = node->getChildNodes();
    int nchildren = children->getLength();

    Teuchos::ParameterList visPL;
    for (int j = 0; j < nchildren; j++) {
      DOMNode* jnode = children->item(j);
      tagname = mm.transcode(jnode->getNodeName());

      if (strcmp(tagname, "base_filename") == 0) {
        text = mm.transcode(jnode->getTextContent());
        visPL.set<std::string>("file name base", output_prefix_ + TrimString_(text));
      } else if (strcmp(tagname, "file_format") == 0) {
        text = mm.transcode(jnode->getTextContent());
        visPL.set<std::string>("file format", TrimString_(text));
      } else if (strcmp(tagname, "cycle_macros") == 0 ||
                 strcmp(tagname, "cycle_macro") == 0) {
        text = mm.transcode(jnode->getTextContent());
        ProcessMacros_("cycles", text, cmPL, visPL);
      } else if (strcmp(tagname, "time_macros") == 0 ||
                 strcmp(tagname, "time_macro") == 0) {
        text = mm.transcode(jnode->getTextContent());
        ProcessMacros_("times", text, tmPL, visPL);
      } else if (strcmp(tagname, "write_regions") == 0) {
        DOMNodeList* kids = jnode->getChildNodes();
        for (int k = 0; k < kids->getLength(); ++k) { 
          DOMNode* knode = kids->item(k);
          if (knode->getNodeType() != DOMNode::ELEMENT_NODE) continue;

          std::string name = GetAttributeValueS_(knode, "name");
          std::string regs = GetAttributeValueS_(knode, "regions");
          std::vector<std::string> regions = CharToStrings_(regs.c_str());
          visPL.sublist("write regions").set<Teuchos::Array<std::string> >(name, regions);
        }
      } else if (strcmp(tagname, "write_partition") == 0) {
        text = mm.transcode(jnode->getTextContent());
        if (strcmp(text, "true") == 0)
          visPL.set<bool>("write partitions", true);
      }
    }
    out_list.sublist("visualization data") = visPL;
  }

  // get output->checkpoint node - this node must exist ONCE
  node = GetUniqueElementByTagsString_("output, checkpoint", flag);

  if (flag && node->getNodeType() == DOMNode::ELEMENT_NODE) {
    if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH)
      *vo_->os() << "Translating output: checkpoint" << std::endl;

    Teuchos::ParameterList chkPL;
    DOMNodeList* children = node->getChildNodes();
    for (int j = 0; j < children->getLength(); j++) {
      DOMNode* jnode = children->item(j) ;
      tagname = mm.transcode(jnode->getNodeName());
      text = mm.transcode(jnode->getTextContent());

      if (strcmp(tagname, "base_filename") == 0) {
        chkPL.set<std::string>("file name base", output_prefix_ + TrimString_(text));
      }
      else if (strcmp(tagname, "num_digits") == 0) {
        chkPL.set<int>("file name digits", std::strtol(text, NULL, 10));
      } else if (strcmp(tagname, "cycle_macros") == 0 ||
                 strcmp(tagname, "cycle_macro") == 0) {
        ProcessMacros_("cycles", text, cmPL, chkPL);
      } else if (strcmp(tagname, "time_macros") == 0 ||
                 strcmp(tagname, "time_macro") == 0) {
        ProcessMacros_("times", text, tmPL, chkPL);
      }
    }
    out_list.sublist("checkpoint data") = chkPL;
  }

  // get output->mesh info - this node must exist ONCE
  node = GetUniqueElementByTagsString_("output, mesh_info", flag);

  if (flag && node->getNodeType() == DOMNode::ELEMENT_NODE) {
    io_mesh_info_ = true;
    if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH)
      *vo_->os() << "Translating output: additional mesh information" << std::endl;

    Teuchos::ParameterList chkPL;
    DOMNodeList* children = node->getChildNodes();
    for (int j = 0; j < children->getLength(); j++) {
      DOMNode* jnode = children->item(j) ;
      tagname = mm.transcode(jnode->getNodeName());
      text = mm.transcode(jnode->getTextContent());

      if (strcmp(tagname, "filename") == 0) {
        chkPL.set<std::string>("filename", output_prefix_ + TrimString_(text));
      }
    }
    out_list.sublist("mesh info") = chkPL;
  }

  // get output->walkabout node - this node must exist ONCE
  node = GetUniqueElementByTagsString_("output, walkabout", flag);

  if (flag && node->getNodeType() == DOMNode::ELEMENT_NODE) {
    io_walkabout_ = true;
    if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH)
      *vo_->os() << "Translating output: walkabout" << std::endl;

    Teuchos::ParameterList chkPL;
    DOMNodeList* children = node->getChildNodes();
    for (int j = 0; j < children->getLength(); j++) {
      DOMNode* jnode = children->item(j);
      tagname = mm.transcode(jnode->getNodeName());
      text = mm.transcode(jnode->getTextContent());

      if (strcmp(tagname, "base_filename") == 0) {
        chkPL.set<std::string>("file name base", output_prefix_ + TrimString_(text));
      } else if (strcmp(tagname, "num_digits") == 0) {
        chkPL.set<int>("file name digits", std::strtol(text, NULL, 10));
      } else if (strcmp(tagname, "cycle_macros") == 0 ||
                 strcmp(tagname, "cycle_macro") == 0) {
        ProcessMacros_("cycles", text, cmPL, chkPL);
      } else if (strcmp(tagname, "time_macros") == 0 ||
                 strcmp(tagname, "time_macro") == 0) {
        ProcessMacros_("times", text, tmPL, chkPL);
      }
    }

    chkPL.sublist("write regions")
        .set<Teuchos::Array<std::string> >("region names", material_regions_)
        .set<Teuchos::Array<std::string> >("material names", material_names_)
        .set<Teuchos::Array<int> >("material ids", material_ids_);

    out_list.sublist("walkabout data") = chkPL;
  }

  // get output->observations node - only one node is allowed
  int nobs_liquid(0), nobs_gas(0);
  node = GetUniqueElementByTagsString_("output, observations", flag);

  if (flag && node->getNodeType() == DOMNode::ELEMENT_NODE) {
    if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH)
      *vo_->os() << "Translating output: observations" << std::endl;

    Teuchos::ParameterList obsPL;
    inode = node->getFirstChild();
    while (inode != NULL) {
      if (inode->getNodeType() != DOMNode::ELEMENT_NODE) { 
        inode = inode->getNextSibling();
        continue;
      }
      tagname = mm.transcode(inode->getNodeName());
      text = mm.transcode(inode->getTextContent());

      if (strcmp(tagname, "filename") == 0) {
        obsPL.set<std::string>("observation output filename", output_prefix_ + TrimString_(text));

      } else if (strcmp(tagname, "units") == 0) {
        unit = GetAttributeValueS_(inode, "time", TYPE_NONE, false, units_.system().time);
        obsPL.set<std::string>("time unit", unit);

        unit = GetAttributeValueS_(inode, "mass", TYPE_NONE, false, units_.system().mass);
        obsPL.set<std::string>("mass unit", unit);

        unit = GetAttributeValueS_(inode, "length", TYPE_NONE, false, units_.system().length);
        obsPL.set<std::string>("length unit", unit);

        unit = GetAttributeValueS_(inode, "concentration", TYPE_NONE, false, units_.system().concentration);
        obsPL.set<std::string>("concentration unit", unit);

      } else if (strcmp(tagname, "liquid_phase") == 0) {
        node = inode->getAttributes()->getNamedItem(mm.transcode("name"));
        if (!node)
            ThrowErrorMissing_("observations", "attribute", "name", "liquid_phase");

        nobs_liquid = 0;
        DOMNodeList* children = inode->getChildNodes();
        for (int j = 0; j < children->getLength(); j++) {
          Teuchos::ParameterList obPL;
          DOMNode* jnode = children->item(j);
          if (jnode->getNodeType() != DOMNode::ELEMENT_NODE) continue;
          char* obs_type = mm.transcode(jnode->getNodeName());

          if (strcmp(obs_type, "aqueous_pressure") == 0) {
            obPL.set<std::string>("variable", "aqueous pressure");
          } else if (strcmp(obs_type, "volumetric_water_content") == 0) {
            obPL.set<std::string>("variable", "volumetric water content");
          } else if (strcmp(obs_type, "gravimetric_water_content") == 0) {
            obPL.set<std::string>("variable", "gravimetric water content");
          } else if (strcmp(obs_type, "x_aqueous_volumetric_flux") == 0) {
            obPL.set<std::string>("variable", "x-aqueous volumetric flux");
          } else if (strcmp(obs_type, "y_aqueous_volumetric_flux") == 0) {
            obPL.set<std::string>("variable", "y-aqueous volumetric flux");
          } else if (strcmp(obs_type, "z_aqueous_volumetric_flux") == 0) {
            obPL.set<std::string>("variable", "z-aqueous volumetric flux");
          } else if (strcmp(obs_type, "material_id") == 0) {
            obPL.set<std::string>("variable", "material id");
          } else if (strcmp(obs_type, "hydraulic_head") == 0) {
            obPL.set<std::string>("variable", "hydraulic head");
          } else if (strcmp(obs_type, "perm_weighted_hydraulic_head") == 0) {
            obPL.set<std::string>("variable", "permeability-weighted hydraulic head");
          } else if (strcmp(obs_type, "aqueous_mass_flow_rate") == 0) {
            obPL.set<std::string>("variable", "aqueous mass flow rate");
          } else if (strcmp(obs_type, "aqueous_volumetric_flow_rate") == 0) {
            obPL.set<std::string>("variable", "aqueous volumetric flow rate");
          } else if (strcmp(obs_type, "aqueous_saturation") == 0) {
            obPL.set<std::string>("variable", "aqueous saturation");
          } else if (strcmp(obs_type, "ph") == 0) {
            obPL.set<std::string>("variable", "pH");
          } else if (strcmp(obs_type, "aqueous_conc") == 0) {
            std::string solute_name = GetAttributeValueS_(jnode, "solute");
            std::stringstream name;
            name << solute_name << " aqueous concentration";
            obPL.set<std::string>("variable", name.str());
          } else if (strcmp(obs_type, "drawdown") == 0) {
            obPL.set<std::string>("variable", "drawdown");
          } else if (strcmp(obs_type, "perm_weighted_drawdown") == 0) {
            obPL.set<std::string>("variable", "permeability-weighted drawdown");
          } else if (strcmp(obs_type, "water_table") == 0) {
            obPL.set<std::string>("variable", "water table");
          } else if (strcmp(obs_type, "solute_volumetric_flow_rate") == 0) {
            std::string solute_name = GetAttributeValueS_(jnode, "solute");
            std::stringstream name;
            name << solute_name << " volumetric flow rate";
            obPL.set<std::string>("variable", name.str());
          } else if (strcmp(obs_type, "fractures_aqueous_volumetric_flow_rate") == 0) {
            obPL.set<std::string>("variable", "fractures aqueous volumetric flow rate");
          }

          std::vector<std::string> regions;

          DOMNodeList* kids = jnode->getChildNodes();
          for (int k = 0; k < kids->getLength(); k++) {
            DOMNode* knode = kids->item(k);
            if (DOMNode::ELEMENT_NODE == knode->getNodeType()) {
              char* elem = mm.transcode(knode->getNodeName());
              char* value = mm.transcode(knode->getTextContent());

              if (strcmp(elem, "assigned_regions") == 0) {
                regions = CharToStrings_(value);
              } else if (strcmp(elem, "functional") == 0) {
                if (strcmp(value, "point") == 0) {
                  obPL.set<std::string>("functional", "observation data: point");
                } else if (strcmp(value, "integral") == 0) {
                  obPL.set<std::string>("functional", "observation data: integral");
                } else if (strcmp(value, "mean") == 0) {
                  obPL.set<std::string>("functional", "observation data: mean");
                }
              } else if (strcmp(elem, "domain_name") == 0) {
                obPL.set<std::string>("domain name", (strcmp(value, "matrix") == 0) ? "domain" : "fracture");
              } else if (strcmp(elem, "interpolation") == 0) {
                if (strcmp(value, "constant") == 0) {
                  obPL.set<std::string>("interpolation", "constant");
                } else if (strcmp(value, "linear")==0) {
                  obPL.set<std::string>("interpolation", "linear");
                }
              } else if (strcmp(elem, "weighting") == 0) {
                if (strcmp(value, "none") == 0) {
                  obPL.set<std::string>("weighting", "none");
                } else if (strcmp(value, "flux_norm")==0) {
                  obPL.set<std::string>("weighting", "flux norm");
                }
              // Keeping singular macro around to help users. This will go away
              } else if (strcmp(elem, "time_macros") == 0 ||
                         strcmp(elem, "time_macro") == 0) {
                ProcessMacros_("times", value, tmPL, obPL);
              } else if (strcmp(elem, "cycle_macros") == 0 ||
                         strcmp(elem, "cycle_macro") == 0) {
                ProcessMacros_("cycles", value, cmPL, obPL);
              }
            }
          }

          // process list of regions
          for (int k = 0; k < regions.size(); ++k) {
            obPL.set<std::string>("region", regions[k]);
            vv_obs_regions_.push_back(regions[k]);

            std::stringstream list_name;
            list_name << "obl " << ++nobs_liquid;
            obsPL.sublist(list_name.str()) = obPL;
          }
        }

      } else if (strcmp(tagname, "gas_phase") == 0) {
        node = inode->getAttributes()->getNamedItem(mm.transcode("name"));
        if (!node)
            ThrowErrorMissing_("observations", "attribute", "name", "gas_phase");

        nobs_gas = 0;
        DOMNodeList* children = inode->getChildNodes();
        for (int j = 0; j < children->getLength(); j++) {
          Teuchos::ParameterList obPL;
          DOMNode* jnode = children->item(j);
          if (jnode->getNodeType() != DOMNode::ELEMENT_NODE) continue;
          char* obs_type = mm.transcode(jnode->getNodeName());

          if (strcmp(obs_type, "gaseous_conc") == 0) {
            std::string solute_name = GetAttributeValueS_(jnode, "solute");
            std::stringstream name;
            name << solute_name << " gaseous concentration";
            obPL.set<std::string>("variable", name.str());
          }

          DOMNodeList* kids = jnode->getChildNodes();
          for (int k = 0; k < kids->getLength(); k++) {
            DOMNode* knode = kids->item(k);
            if (DOMNode::ELEMENT_NODE == knode->getNodeType()) {
              char* elem = mm.transcode(knode->getNodeName());
              char* value = mm.transcode(knode->getTextContent());

              if (strcmp(elem, "assigned_regions") == 0) {
                obPL.set<std::string>("region", TrimString_(value));
                vv_obs_regions_.push_back(TrimString_(value));
              } else if (strcmp(elem, "functional") == 0) {
                if (strcmp(value, "point") == 0) {
                  obPL.set<std::string>("functional", "observation data: point");
                } else if (strcmp(value, "integral") == 0) {
                  obPL.set<std::string>("functional", "observation data: integral");
                } else if (strcmp(value, "mean") == 0) {
                  obPL.set<std::string>("functional", "observation data: mean");
                }
              // Keeping singular macro around to help users. This will go away
              } else if (strcmp(elem, "time_macros") == 0 ||
                         strcmp(elem, "time_macro") == 0) {
                ProcessMacros_("times", value, tmPL, obPL);
              } else if (strcmp(elem, "cycle_macros") == 0 ||
                         strcmp(elem, "cycle_macro") == 0) {
                ProcessMacros_("cycles", value, cmPL, obPL);
              }
            }
          }
          std::stringstream list_name;
          list_name << "obg " << ++nobs_gas;
          obsPL.sublist(list_name.str()) = obPL;
        }
      }
      inode = inode->getNextSibling();
    }

    out_list.sublist("observation data") = obsPL;
  }

  if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH) {
    *vo_->os() << "Found " << nobs_liquid << " liquid observations" << std::endl;
    *vo_->os() << "Found " << nobs_gas << " gas observations" << std::endl;
  }

  return out_list;
}


/* ******************************************************************
* Converts macros and macro to multiple SPS parameters or the single
* value parameter.
****************************************************************** */
void InputConverterU::ProcessMacros_(
    const std::string& prefix, char* text_content,
    Teuchos::ParameterList& mPL, Teuchos::ParameterList& outPL)
{
  std::vector<std::string> macro = CharToStrings_(text_content);
  Teuchos::Array<int> cycles;
  Teuchos::Array<double> times;

  std::list<int> cm_list;
  std::list<double> tm_list;

  int k(0);
  bool flag(false);
  for (int i = 0; i < macro.size(); i++) {
    if (!mPL.isSublist(macro[i])) {
      Errors::Message msg;
      msg << "Macro \"" << macro[i] << "\" was not defined.\n";
      msg << "Please correct and try again.\n";
      Exceptions::amanzi_throw(msg);
    }

    Teuchos::ParameterList& mlist = mPL.sublist(macro[i]);
    if (mlist.isParameter("sps")) {
      std::stringstream ss;
      ss << prefix << " start period stop " << k;
      k++;
      if (prefix == "cycles") {
        outPL.set<Teuchos::Array<int> >(ss.str(), mlist.get<Teuchos::Array<int> >("sps"));
      } else {
        outPL.set<Teuchos::Array<double> >(ss.str(), mlist.get<Teuchos::Array<double> >("sps"));
      }
    } else if (mlist.isParameter("values")) {
      flag = true;
      if (prefix == "cycles") {
        cycles = mlist.get<Teuchos::Array<int> >("values");
        for (int n = 0; n < cycles.size(); n++) cm_list.push_back(cycles[n]);
      } else {
        times = mlist.get<Teuchos::Array<double> >("values");
        for (int n = 0; n < times.size(); n++) tm_list.push_back(times[n]);
      }
    }
  }

  if (flag) {
    if (prefix == "cycles") {
      cm_list.sort();
      cm_list.unique();

      cycles.clear();
      for (std::list<int>::iterator it = cm_list.begin(); it != cm_list.end(); ++it) 
          cycles.push_back(*it);

      outPL.set<Teuchos::Array<int> >("values", cycles);
    } else {
      tm_list.sort();
      tm_list.unique();

      times.clear();
      for (std::list<double>::iterator it = tm_list.begin(); it != tm_list.end(); ++it) 
          times.push_back(*it);

      outPL.set<Teuchos::Array<double> >("times", times);
    }
  }
}

}  // namespace AmanziInput
}  // namespace Amanzi
