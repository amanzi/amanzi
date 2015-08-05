/*
  This is the input component of the Amanzi code. 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <sstream>
#include <fstream>
#include <string>
#include <algorithm>
#include <boost/lambda/lambda.hpp>
#include <boost/bind.hpp>
#include <boost/algorithm/string.hpp>

#include "errors.hh"
#include "exceptions.hh"
#include "dbc.hh"

#define  BOOST_FILESYTEM_NO_DEPRECATED
#include "boost/filesystem/operations.hpp"
#include "boost/filesystem/path.hpp"
#include "boost/format.hpp"
#include "boost/lexical_cast.hpp"

// TPLs
#include <xercesc/dom/DOM.hpp>
#include <xercesc/util/XMLString.hpp>
#include <xercesc/parsers/DOMLSParserImpl.hpp>

#include "Teuchos_XMLParameterListHelpers.hpp"

// Amanzi's
#include "ErrorHandler.hpp"
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

  if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH) {
    *vo_->os() << "Translating output" << std::endl;
  }

  DOMNamedNodeMap* attr_map;
  DOMNode* tmp_node;
  DOMNode* node_attr;
  char* text_content;
  char* text_content2;

  // get definitions node - this node MAY exist ONCE
  // this contains any time macros and cycle macros
  // they are stored in the outputs of the old format
  Teuchos::ParameterList tmPL, cmPL;
  DOMNodeList* macro_list = doc_->getElementsByTagName(XMLString::transcode("macros"));

  if (macro_list->getLength() > 0) {
    DOMNodeList* children = macro_list->item(0)->getChildNodes();

    for (int i = 0; i < children->getLength(); i++) {
      DOMNode* inode = children->item(i) ;
      if (DOMNode::ELEMENT_NODE == inode->getNodeType()) {
        char* tagname = XMLString::transcode(inode->getNodeName());

        // process time macros
        if (strcmp(tagname, "time_macro") == 0) {
          Teuchos::ParameterList tm_parameter;
          attr_map = inode->getAttributes();
          node_attr = attr_map->getNamedItem(XMLString::transcode("name"));
          if (node_attr) {
            text_content = XMLString::transcode(node_attr->getNodeValue());
          } else {
            ThrowErrorMissattr_("definitions", "attribute", "name", "time_macro");
          }

          // deal differently if "times" or "start-inter-stop"
          DOMNodeList* children = inode->getChildNodes();
          bool isTime(false);

          for (int j = 0; j < children->getLength(); j++) {
            DOMNode* time_node = children->item(j) ;
            if (DOMNode::ELEMENT_NODE == time_node->getNodeType()) {
              if (strcmp(XMLString::transcode(time_node->getNodeName()), "time") == 0) isTime = true;
            }   
          }

          if (isTime) {
            Teuchos::Array<double> times;
            for (int j = 0; j < children->getLength(); j++) {
              DOMNode* time_node = children->item(j) ;
              if (DOMNode::ELEMENT_NODE == time_node->getNodeType()) {
                char* node_txt = XMLString::transcode(time_node->getTextContent());
                times.append(GetTimeValue_(node_txt));
                XMLString::release(&node_txt);
              }
            }
            tm_parameter.set<Teuchos::Array<double> >("values", times);

          } else {
            DOMElement* element = static_cast<DOMElement*>(inode);
            DOMNodeList* list = element->getElementsByTagName(XMLString::transcode("start"));
            tmp_node = list->item(0);

            char* node_txt = XMLString::transcode(tmp_node->getTextContent());
            Teuchos::Array<double> sps;
            sps.append(GetTimeValue_(node_txt));
            XMLString::release(&node_txt);
            list = element->getElementsByTagName(XMLString::transcode("timestep_interval"));

            if (list->getLength() > 0) {
              tmp_node = list->item(0);
              node_txt = XMLString::transcode(tmp_node->getTextContent());
              sps.append(GetTimeValue_(node_txt));
              XMLString::release(&node_txt);

              list = element->getElementsByTagName(XMLString::transcode("stop"));
              if (list->getLength() > 0) {
                tmp_node = list->item(0);
                node_txt = XMLString::transcode(tmp_node->getTextContent());
                sps.append(GetTimeValue_(node_txt));
                XMLString::release(&node_txt);
              } else {
                sps.append(-1.0);
              }

              tm_parameter.set<Teuchos::Array<double> >("sps", sps);
            } else {
              tm_parameter.set<Teuchos::Array<double> >("values", sps);
            }
          }
          tmPL.sublist(text_content) = tm_parameter;
          XMLString::release(&text_content);
        }

        // process cycle macros
        else if (strcmp(tagname,"cycle_macro") == 0) {
          Teuchos::ParameterList cm_parameter;
          attr_map = inode->getAttributes();
          node_attr = attr_map->getNamedItem(XMLString::transcode("name"));
          if (node_attr) {
            text_content = XMLString::transcode(node_attr->getNodeValue());
          } else {
            ThrowErrorMissattr_("definitions", "attribute", "name", "cycle_macro");
          }

          DOMElement* element = static_cast<DOMElement*>(inode);
          DOMNodeList* list = element->getElementsByTagName(XMLString::transcode("start"));
          tmp_node = list->item(0);

          char* node_txt = XMLString::transcode(tmp_node->getTextContent());
          Teuchos::Array<int> sps;
          sps.append(atoi(node_txt));
          XMLString::release(&node_txt);

          list = element->getElementsByTagName(XMLString::transcode("timestep_interval"));
          if (list->getLength() > 0) {
            tmp_node = list->item(0);
            node_txt = XMLString::transcode(tmp_node->getTextContent());
            sps.append(atoi(node_txt));
            XMLString::release(&node_txt);

            list = element->getElementsByTagName(XMLString::transcode("stop"));
            if (list->getLength() > 0) {
              tmp_node = list->item(0);
              node_txt = XMLString::transcode(tmp_node->getTextContent());
              sps.append(atoi(node_txt));
              XMLString::release(&node_txt);
            } else {
              sps.append(-1);
            }
            cm_parameter.set<Teuchos::Array<int> >("sps", sps);
          } else {
            cm_parameter.set<Teuchos::Array<int> >("values", sps);
          }
          cmPL.sublist(text_content) = cm_parameter;
          XMLString::release(&text_content);
        }
      }
    }
  }

  // get output node - this node must exist ONCE
  DOMNodeList* node_list = doc_->getElementsByTagName(XMLString::transcode("output"));
  DOMNode* out_node = node_list->item(0);

  if (DOMNode::ELEMENT_NODE == out_node->getNodeType()) {
    DOMNodeList* out_child_list = out_node->getChildNodes();
    for (int m = 0; m < out_child_list->getLength(); m++) {
      DOMNode* mnode = out_child_list->item(m);

      if (DOMNode::ELEMENT_NODE == mnode->getNodeType()) {
        char* outName = XMLString::transcode(mnode->getNodeName());
        if (strcmp(outName, "vis") == 0) {
          if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH) {
            *vo_->os() << "Translating visualization" << std::endl;
          }

          // get list of vis - this node MAY exist ONCE
          DOMNodeList* children = mnode->getChildNodes();
          Teuchos::ParameterList visPL;
          for (int j = 0; j < children->getLength(); j++) {
            DOMNode* jnode = children->item(j);
            text_content = XMLString::transcode(jnode->getNodeName());

            if (strcmp(text_content, "base_filename") == 0) {
              text_content2 = XMLString::transcode(jnode->getTextContent());
              visPL.set<std::string>("file name base", TrimString_(text_content2));
              XMLString::release(&text_content2);
            }
            else if (strcmp(text_content, "num_digits") == 0) {
              text_content2 = XMLString::transcode(jnode->getTextContent());
              visPL.set<int>("file name digits", atoi(text_content2));
              XMLString::release(&text_content2);
            }
            else if (strcmp(text_content, "cycle_macros") == 0 ||
                     strcmp(text_content, "cycle_macro") == 0) {
              text_content2 = XMLString::transcode(jnode->getTextContent());
              ProcessMacros_("cycles", text_content2, cmPL, visPL);
              XMLString::release(&text_content2);
            }
            else if (strcmp(text_content, "time_macros") == 0 ||
                     strcmp(text_content, "time_macro") == 0) {
              text_content2 = XMLString::transcode(jnode->getTextContent());
              ProcessMacros_("times", text_content2, tmPL, visPL);
              XMLString::release(&text_content2);
            }
            else if (strcmp(text_content, "write_regions") == 0) {
              text_content2 = XMLString::transcode(jnode->getTextContent());
              visPL.set<Teuchos::Array<std::string> >("write regions", CharToStrings_(text_content2));
              XMLString::release(&text_content2);
            }
            XMLString::release(&text_content);
          }
          out_list.sublist("Visualization Data") = visPL;
        }

        else if (strcmp(outName, "checkpoint") == 0) {
          if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH) {
            *vo_->os() << "Translating checkpoint" << std::endl;
          }
          // get list of checkpoint - this node MAY exist ONCE
          Teuchos::ParameterList chkPL;
          DOMNodeList* children = mnode->getChildNodes();
          for (int j = 0; j < children->getLength(); j++) {
            DOMNode* jnode = children->item(j) ;
            text_content = XMLString::transcode(jnode->getNodeName());

            if (strcmp(text_content, "base_filename") == 0) {
              text_content2 = XMLString::transcode(jnode->getTextContent());
              chkPL.set<std::string>("file name base", TrimString_(text_content2));
              XMLString::release(&text_content2);
            }
            else if (strcmp(text_content, "num_digits") == 0) {
              text_content2 = XMLString::transcode(jnode->getTextContent());
              chkPL.set<int>("file name digits", atoi(text_content2));
              XMLString::release(&text_content2);
            }
            else if (strcmp(text_content, "cycle_macros") == 0 ||
                     strcmp(text_content, "cycle_macro") == 0) {
              text_content2 = XMLString::transcode(jnode->getTextContent());
              ProcessMacros_("cycles", text_content2, cmPL, chkPL);
              XMLString::release(&text_content2);
            }
            XMLString::release(&text_content);
          }
          out_list.sublist("Checkpoint Data") = chkPL;
        }

        else if (strcmp(outName,"walkabout") == 0) {
          if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH) {
            *vo_->os() << "Translating walkabout" << std::endl;
          }

          // get list of walkabout - this node MAY exist ONCE
          Teuchos::ParameterList chkPL;
          DOMNodeList* children = mnode->getChildNodes();
          for (int j = 0; j < children->getLength(); j++) {
            DOMNode* jnode = children->item(j);
            text_content = XMLString::transcode(jnode->getNodeName());
            if (strcmp(text_content,"base_filename") == 0) {
              text_content2 = XMLString::transcode(jnode->getTextContent());
              chkPL.set<std::string>("file name base", TrimString_(text_content2));
              XMLString::release(&text_content2);
            }
            else if (strcmp(text_content, "num_digits") == 0) {
              text_content2 = XMLString::transcode(jnode->getTextContent());
              chkPL.set<int>("file name digits", atoi(text_content2));
              XMLString::release(&text_content2);
            }
            else if (strcmp(text_content, "cycle_macros") == 0 ||
                     strcmp(text_content, "cycle_macro") == 0) {
              text_content2 = XMLString::transcode(jnode->getTextContent());
              ProcessMacros_("cycles", text_content2, cmPL, chkPL);
              XMLString::release(&text_content2);
            }
            XMLString::release(&text_content);
          }
          out_list.sublist("Walkabout Data") = chkPL;
        }

        else if (strcmp(outName,"observations") == 0) {
          if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH) {
            *vo_->os() << "Translating observations" << std::endl;
          }

          Teuchos::ParameterList obsPL;
          DOMNodeList* OBList = mnode->getChildNodes();
          for (int i = 0; i < OBList->getLength(); i++) {
            DOMNode* inode = OBList->item(i) ;
            if (DOMNode::ELEMENT_NODE == inode->getNodeType()) {
              text_content  = XMLString::transcode(inode->getNodeName());
              if (strcmp(text_content, "filename") == 0) {
                text_content2 = XMLString::transcode(inode->getTextContent());
                obsPL.set<std::string>("Observation Output Filename", TrimString_(text_content2));
                XMLString::release(&text_content2);
              }
              else if (strcmp(text_content, "liquid_phase") == 0) {
                attr_map = inode->getAttributes();
                node_attr = attr_map->getNamedItem(XMLString::transcode("name"));
                std::string phaseName;

                if (node_attr) {
                  text_content2 = XMLString::transcode(node_attr->getNodeValue());
                  phaseName = std::string(text_content2);
                  if (phaseName=="water") {
                    phaseName = "Water";
                  }
                  XMLString::release(&text_content2);
                }
                else {
                  ThrowErrorMissattr_("observations", "attribute", "name", "liquid_phase");
                }

                // loop over observations
                DOMNodeList* children = inode->getChildNodes();
                for (int j = 0; j < children->getLength(); j++) {
                  Teuchos::ParameterList obPL;
                  DOMNode* jnode = children->item(j);
                  if (DOMNode::ELEMENT_NODE == jnode->getNodeType()) {
                    char* obsType = XMLString::transcode(jnode->getNodeName());
                    if (strcmp(obsType, "aqueous_pressure") == 0) {
                      obPL.set<std::string>("variable", "Aqueous pressure");
                    }
                    else if (strcmp(obsType, "integrated_mass") == 0) {
                      // TODO: can't find matching version
                    }
                    else if (strcmp(obsType, "volumetric_water_content") == 0) {
                      obPL.set<std::string>("variable", "Volumetric water content");
                    }
                    else if (strcmp(obsType, "gravimetric_water_content") == 0) {
                      obPL.set<std::string>("variable", "Gravimetric water content");
                    }
                    else if (strcmp(obsType, "x_aqueous_volumetric_flux") == 0) {
                      obPL.set<std::string>("variable", "X-Aqueous volumetric flux");
                    }
                    else if (strcmp(obsType, "y_aqueous_volumetric_flux") == 0) {
                      obPL.set<std::string>("variable", "Y-Aqueous volumetric flux");
                    }
                    else if (strcmp(obsType, "z_aqueous_volumetric_flux") == 0) {
                      obPL.set<std::string>("variable", "Z-Aqueous volumetric flux");
                    }
                    else if (strcmp(obsType, "material_id") == 0) {
                      obPL.set<std::string>("variable", "MaterialID");
                    }
                    else if (strcmp(obsType, "hydraulic_head") == 0) {
                      obPL.set<std::string>("variable", "Hydraulic Head");
                    }
                    else if (strcmp(obsType, "aqueous_mass_flow_rate") == 0) {
                      obPL.set<std::string>("variable", "Aqueous mass flow rate");
                    }
                    else if (strcmp(obsType, "aqueous_volumetric_flow_rate") == 0) {
                      obPL.set<std::string>("variable", "Aqueous volumetric flow rate");
                    }
                    else if (strcmp(obsType, "aqueous_saturation") == 0) {
                      obPL.set<std::string>("variable", "Aqueous saturation");
                    }
                    else if (strcmp(obsType, "aqueous_conc") == 0) {
                      // get solute name
                      attr_map = jnode->getAttributes();
                      node_attr = attr_map->getNamedItem(XMLString::transcode("solute"));
                      char* soluteName;
                      if (node_attr) {
                        soluteName = XMLString::transcode(node_attr->getNodeValue());
                      } else {
                        ThrowErrorMissattr_("observations", "attribute", "solute", "aqueous_conc");
                      }

                      std::stringstream name;
                      name<< soluteName << " Aqueous concentration";
                      obPL.set<std::string>("variable", name.str());
                    }
                    else if (strcmp(obsType, "drawdown") == 0) {
                      obPL.set<std::string>("variable", "Drawdown");
                    }
                    else if (strcmp(obsType, "solute_volumetric_flow_rate") == 0) {
                      // get solute name
                      attr_map = jnode->getAttributes();
                      node_attr = attr_map->getNamedItem(XMLString::transcode("solute"));
                      char* soluteName;
                      if (node_attr) {
                        soluteName = XMLString::transcode(node_attr->getNodeValue());
                      } else {
                        ThrowErrorMissattr_("observations", "attribute", "solute", "solute_volumetric_flow_rate");
                      }
                      
                      std::stringstream name;
                      name<< soluteName << " volumetric flow rate";
                      obPL.set<std::string>("variable", name.str());
                    }

                    DOMNodeList* kidList = jnode->getChildNodes();
                    for (int k = 0; k < kidList->getLength(); k++) {
                      DOMNode* knode = kidList->item(k) ;
                      if (DOMNode::ELEMENT_NODE == knode->getNodeType()) {
                        char* elem = XMLString::transcode(knode->getNodeName());
                        char* value = XMLString::transcode(knode->getTextContent());

                        if (strcmp(elem, "assigned_regions") == 0) {
                          //TODO: EIB - really a note, REGION != ASSIGNED REGIONS, this isn't consistent!!!
                          obPL.set<std::string>("region", TrimString_(value));
                        }
                        else if (strcmp(elem, "functional") == 0) {
                          if (strcmp(value, "point") == 0) {
                            obPL.set<std::string>("functional", "Observation Data: Point");
                          }
                          else if (strcmp(elem, "integral") == 0) {
                            obPL.set<std::string>("functional", "Observation Data: Integral");
                          }
                          else if (strcmp(elem, "mean") == 0) {
                            obPL.set<std::string>("functional", "Observation Data: Mean");
                          }
                        }

                        // Keeping singular macro around to help users. This will go away
                        else if (strcmp(elem, "time_macros") == 0 ||
                                 strcmp(elem, "time_macro") == 0) {
                          ProcessMacros_("times", value, tmPL, obPL);
                        }
                        else if (strcmp(elem, "cycle_macros") == 0 ||
                                 strcmp(elem, "cycle_macro") == 0) {
                          ProcessMacros_("cycles", value, cmPL, obPL);
                        }
                        XMLString::release(&elem);
                        XMLString::release(&value);
                      }
                    }
                    std::stringstream list_name;
                    list_name << "observation-" << j+1 << ":" << phaseName;
                    obsPL.sublist(list_name.str()) = obPL;
                  }
                }
              }
              XMLString::release(&text_content);
              out_list.sublist("Observation Data") = obsPL;
            }
          }
        }
        XMLString::release(&outName);
      }
    }
  }

  return out_list;
}


/* ******************************************************************
* Converts macros and macro to multiple SPS parameters orthe single
* values parameter.
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

      outPL.set<Teuchos::Array<double> >("values", times);
    }
  }
}

}  // namespace AmanziInput
}  // namespace Amanzi
