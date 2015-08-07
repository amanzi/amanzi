/*
  This is the input component of the Amanzi code. 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <string>

// TPLs
#include <boost/algorithm/string.hpp>

#include "Teuchos_ParameterList.hpp"

#include "InputConverterU.hh"

namespace Amanzi {
namespace AmanziInput {

XERCES_CPP_NAMESPACE_USE

/* ******************************************************************
* Empty
****************************************************************** */
Teuchos::ParameterList InputConverterU::Translate()
{
  Teuchos::ParameterList out_list;
  
  // grab verbosity early
  verb_list_ = TranslateVerbosity_();
  Teuchos::ParameterList tmp_list(verb_list_);
  vo_ = new VerboseObject("InputConverter", tmp_list);
  Teuchos::OSTab tab = vo_->getOSTab();

  // parsing of miscalleneous lists
  ParseSolutes_();
  
  out_list.sublist("Mesh") = TranslateMesh_();
  out_list.sublist("Domain").set<int>("Spatial Dimension", dim_);
  out_list.sublist("Regions") = TranslateRegions_();
  out_list.sublist("Output") = TranslateOutput_();
  out_list.sublist("Solvers") = TranslateSolvers_();
  out_list.sublist("Preconditioners") = TranslatePreconditioners_();
  out_list.sublist("State") = TranslateState_();
  out_list.sublist("Cycle Driver") = TranslateCycleDriver_();

  Teuchos::ParameterList& cd_list = out_list.sublist("Cycle Driver");
  out_list.sublist("PKs") = TranslatePKs_(cd_list);
  
  return out_list;
}


/* ******************************************************************
* Extract information of solute components.
****************************************************************** */
void InputConverterU::ParseSolutes_()
{
  bool flag;
  char* tagname;
  char* text_content;

  DOMNode* inode = doc_->getElementsByTagName(XMLString::transcode("phases"))->item(0);
  DOMNode* node = getUniqueElementByTagsString_(inode, "liquid_phase, dissolved_components, solutes", flag);

  DOMNodeList* children = node->getChildNodes();
  for (int i = 0; i < children->getLength(); ++i) {
    inode = children->item(i);
    tagname = XMLString::transcode(inode->getNodeName());
    text_content = XMLString::transcode(inode->getTextContent());

    if (strcmp(tagname, "solute") == 0) {
      phases_["water"].push_back(TrimString_(text_content));
    }

    XMLString::release(&text_content);
    XMLString::release(&tagname);
  }
}


/* ******************************************************************
* Extract generic verbosity object for all sublists.
****************************************************************** */
Teuchos::ParameterList InputConverterU::TranslateVerbosity_()
{
  Teuchos::ParameterList vlist;

  DOMNodeList* node_list;
  DOMNode* node_attr;
  DOMNamedNodeMap* attr_map;
  char* text_content;
    
  // get execution contorls node
  XMLCh* xstr = XMLString::transcode("execution_controls");
  node_list = doc_->getElementsByTagName(xstr);
  XMLString::release(&xstr);
  
  for (int i = 0; i < node_list->getLength(); i++) {
    DOMNode* inode = node_list->item(i);

    if (DOMNode::ELEMENT_NODE == inode->getNodeType()) {
      DOMNodeList* children = inode->getChildNodes();

      for (int j = 0; j < children->getLength(); j++) {
        DOMNode* jnode = children->item(j);

        if (DOMNode::ELEMENT_NODE == jnode->getNodeType()) {
          char* tagname = XMLString::transcode(jnode->getNodeName());
          if (std::string(tagname) == "verbosity") {
            attr_map = jnode->getAttributes();
            node_attr = attr_map->getNamedItem(XMLString::transcode("level"));
            if (node_attr) {
              text_content = XMLString::transcode(node_attr->getNodeValue());
              vlist.sublist("VerboseObject").set<std::string>("Verbosity Level", TrimString_(text_content));
              XMLString::release(&text_content);
              break;
            } else {
              ThrowErrorIllformed_("verbosity", "value", "level");
            }
            XMLString::release(&text_content);
          }
        }
      }
    }
  }
  return vlist;
}

}  // namespace AmanziInput
}  // namespace Amanzi
