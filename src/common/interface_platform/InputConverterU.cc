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
  verb_list_ = GetVerbosity_();
  vo_ = new VerboseObject("InputTranslator", verb_list_);
  Teuchos::OSTab tab = vo_->getOSTab();
  
  out_list.sublist("Mesh") = TranslateMesh_();
  out_list.sublist("Domain").set<int>("Spatial Dimension", dim_);
  out_list.sublist("Regions") = TranslateRegions_();
  out_list.sublist("Output") = TranslateOutput_();
  out_list.sublist("Solvers") = TranslateSolvers_();
  out_list.sublist("Preconditioners") = TranslatePreconditioners_();
  out_list.sublist("State") = TranslateState_();
  
  return out_list;
}


/* ******************************************************************
* Extract generic verbosity object for all sublists.
****************************************************************** */
Teuchos::ParameterList InputConverterU::GetVerbosity_()
{
  Teuchos::ParameterList vlist;

  DOMNodeList* node_list;
  DOMNode* node_attr;
  DOMNamedNodeMap* attr_map;
  char* text_content;
    
  // get execution contorls node
  node_list = doc_->getElementsByTagName(XMLString::transcode("execution_controls"));
  
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
              vlist.set<std::string>("verbosity", TrimString_(text_content));
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


/* ******************************************************************
* Empty
****************************************************************** */
Teuchos::Array<double> InputConverterU::MakeCoordinates_(char* char_array)
{
  Teuchos::Array<double> coords;
  char* tmp;
  tmp = strtok(char_array, "(, ");

  while (tmp != NULL) {
    std::string str(tmp);
    boost::algorithm::trim(str);
    coords.append(std::strtod(str.c_str(), NULL));
    tmp = strtok(NULL, ",");
  }

  return coords;
}

}  // namespace AmanziInput
}  // namespace Amanzi
