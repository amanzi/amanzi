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
#include <xercesc/util/PlatformUtils.hpp>
#include <xercesc/parsers/DOMLSParserImpl.hpp>
#include <xercesc/framework/StdOutFormatTarget.hpp>
#include <xercesc/util/OutOfMemoryException.hpp>

#include "Teuchos_StandardParameterEntryValidators.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

// Amanzi's
#include "ErrorHandler.hpp"
#include "InputConverter.hh"


namespace Amanzi {
namespace AmanziInput {

XERCES_CPP_NAMESPACE_USE

/* ******************************************************************
* Initialization of xercecs document.
****************************************************************** */
void InputConverter::Init(const std::string& xmlfilename)
{
  Teuchos::ParameterList out_list;
  
  XercesDOMParser *parser = new XercesDOMParser();
  parser->setExitOnFirstFatalError(true);
  parser->setValidationConstraintFatal(true);
  parser->setValidationScheme(XercesDOMParser::Val_Never);
  parser->setDoNamespaces(true);

  AmanziErrorHandler* errorHandler = new AmanziErrorHandler();
  parser->setErrorHandler(errorHandler);
  parser->useCachedGrammarInParse(true);
 
  bool errorsOccured = false;

  try {
    parser->parse(xmlfilename.c_str());
  }
  catch (const OutOfMemoryException& e) {
    std::cerr << "OutOfMemoryException" << std::endl;
    errorsOccured = true;
    Exceptions::amanzi_throw(Errors::Message("Ran out of memory while parsing the input file. Aborting."));
  }
  catch (...) {
    errorsOccured = true;
    Exceptions::amanzi_throw(Errors::Message("Errors occured while parsing the input file. Aborting."));
  }

  doc_ = parser->getDocument();

  // grab verbosity early
  verb_list_ = GetVerbosity_();
  vo_ = new VerboseObject("InputTranslator", verb_list_);
  Teuchos::OSTab tab = vo_->getOSTab();
  
  delete errorHandler;
}


/* ******************************************************************
* Extract generic verbosity object for all sublists.
****************************************************************** */
Teuchos::ParameterList InputConverter::GetVerbosity_()
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
Teuchos::Array<double> InputConverter::MakeCoordinates_(char* char_array)
{
  Teuchos::Array<double> coords;
  char* tmp;
  tmp = strtok(char_array, "(, ");

  while (tmp != NULL) {
    std::string str(tmp);
    boost::algorithm::trim(str);
    coords.append(atof(str.c_str()));
    tmp = strtok(NULL, ",");
  }

  return coords;
}


/* ******************************************************************
* Empty
****************************************************************** */
std::string InputConverter::TrimString_(char* tmp)
{
  std::string str(tmp);
  boost::algorithm::trim(str);
  return str;
}


/* *******************************************************************
* Generate unified error message for ill-formed element
******************************************************************* */
void InputConverter::ThrowErrorIllformed_(
    std::string section, std::string element_type, std::string ill_formed)
{
  Errors::Message msg;
  msg << "Amanzi::InputConverter: An error occurred during parsing " << section << "- ";
  msg << "  Missing or Ill-formed '" << element_type << "' for '" << ill_formed << "'. \n";
  msg << "  Please correct and try again \n";
  Exceptions::amanzi_throw(msg);
}


/* *****************************************************************************
* Generate unified error message for ill-formed element with options provided
***************************************************************************** */
void InputConverter::ThrowErrorIllformed_(
    std::string section, std::string element_type, std::string ill_formed, std::string options)
{
  Errors::Message msg;
  msg << "Amanzi::InputTranslator: An error occurred during parsing " << section << " - " ;
  msg << "  Missing or Ill-formed '" << element_type << "' for '" << ill_formed << "'. Valid options are: " << options << "\n" ;
  msg << "  Please correct and try again \n" ;
  Exceptions::amanzi_throw(msg);
}


/* *******************************************************************
* Generate unified error message for missing item
******************************************************************* */
void InputConverter::ThrowErrorMissattr_(
    std::string section, std::string att_elem_type, std::string missing, std::string elem_name)
{
  Errors::Message msg;
  msg << "Amanzi::InputConverter: An error occurred during parsing " << section << " - \n";
  msg << "  No " << att_elem_type << " " << missing << " found for " << elem_name << ". \n";
  msg << "  Please correct and try again \n";
  Exceptions::amanzi_throw(msg);
}

}  // namespace AmanziInput
}  // namespace Amanzi
