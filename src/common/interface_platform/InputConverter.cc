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

const Teuchos::Array<std::string> meshfileStrings = 
   Teuchos::tuple<std::string>("exodus ii", "exodus II", "Exodus II", "Exodus ii", "H5M", "h5m");

/* ******************************************************************
* Empty
****************************************************************** */
Teuchos::ParameterList InputConverter::Translate(const std::string& xmlfilename)
{
  Teuchos::ParameterList out_list;
  
  XMLPlatformUtils::Initialize();

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

  DOMDocument *doc = parser->getDocument();

  // grab verbosity early
  Teuchos::ParameterList verb = GetVerbosity_(doc);
  vo_ = new VerboseObject("InputTranslator", verb);
  Teuchos::OSTab tab = vo_->getOSTab();
  
  // translating...
  out_list.sublist("Mesh") = TranslateMesh_(doc);
  
  delete errorHandler;
  XMLPlatformUtils::Terminate();

  return out_list;
}


/* ******************************************************************
* Empty
****************************************************************** */
Teuchos::ParameterList InputConverter::TranslateMesh_(DOMDocument* doc)
{
  Teuchos::ParameterList out_list;

  bool generate = true;
  bool read = false;
  char *framework;
  Teuchos::ParameterList mesh_list;
  bool all_good = false;
  Errors::Message msg;
  std::stringstream helper;
    
  if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH) {
    *vo_->os() << "Translating mesh" << std::endl;
  }

  Teuchos::RCP<Teuchos::StringValidator> meshfileValidator = Teuchos::rcp(new Teuchos::StringValidator(meshfileStrings));
  
  // read in new stuff
  XMLCh* tag = XMLString::transcode("mesh");
  DOMNodeList* nodeList = doc->getElementsByTagName(tag);
  XMLString::release(&tag);

  // read the attribute to set the framework sublist
  const XMLSize_t nodeCount = nodeList->getLength() ;

  if (nodeList->getLength() > 0) {
    DOMNode* nodeMesh = nodeList->item(0);
    DOMElement* elementMesh = static_cast<DOMElement*>(nodeMesh);
    // if unstructure, look for framework attribute
    if(elementMesh->hasAttribute(XMLString::transcode("framework"))) {
      framework = XMLString::transcode(elementMesh->getAttribute(XMLString::transcode("framework")));
    } else { 
      msg << "Amanzi::InputConverter: An error occurred during parsing mesh - "
          << "framework was missing or ill-formed.\n" 
          << "Use default framework='mstk' if unsure. Please correct and try again.\n";
      Exceptions::amanzi_throw(msg);
    }

    // first figure out what the dimension is
    DOMNodeList* children = nodeMesh->getChildNodes();
    all_good = false;

    for (int i = 0; i < children->getLength(); i++) {
      DOMNode* inode = children->item(i) ;
      if (DOMNode::ELEMENT_NODE == inode->getNodeType()) {
	char* tagname = XMLString::transcode(inode->getNodeName());
	if (strcmp(tagname,"dimension") == 0) {
	  char* temp = XMLString::transcode(inode->getTextContent());
	  if (strlen(temp) > 0) {
	      dim_ = atoi(temp);
	      all_good = true;
	  }
          XMLString::release(&temp);
	}
        XMLString::release(&tagname);
      }
    }

    if (!all_good) {
      ThrowErrorIllformed_("mesh", "element", "dimension");
    }

    // now we can properly parse the generate/read list
    all_good = false;
    for (int i = 0; i < children->getLength(); i++) {
      DOMNode* inode = children->item(i) ;
      if (DOMNode::ELEMENT_NODE == inode->getNodeType()) {
	char* tagname = XMLString::transcode(inode->getNodeName());   

	if (strcmp(tagname,"generate")==0) {
	  all_good = true;
	  generate = true;
	  read = false;
	  DOMElement* elementGen = static_cast<DOMElement*>(inode);

	  // get Number of Cells
	  DOMNodeList* node_list = elementGen->getElementsByTagName( XMLString::transcode("number_of_cells"));
	  DOMNode* node = node_list->item(0);
	  DOMElement* element_node = static_cast<DOMElement*>(node);
	  DOMNamedNodeMap *attr_map = node->getAttributes();

	  Teuchos::Array<int> ncells; 
          DOMNode* node_attr;
          char* attr_name;
	  char* temp;

	  // make sure number of attributes equals dimension
	  if (attr_map->getLength() == dim_) {
	    // loop over attributes to get nx, ny, nz as needed
	    for (int j = 0; j < attr_map->getLength(); j++) {
	      node_attr = attr_map->item(j);
	      attr_name =XMLString::transcode(node_attr->getNodeName());

	      if (attr_name) {
	        temp = XMLString::transcode(node_attr->getNodeValue());
	        if (strlen(temp) > 0) {
	          ncells.append(atoi(temp));
		} else {
		  all_good = false;
	          helper << "number_of_cells "<<attr_name;
	        }
                XMLString::release(&temp);
	      } else {
		all_good = false;
	        helper << "number_of_cells "<<attr_name;
	      }
	      XMLString::release(&attr_name);
	    }
            mesh_list.set<Teuchos::Array<int> >("Number of Cells",ncells);
	  } else {
	    helper << "number_of_cells";
	    all_good = false;
	  }

	  // get Box - generalize
	  //char* char_array;
	  nodeList = elementGen->getElementsByTagName( XMLString::transcode("box"));
	  node = nodeList->item(0);
	  element_node = static_cast<DOMElement*>(node);
	  temp = XMLString::transcode(element_node->getAttribute( XMLString::transcode("low_coordinates")));
	  if (strlen(temp) > 0) {
	    // translate to array
	    Teuchos::Array<double> low = MakeCoordinates_(temp);
            mesh_list.set<Teuchos::Array<double> >("Domain Low Coordinate", low);
	    if (low.length() != dim_) {
	      helper << "low_coordinates";
	      all_good = false;
	    }
	  } else {
	    helper << "low_coordinates";
	    all_good = false;
	  }
	  XMLString::release(&temp);
	  temp = XMLString::transcode(element_node->getAttribute( XMLString::transcode("high_coordinates")));
	  if (strlen(temp) > 0) {
	    // translate to array
	    Teuchos::Array<double> high = MakeCoordinates_(temp);
            mesh_list.set<Teuchos::Array<double> >("Domain High Coordinate",high);
	    if (high.length() != dim_) {
	      helper << "high_coordinates";
	      all_good = false;
	    }
	  }
          else {
	    helper << "high_coordinates";
	    all_good = false;
	  }
	  XMLString::release(&temp);
	}

	else if (strcmp(tagname,"read")==0) {
	  read = true;
	  generate = false;
	  bool goodtype = false;
	  bool goodname = false;
	  DOMElement* elementRead = static_cast<DOMElement*>(inode);

	  char* value = XMLString::transcode(elementRead->getElementsByTagName(
              XMLString::transcode("format"))->item(0)->getTextContent());
          std::string format(TrimString_(value));

          if (boost::iequals(format, "exodus ii")) {
            mesh_list.set<std::string>("Format", "Exodus II");
            goodtype = true;
	  } else if (boost::iequals(format, "h5m")) {
            mesh_list.set<std::string>("Format", "H5M");
            goodtype = true;
          } else {
            mesh_list.set<std::string>("Format", format, "Format of meshfile", meshfileValidator);
          }
	  char* filename = XMLString::transcode(elementRead->getElementsByTagName(
				  XMLString::transcode("file"))->item(0)->getTextContent());
	  if (strlen(filename) > 0) {
            mesh_list.set<std::string>("File", TrimString_(filename));
	    goodname = true;
	  }
	  XMLString::release(&value);
	  XMLString::release(&filename);
	  if (goodtype && goodname) all_good = true;
	}

	//EIB - handles legacy files
	//TODO:: EIB - remove to conform to updated schema
	else if (strcmp(tagname,"file")==0) {
	  read = true;
	  generate = false;
	  char* filename = XMLString::transcode(inode->getTextContent());
	  if (strlen(filename) > 0) {
            mesh_list.set<std::string>("File", TrimString_(filename));
	    mesh_list.set<std::string>("Format","Exodus II");
	    all_good = true;
	    std::cout << "Amanzi::InputTranslator: Warning - " << std::endl;
	    std::cout << "    Please note - the XML Schema for specifing a mesh file to read has been updated." << std::endl;
	    std::cout << "    See the Amanzi User Guide for the latest information." << std::endl;
	    std::cout <<"     The legacy format is being handled for now.  Please update input files for future versions." << std::endl;
	  }
	  XMLString::release(&filename);
	}

      }
    }
    if (!all_good) {
      ThrowErrorIllformed_("mesh", helper.str(), "generate/read");
    }
    
    if (generate || read) {
      if (strcmp(framework,"mstk")==0) {
        out_list.sublist("Unstructured").sublist("Expert").set<std::string>("Framework","MSTK");
      } else if (strcmp(framework,"moab")==0) {
        out_list.sublist("Unstructured").sublist("Expert").set<std::string>("Framework","MOAB");
      } else if (strcmp(framework,"simple")==0) {
        out_list.sublist("Unstructured").sublist("Expert").set<std::string>("Framework","Simple");
      } else if (strcmp(framework,"stk::mesh")==0) {
        out_list.sublist("Unstructured").sublist("Expert").set<std::string>("Framework","stk::mesh");
      } else {
        //list.sublist("Unstructured").sublist("Expert").set<std::string>("Framework","MSTK");
        msg << "Amanzi::InputTranslator: An error occurred during parsing mesh - "
            << "unknown framework=" << framework << ". See the schema for acceptable types.\n"  
            << "Please correct and try again \n";
        Exceptions::amanzi_throw(msg);
      }
    }

    if (generate) {
      out_list.sublist("Unstructured").sublist("Generate Mesh").sublist("Uniform Structured") = mesh_list;
    } else if (read) {
      out_list.sublist("Unstructured").sublist("Read Mesh File") = mesh_list;
    } else {
      // bad mesh, again if validated shouldn't need this
      msg << "Amanzi::InputTranslator: An error occurred during parsing mesh.\n"
          << "Please correct and try again.\n";
      Exceptions::amanzi_throw(msg);
    }
  } else {
    // no mesh sub-elements
    ThrowErrorIllformed_("mesh", "element", "framework");
  }

  return out_list;
}


/* ******************************************************************
* Empty
****************************************************************** */
Teuchos::ParameterList InputConverter::GetVerbosity_(DOMDocument* doc)
{
  Teuchos::ParameterList vlist;

  DOMNodeList* node_list;
  DOMNode* node_attr;
  DOMNamedNodeMap* attr_map;
  char* text_content;
    
  // get execution contorls node
  node_list = doc->getElementsByTagName(XMLString::transcode("execution_controls"));
  
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
  msg << "Amanzi::InputConverter: An error occurred during parsing " << section << "- " ;
  msg << "  Missing or Ill-formed '" << element_type << "' for '" << ill_formed << "'. \n" ;
  msg << "  Please correct and try again \n" ;
  Exceptions::amanzi_throw(msg);
}

}  // namespace AmanziInput
}  // namespace Amanzi
